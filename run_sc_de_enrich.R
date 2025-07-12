#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# run_sc_de_enrich.R
# ---------------------------------------------------------------------------
# End‑to‑end pipeline for:
#   1) single‑cell differential expression (Seurat + MAST)
#   2) GO (BP) over‑representation analysis (clusterProfiler)
#   3) GSEA (fgsea / MSigDB C5 BP)
# ---------------------------------------------------------------------------
# All parameters are supplied via the command line.  Typical usage:
#
#   Rscript run_sc_de_enrich.R \
#       --rds_file  NKTsub_recluster.anno.rds \
#       --out_dir   results/NKT \
#       --prefix    NKT \
#       --geneset_rds msigdbr.C5.rds
#
# Or, if DE has already been computed:
#
#   Rscript run_sc_de_enrich.R \
#       --de_csv    NKT_DE_case_vs_control_byCellType_SCTB.csv \
#       --out_dir   results/NKT \
#       --prefix    NKT
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)          # v5.1+
  library(future.apply)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(fs)
  library(glue)
  library(tibble)
})

plan(multisession, workers = 4)
options(future.globals.maxSize = +Inf, future.seed = TRUE)

# ---------------------------------------------------------------------------
# Helpers -------------------------------------------------------------------

save_plot <- function(p, path, w = 14, h = 10) {
  message(" → Saving plot `", path, "`")
  ggsave(path, plot = p, width = w, height = h, units = "in")
}

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ... , "\n"))

# ---------------------------------------------------------------------------
# CLI -----------------------------------------------------------------------

opt_list <- list(
  make_option(c("--rds_file"),  type = "character", help = "Path to reclustered Seurat RDS object."),
  make_option(c("--de_csv"),    type = "character", default = NULL, help = "Optional: existing DE CSV (skip DE step)."),
  make_option(c("--out_dir"),   type = "character", help = "Output root directory."),
  make_option(c("--prefix"),    type = "character", default = "Sample", help = "File prefix [default %default]."),
  make_option(c("--geneset_rds"), type = "character", default = "msigdbr.C5.rds", help = "RDS containing MSigDB gene sets [default %default]."),
  make_option(c("--min_genes"), type = "integer", default = 15, help = "Minimum genes for GO/GSEA [default %default]."),
  make_option(c("--min_cells"), type = "integer", default = 3, help = "Minimum cells per group to run DE [default %default].")
)
opts <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opts$out_dir)) stop("--out_dir is required")
if (is.null(opts$de_csv) && is.null(opts$rds_file)) stop("Provide either --rds_file or --de_csv (or both)")

dir_create(opts$out_dir)
plots_dir   <- path(opts$out_dir, "plots");   dir_create(plots_dir)
objects_dir <- path(opts$out_dir, "rds");     dir_create(objects_dir)
go_dir      <- path(opts$out_dir, "go");      dir_create(go_dir)

# ---------------------------------------------------------------------------
# Step 1 – Differential expression ------------------------------------------

de_csv_path <- if (!is.null(opts$de_csv)) opts$de_csv else path(objects_dir, paste0(opts$prefix, "_DE_case_vs_control_byCellType_SCTB.csv"))

if (!is.null(opts$de_csv) && file.exists(opts$de_csv)) {
  msg("Reading existing DE file …")
  de_all <- read.csv(opts$de_csv, stringsAsFactors = FALSE)
} else {
  if (is.null(opts$rds_file)) stop("Need --rds_file to compute DE")
  msg("Loading Seurat object …")
  NKT <- readRDS(opts$rds_file)

  # remove legacy SCT if present and set default assay
  if ("SCT" %in% Assays(NKT)) NKT[["SCT"]] <- NULL
  DefaultAssay(NKT) <- "SCT_B"
  NKT <- PrepSCTFindMarkers(NKT, assay = "SCT_B")

  Idents(NKT) <- "group"  # assumes meta column 'group' exists
  ctypes <- sort(unique(NKT$cell_type))

  msg("Running DE across ", length(ctypes), " cell types …")
  de_list <- future_lapply(ctypes, future.seed = TRUE, FUN = function(ct) {
    sub <- subset(NKT, cell_type == ct)
    if (length(unique(sub$group)) < 2 || min(table(sub$group)) < opts$min_cells) return(NULL)

    res <- FindMarkers(
      sub,
      ident.1 = "case", ident.2 = "control",
      assay  = "SCT_B", slot = "data",
      test.use = "MAST",
      latent.vars = "batch",
      logfc.threshold = 0,
      min.pct = 0.10,  #### 添加调整；
      recorrect_umi = FALSE)
    res$gene      <- rownames(res)
    res$cell_type <- ct
    res
  })
  de_all <- bind_rows(de_list)
  write.csv(de_all, de_csv_path, row.names = FALSE)
  msg("DE table saved → ", de_csv_path)
}

ctypes <- sort(unique(de_all$cell_type))

# ---------------------------------------------------------------------------
# Step 2 – GO and GSEA -------------------------------------------------------

msg("Loading gene sets for GSEA …")
if (!file.exists(opts$geneset_rds)) stop("Gene set RDS not found: ", opts$geneset_rds)
go_bp_df   <- readRDS(opts$geneset_rds)
go_bp_sets <- split(go_bp_df$gene_symbol, go_bp_df$gs_name)

for (ct in ctypes) {
  msg("[", ct, "] Enrichment")
  res_ct <- de_all %>% filter(cell_type == ct)
  if (nrow(res_ct) == 0) next

  # background & DEG sets ----------------------------------------------------
  bg_symbol  <- res_ct$gene
  deg_symbol <- res_ct %>% filter(p_val_adj < 0.05, abs(avg_log2FC) >= 0.25) %>% pull(gene)
  up_symbol  <- res_ct %>% filter(p_val_adj < 0.05,  avg_log2FC >=  0.25) %>% pull(gene)
  dn_symbol  <- res_ct %>% filter(p_val_adj < 0.05,  avg_log2FC <= -0.25) %>% pull(gene)

  gene_sets <- list(all = deg_symbol, up = up_symbol, down = dn_symbol)

  for (grp in names(gene_sets)) {
    genes <- gene_sets[[grp]]
    if (length(genes) < opts$min_genes) {
      msg("  • ", grp, " set (", length(genes), ") < min_genes → skip GO/GSEA")
      next
    }

    # ---- GO ORA -----------------------------------------------------------
    ego <- tryCatch({
      enrichGO(gene = genes,
               universe = bg_symbol,
               OrgDb = org.Hs.eg.db,
               keyType = "SYMBOL",
               ont = "BP",
               pAdjustMethod = "BH",
               qvalueCutoff = 1,
               readable = TRUE)
    }, error = function(e) NULL)

    if (!is.null(ego) && nrow(ego) > 0) {
      saveRDS(ego, file = path(go_dir, glue("{ct}_GO_{grp}.rds")))
      write.csv(as.data.frame(ego), path(go_dir, glue("{ct}_GO_{grp}.csv")), row.names = FALSE)

      p_bar <- barplot(ego, showCategory = 20, title = paste(ct, "GO BP", grp))
      save_plot(p_bar, path(go_dir, glue("{ct}_GO_{grp}_top20.pdf")), w = 14, h = 9)
    } else {
      msg("  • ", grp, " GO returned no terms → skipped plots")
    }
  }

  # ---- GSEA ---------------------------------------------------------------
  stat_vec <- sign(res_ct$avg_log2FC) * -log10(res_ct$p_val_adj)
  stat_vec[is.infinite(stat_vec)] <- max(stat_vec[is.finite(stat_vec)]) * 1.1
  names(stat_vec) <- res_ct$gene
  stat_vec <- sort(stat_vec, decreasing = TRUE)

  if (sum(!is.na(stat_vec)) < opts$min_genes) {
    msg("  • Not enough ranked genes for GSEA → skip")
    next
  }

  gsea_res <- fgsea(pathways = go_bp_sets,
                    stats    = stat_vec,
                    minSize  = opts$min_genes,
                    maxSize  = 500,
                    eps      = 0,
                    sampleSize = 500)
  saveRDS(gsea_res, path(go_dir, glue("{ct}_GSEA_GO_BP.rds")))
  out <- as.data.frame(gsea_res)
  out$leadingEdge <- vapply(out$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  write.csv(out[order(out$padj), ], path(go_dir, glue("{ct}_GSEA_GO_BP.csv")), row.names = FALSE)

  # plot NES top20
  if (nrow(out) > 0) {
    top20 <- head(out[order(out$padj), ], 20)
    p_gsea <- ggplot(top20, aes(x = reorder(pathway, NES), y = NES)) +
                geom_col() + coord_flip() +
                labs(title = paste(ct, "GSEA GO:BP (Top 20)"), x = NULL, y = "NES") +
                theme_bw()
    save_plot(p_gsea, path(go_dir, glue("{ct}_GSEA_GO_BP_top20.pdf")), w = 14, h = 9)
  }
}

msg("Pipeline complete ✔")
