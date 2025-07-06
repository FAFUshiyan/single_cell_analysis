#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# sc_preprocess_pipeline.R 
# ---------------------------------------------------------------------------
# Workflow:
#   1) Read each 10x dataset, basic QC
#   2) Merge all samples
#   3) SCTransform normalisation + Harmony batch correction (user‑settable npc)
#   4) *Optional* FindClusters on the integrated object (user‑settable resolution)
#   5) *Optional* DotPlot of marker genes (TSV list)
# ---------------------------------------------------------------------------
# CLI EXAMPLE ---------------------------------------------------------------
# Rscript sc_preprocess_pipeline.R \
#   --meta_csv   samples_meta.csv \
#   --out_dir    results/PBMC \
#   --prefix     PBMC \
#   --npc        40 \
#   --resolution 0.8 \
#   --markers_tsv markers.tsv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(purrr)
  library(dplyr)
  library(fs)
  library(future.apply)
  library(tibble)
  library(stringr)
})

plan(multisession, workers = 4)
options(future.globals.maxSize = +Inf, future.seed = TRUE)

# ---------------------------------------------------------------------------
# Utility wrappers -----------------------------------------------------------

save_plot <- function(p, path, w = 8, h = 6) {
  message(" → Saving `", path, "`")
  ggsave(path, plot = p, width = w, height = h, units = "in")
}

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# QC helpers -----------------------------------------------------------------

auto_cut <- function(v, lower.bound = 0, upper.cap = Inf) {
  med <- median(v); mad3 <- 3 * mad(v)
  c(low = max(lower.bound, med - mad3), up = min(upper.cap, med + mad3))
}

plot_qc <- function(sample_list, out_file) {
  qc_df <- map_dfr(sample_list, ~ tibble(
    sample   = unique(.x$orig.ident),
    nFeature = .x$nFeature_RNA,
    pctMT    = .x$percent.mt))

  p1 <- ggplot(qc_df, aes(sample, nFeature)) +
          geom_violin(fill = "#6BAED6", scale = "width") +
          geom_boxplot(width = .12, outlier.shape = NA) +
          labs(y = "nFeature_RNA", x = NULL) +
          theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2 <- ggplot(qc_df, aes(sample, pctMT)) +
          geom_violin(fill = "#FB6A4A", scale = "width") +
          geom_boxplot(width = .12, outlier.shape = NA) +
          labs(y = "percent.mt (%)", x = NULL) +
          theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(p1 | p2, out_file, w = 9, h = 5)
}

# ---------------------------------------------------------------------------
# Core functions -------------------------------------------------------------

create_obj <- function(row, gene_cap = 15000, mt_cap = 20) {
  mtx <- Read10X(row$data_dir)
  obj <- CreateSeuratObject(counts = mtx, project = row$sample,
                            min.cells = 3, min.features = 200)
  obj$batch <- row$batch; obj$group <- row$group
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  lim <- list(g  = auto_cut(obj$nFeature_RNA, lower.bound = 200, upper.cap = gene_cap),
              mt = auto_cut(obj$percent.mt,   lower.bound = 0,   upper_cap = mt_cap))
  keep <- with(lim, obj$nFeature_RNA > g["low"] & obj$nFeature_RNA < g["up"] & obj$percent.mt < mt["up"])
  obj <- obj[, keep]
  msg(row$sample, " kept ", sum(keep), " / ", length(keep), " cells after QC")
  obj
}

run_sct_harmony <- function(obj, npc = 50, harmony_dims = NULL, slot_name = "SCT") {
  harmony_dims <- harmony_dims %||% 1:npc
  obj <- SCTransform(obj, vars.to.regress ="percent.mt", vst.flavor = "v2",
                     verbose = FALSE, new.assay.name = slot_name)
  obj <- RunPCA(obj, npcs = npc, assay = slot_name, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "batch", assay.use = slot_name,
                    reduction.use = "pca", dims.use = harmony_dims,
                    plot_convergence = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = harmony_dims)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = harmony_dims)
  obj
}

plot_markers <- function(obj, markers_tsv, assay = "SCT", out_path, group.by = "seurat_clusters", title = NULL) {
  mk_tbl <- read.table(markers_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  genes  <- unique(mk_tbl$gene)
  p <- DotPlot(obj, features = genes, assay = assay, group.by = group.by, dot.scale = 6,
               cols = c("grey90", "royalblue")) + RotatedAxis() +
       labs(title = title %||% "Marker genes", x = group.by, y = NULL)
  save_plot(p, out_path, w = 14, h = 6)
}

# ---------------------------------------------------------------------------
# CLI -----------------------------------------------------------------------

op <- list(
  make_option("--meta_csv",    type = "character", help = "CSV/TSV: sample,data_dir,batch,group"),
  make_option("--out_dir",     type = "character", help = "Root output directory"),
  make_option("--prefix",      type = "character", default = "Sample", help = "File prefix"),
  make_option("--npc",         type = "integer",  default = 50, help = "Number of PCs for integration"),
  make_option("--resolution",  type = "double",   default = NA, help = "Resolution for FindClusters; NA = skip"),
  make_option("--markers_tsv", type = "character", default = NULL, help = "TSV (celltype & gene) for DotPlot"),
  make_option("--gene_cap",    type = "integer",  default = 15000, help = "nFeature upper cap"),
  make_option("--mt_cap",      type = "integer",  default = 20,    help = "percent.mt upper cap")
)
opts <- parse_args(OptionParser(option_list = op))

for (req in c("meta_csv", "out_dir")) if (is.null(opts[[req]])) stop("--", req, " is required")
if (!file.exists(opts$meta_csv)) stop("meta_csv not found: ", opts$meta_csv)

dir_create(opts$out_dir, recurse = TRUE)
plots_dir   <- path(opts$out_dir, "plots");   dir_create(plots_dir)
objects_dir <- path(opts$out_dir, "rds");     dir_create(objects_dir)

# ---------------------------------------------------------------------------
# Step 1 – Load & QC ---------------------------------------------------------

meta_tbl <- read.delim(opts$meta_csv, sep = ifelse(grepl("\.tsv$", opts$meta_csv), "\t", ","),
                       header = TRUE, stringsAsFactors = FALSE)

sample_list <- map(seq_len(nrow(meta_tbl)), function(i) {
  row <- meta_tbl[i, ]; create_obj(row, gene_cap = opts$gene_cap, mt_cap = opts$mt_cap)
})
names(sample_list) <- meta_tbl$sample
plot_qc(sample_list, path(plots_dir, "QC_violin.pdf"))

# ---------------------------------------------------------------------------
# Step 2 – Merge & integrate -------------------------------------------------

merged <- reduce(sample_list, merge, add.cell.ids = names(sample_list),
                 project = paste0(opts$prefix, "_merged"), merge.data = TRUE)
merged <- run_sct_harmony(merged, npc = opts$npc)

# ---------------------------------------------------------------------------
# Step 3 – optional clustering & markers ------------------------------------

if (!is.na(opts$resolution)) {
  msg("FindClusters at resolution ", opts$resolution)
  merged <- FindClusters(merged, resolution = opts$resolution)
}

saveRDS(merged, path(objects_dir, paste0(opts$prefix, "_integrated.rds")))

if (!is.null(opts$markers_tsv)) {
  group_by <- if (!is.na(opts$resolution)) "seurat_clusters" else "batch"
  plot_markers(merged, opts$markers_tsv, assay = "SCT",
               out_path = path(plots_dir, paste0(opts$prefix, "_markers.pdf")),
               group.by = group_by,
               title = paste(opts$prefix, "marker genes"))
}

msg("Pipeline finished ✓")
