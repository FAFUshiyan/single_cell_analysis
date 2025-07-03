#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# annotate_clusters.R  (comma‑only version)
# ---------------------------------------------------------------------------
# Re‑annotate clusters in a Seurat object using an external CSV mapping file
# with two columns: <cluster>,<annotation> (comma‑separated, **no spaces** in
# the delimiter).  The script saves diagnostic plots and the updated object.
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
  library(future)
  library(tibble)
  library(readr)
})

plan(multisession, workers = 4)
options(future.globals.maxSize = +Inf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper: save_plot ----------------------------------------------------------
save_plot <- function(p, filepath, w = 14, h = 10) {
  message(" → Saving `", filepath, "`")
  ggsave(filepath, plot = p, width = w, height = h, units = "in")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main ----------------------------------------------------------------------
main <- function() {

  option_list <- list(
    make_option(c("-a", "--anno_file"), type = "character", metavar = "FILE",
                help = "CSV mapping file: <cluster>,<annotation> (comma‑separated only)."),
    make_option(c("-r", "--rds_file"),  type = "character", metavar = "FILE",
                help = "Path to clustered Seurat *.rds object."),
    make_option(c("-o", "--out_dir"),   type = "character", metavar = "DIR",
                help = "Output directory (plots/ and rds/ will be created inside)."),
    make_option(c("-p", "--prefix"),    type = "character", default = "Sample",
                help = "Prefix for output files [default %default]."),
    make_option(c("-g", "--genes"),     type = "character", default = NULL,
                help = "Comma‑separated list of marker genes for expression plots (optional)."),
    make_option(c("-s", "--split_by"),  type = "character", default = "group",
                help = "Metadata column for split.by in gene plots [default %default].")
  )

  opts <- parse_args(OptionParser(option_list = option_list))

  # ---- sanity checks -------------------------------------------------------
  for (fld in c("anno_file", "rds_file", "out_dir")) {
    if (is.null(opts[[fld]])) stop("Missing required option --", fld)
  }
  if (!file.exists(opts$anno_file)) stop("Annotation file not found: ", opts$anno_file)
  if (!file.exists(opts$rds_file))  stop("RDS file not found: ", opts$rds_file)

  dir_create(opts$out_dir)
  plots_dir   <- path(opts$out_dir, "plots");   dir_create(plots_dir)
  objects_dir <- path(opts$out_dir, "rds");     dir_create(objects_dir)

  # ---- read annotation table (comma‑only) ---------------------------------
  if (!grepl(",", read_lines(opts$anno_file, n_max = 1))) {
    stop("The annotation file appears not to be comma‑separated. Please format it as <cluster>,<annotation> on each line.")
  }
  anno_tbl <- read_delim(opts$anno_file, delim = ",",
                         col_names = c("cluster", "anno"),
                         trim_ws = TRUE, show_col_types = FALSE)
  if (!is.numeric(anno_tbl$cluster)) anno_tbl$cluster <- as.integer(anno_tbl$cluster)
  new.ids <- setNames(anno_tbl$anno, anno_tbl$cluster)

  # ---- load Seurat object --------------------------------------------------
  obj <- readRDS(opts$rds_file)
  obj$cluster_orig <- Idents(obj)

  # ---- re‑annotate ---------------------------------------------------------
  Idents(obj) <- factor(Idents(obj), levels = names(new.ids), labels = new.ids)
  obj$cell_type <- Idents(obj)
  message("Annotation done. Distribution:")
  print(table(obj$cell_type))

  # ---- UMAP plot -----------------------------------------------------------
  p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
            ggtitle(paste0(opts$prefix, " UMAP (annotated)"))
  save_plot(p_umap, path(plots_dir, paste0(opts$prefix, "_umap.annotated.pdf")), w = 7, h = 5)

  # ---- optional gene plots -------------------------------------------------
  if (!is.null(opts$genes)) {
    genes <- stringr::str_split(opts$genes, ",", simplify = TRUE) |>
             as.vector() |>
             stringr::str_trim()
    walk(genes, function(g) {
      if (!g %in% rownames(obj)) {
        warning("Gene `", g, "` not present in the object; skipping.")
        return(NULL)
      }
      p_vln <- VlnPlot(obj, features = g, group.by = "cell_type", split.by = opts$split_by) +
               ggtitle(paste0(g, " expression"))
      save_plot(p_vln, path(plots_dir, paste0(opts$prefix, "_", g, ".vln.pdf")), w = 7, h = 5)

      p_feat <- FeaturePlot(obj, features = g, split.by = opts$split_by, reduction = "umap") +
                ggtitle(paste0(g, " expression (split)"))
      save_plot(p_feat, path(plots_dir, paste0(opts$prefix, "_", g, ".feature.pdf")), w = 14, h = 5)
    })
  }

  # ---- save RDS ------------------------------------------------------------
  saveRDS(obj, path(objects_dir, paste0(opts$prefix, ".annotated.rds")))
  message("All done! Results written to ", opts$out_dir)
}

# run -----------------------------------------------------------------------
if (sys.nframe() == 0) main()
