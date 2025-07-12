#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# sc_subset_recluster.R
# ---------------------------------------------------------------------------
# Extract a subset of clusters from an integrated Seurat object, then
# re‑normalise (SCTransform), batch‑correct (Harmony), and optionally cluster
# and generate marker DotPlots.
# ---------------------------------------------------------------------------
# CLI EXAMPLE ---------------------------------------------------------------
# Rscript sc_subset_recluster.R \
#     --rds_file  PBMC_integrated.rds \
#     --clusters  6,17,19,29 \
#     --subset    Bcells \
#     --out_dir   results/Bcells \
#     --npc       30 \
#     --resolution 1.2 \
#     --markers_tsv markers.tsv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(fs)
  library(stringr)
})

plan(multisession, workers = 4)
options(future.globals.maxSize = +Inf, future.seed = TRUE)

# ---------------------------------------------------------------------------
# Helper --------------------------------------------------------------------

save_plot <- function(p, path, w = 8, h = 6) {
  message(" → Saving `", path, "`")
  ggsave(path, plot = p, width = w, height = h, units = "in")
}

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# CLI opts ------------------------------------------------------------------

opt_list <- list(
  make_option("--rds_file",    type = "character", help = "Integrated Seurat *.rds* file"),
  make_option("--clusters",    type = "character", help = "Comma‑separated cluster IDs to extract"),
  make_option("--subset",      type = "character", default = "subset", help = "Name/tag for this subset"),
  make_option("--out_dir",     type = "character", help = "Output directory for subset results"),
  make_option("--npc",         type = "integer",  default = 20,    help = "Number of PCs for reclustering"),
  make_option("--resolution",  type = "double",   default = NA,   help = "Resolution for FindClusters; NA = skip"),
  make_option("--markers_tsv", type = "character", default = NULL, help = "Marker TSV (celltype & gene) for DotPlot")
)
opts <- parse_args(OptionParser(option_list = opt_list))

for (req in c("rds_file", "clusters", "out_dir")) if (is.null(opts[[req]])) stop("--", req, " is required")
if (!file.exists(opts$rds_file)) stop("rds_file not found: ", opts$rds_file)

dir_create(opts$out_dir, recurse = TRUE)
plots_dir   <- path(opts$out_dir, "plots");   dir_create(plots_dir)
objects_dir <- path(opts$out_dir, "rds");     dir_create(objects_dir)

# ---------------------------------------------------------------------------
# Load object & subset ------------------------------------------------------

obj <- readRDS(opts$rds_file)
clust_vec <- as.integer(str_split(opts$clusters, ",", simplify=TRUE))
msg("Subsetting clusters: ", paste(clust_vec, collapse=","))
sub <- subset(obj, idents = clust_vec)

# ---------------------------------------------------------------------------
# Re‑normalise & integrate ---------------------------------------------------

assay_nm <- paste0("SCT_", opts$subset)
sub <- SCTransform(sub, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE,
                   new.assay.name = assay_nm, return.only.var.genes = FALSE)
sub <- RunPCA(sub, npcs = opts$npc, assay = assay_nm, verbose = FALSE)
sub <- RunHarmony(sub, group.by.vars = "batch", assay.use = assay_nm,
                  reduction.use = "pca", dims.use = 1:opts$npc)
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:opts$npc)
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:opts$npc)

if (!is.na(opts$resolution)) {
  msg("FindClusters at resolution ", opts$resolution)
  sub <- FindClusters(sub, resolution = opts$resolution)
}

# ---------------------------------------------------------------------------
# Save objects & plots -------------------------------------------------------

saveRDS(sub, path(objects_dir, paste0(opts$subset, "_recluster.rds")))

p_umap <- DimPlot(sub, reduction = "umap", label = !is.na(opts$resolution)) +
          ggtitle(paste(opts$subset, "subset"))
save_plot(p_umap, path(plots_dir, paste0(opts$subset, "_umap.pdf")), w = 7, h = 5)

if (!is.null(opts$markers_tsv)) {
  mk_tbl <- read.table(opts$markers_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  genes  <- unique(mk_tbl$gene)
  group_by <- if (!is.na(opts$resolution)) "seurat_clusters" else "batch"
  p_dot <- DotPlot(sub, features = genes, assay = assay_nm, group.by = group_by, dot.scale = 6,
                   cols = c("grey90", "royalblue")) + RotatedAxis() +
           labs(title = paste(opts$subset, "marker genes"), x = group_by, y = NULL)
  save_plot(p_dot, path(plots_dir, paste0(opts$subset, "_markers.pdf")), w = 14, h = 6)
}

msg("Subset pipeline finished ✓")