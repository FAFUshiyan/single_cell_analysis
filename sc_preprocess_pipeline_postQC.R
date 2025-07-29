#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# sc_preprocess_pipeline_postQC.R
#   * 读取每个样本已做 QC / 去双细胞的 Seurat 对象 (.rds)
#   * SCTransform 归一化 + Harmony 批次校正
#   * 可选聚类、Marker DotPlot
#   * **新增**: 生成 DimPlot (UMAP) PDF，并在标题中标注 npc、harmony_dims、resolution 参数
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

msg  <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# 读取已 QC 的对象 ------------------------------------------------------------

load_obj <- function(row) {
  obj <- readRDS(row$rds_file)
  obj$orig.ident <- row$sample   # 补充元数据
  obj$batch      <- row$batch
  obj$group      <- row$group
  msg("Loaded ", row$sample, " with ", ncol(obj), " cells")
  obj
}

# ---------------------------------------------------------------------------
# Harmony + SCT + 后续分析 ----------------------------------------------------

run_sct_harmony <- function(obj, npc = 50, harmony_dims = NULL, slot_name = "SCT") {
  harmony_dims <- harmony_dims %||% 1:npc
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", vst.flavor = "v2",
                     verbose = FALSE, new.assay.name = slot_name)
  obj <- RunPCA(obj, npcs = npc, assay = slot_name, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "batch", assay.use = slot_name,
                    reduction.use = "pca", dims.use = harmony_dims,
                    plot_convergence = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = harmony_dims)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = harmony_dims)
  obj
}

# ---------------------------------------------------------------------------
# DotPlot helper -------------------------------------------------------------

plot_markers <- function(obj, markers_tsv, assay = "SCT", out_path,
                         group.by = "seurat_clusters", title = NULL) {
  mk_tbl <- read.table(markers_tsv, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
  genes  <- unique(mk_tbl$gene)
  p <- DotPlot(obj, features = genes, assay = assay, group.by = group.by,
               dot.scale = 6, cols = c("grey90", "royalblue")) +
       RotatedAxis() +
       labs(title = title %||% "Marker genes", x = group.by, y = NULL)
  ggsave(out_path, plot = p, width = 16, height = 12, units = "in")
}

# ---------------------------------------------------------------------------
# CLI 设置 --------------------------------------------------------------------

opt_list <- list(
  make_option("--meta_csv",    type = "character",
              help = "CSV/TSV: sample,rds_file,batch,group (必填)"),
  make_option("--out_dir",     type = "character",
              help = "输出根目录 (必填)"),
  make_option("--prefix",      type = "character", default = "Sample",
              help = "文件名前缀"),
  make_option("--npc",         type = "integer",  default = 50,
              help = "PCA 维度数"),
  make_option("--resolution",  type = "double",   default = NA,
              help = "FindClusters 分辨率；NA = 跳过"),
  make_option("--markers_tsv", type = "character", default = NULL,
              help = "TSV: gene 列；用于 DotPlot")
)
opts <- parse_args(OptionParser(option_list = opt_list))
for (req in c("meta_csv", "out_dir"))
  if (is.null(opts[[req]])) stop("--", req, " is required")
if (!file.exists(opts$meta_csv)) stop("meta_csv not found: ", opts$meta_csv)

dir_create(opts$out_dir, recurse = TRUE)
plots_dir   <- path(opts$out_dir, "plots");   dir_create(plots_dir)
objects_dir <- path(opts$out_dir, "rds");     dir_create(objects_dir)

# ---------------------------------------------------------------------------
# Step 1 – 读取所有样本 --------------------------------------------------------

meta_tbl <- read.delim(opts$meta_csv,
                       sep    = ifelse(grepl("\\.tsv$", opts$meta_csv), "\t", ","),
                       header = TRUE, stringsAsFactors = FALSE)

sample_list <- map(seq_len(nrow(meta_tbl)), ~ load_obj(meta_tbl[., ]))
names(sample_list) <- meta_tbl$sample

# 如果只有一个样本，直接使用；否则 merge（需 add.cell.ids）
if (length(sample_list) == 1) {
  merged <- sample_list[[1]]
  msg("Single sample detected → skip merge() step")
} else {
  merged <- merge(x = sample_list[[1]],
                  y = sample_list[-1],
                  add.cell.ids = names(sample_list),
                  project = paste0(opts$prefix, "_merged"),
                  merge.data = TRUE)
}

# ---------------------------------------------------------------------------
# Step 2 – Harmony 整合 ------------------------------------------------------

merged <- run_sct_harmony(merged, npc = opts$npc)

# ---------------------------------------------------------------------------
# Step 3 – 可选聚类 & DotPlot -------------------------------------------------

if (!is.na(opts$resolution)) {
  msg("FindClusters at resolution ", opts$resolution)
  merged <- FindClusters(merged, resolution = opts$resolution)
}

saveRDS(merged, path(objects_dir, paste0(opts$prefix, "_integrated.rds")))

# 决定用于分组的列
cluster_group <- if (!is.na(opts$resolution)) "seurat_clusters" else "batch"

# DotPlot (可选)
if (!is.null(opts$markers_tsv)) {
  plot_markers(merged, opts$markers_tsv, assay = "SCT",
               out_path = path(plots_dir, paste0(opts$prefix, "_markers.pdf")),
               group.by = cluster_group,
               title = paste(opts$prefix, "marker genes"))
}

# ---------------------------------------------------------------------------
# Step 4 – 新增 DimPlot -------------------------------------------------------

plot_title <- sprintf("UMAP (npc=%d, harmony_dims=1:%d, res=%s)",
                      opts$npc, opts$npc,
                      ifelse(is.na(opts$resolution), "NA", opts$resolution))

p_dim <- DimPlot(merged, reduction = "umap", group.by = cluster_group,
                 label = TRUE, repel = TRUE) + ggtitle(plot_title)
p_dim_b <- DimPlot(merged, reduction = "umap", group.by = "batch",
                 label = TRUE, repel = TRUE) + ggtitle(plot_title)
ggsave(path(plots_dir, paste0(opts$prefix, "_UMAP_dimplot.cluster.pdf")),
       plot = p_dim, width = 16, height = 12, units = "in")
ggsave(path(plots_dir, paste0(opts$prefix, "_UMAP_dimplot.batch.pdf")),plot = p_dim_b, width = 16, height = 12, units = "in")
msg("Pipeline finished ✓")

