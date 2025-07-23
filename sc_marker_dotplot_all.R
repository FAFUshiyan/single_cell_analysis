#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# sc_marker_dotplot_all.R
#   * 一次读取 RDS，可针对 **多组 resolution** 依次聚类 → 生成 UMAP / DotPlot / (可选) FindAllMarkers
#   * 默认自动优先使用 `harmony` 降维；亦可 `--reduction pca` 等手动指定
#   * 支持 SCT + Harmony 管线，确保仅调分辨率不重算 UMAP
# ---------------------------------------------------------------------------
# 调用示例：
#   Rscript sc_marker_dotplot_all.R \
#       --rds_file integrated.rds \
#       --markers_tsv markers.tsv \
#       --out_dir results \
#       --prefix PBMC \
#       --resolution 0.4,0.8,1.2   # 多分辨率逗号分隔
#       --run_markers              # 如需差异表达
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(fs)
})

save_plot <- function(p, path, w = 8, h = 6) {
  message(" → Saving `", path, "`")
  ggsave(path, plot = p, width = w, height = h, units = "in")
}

opt_list <- list(
  make_option("--rds_file",    type = "character", help = "Seurat *.rds* 文件"),
  make_option("--markers_tsv", type = "character", help = "两列 TSV (celltype gene)"),
  make_option("--out_dir",     type = "character", default = "dotplots", help = "输出目录"),
  make_option("--prefix",      type = "character", default = "", help = "输出文件前缀 (可选)"),
  make_option("--group_by",    type = "character", default = "seurat_clusters", help = "DotPlot 分组元数据列"),
  make_option("--assay",       type = "character", default = NULL, help = "使用的 assay；缺省用默认"),
  make_option("--dot_scale",   type = "double",   default = 6, help = "DotPlot dot.scale 参数"),
  make_option("--palette",     type = "character", default = "grey90,royalblue", help = "两端色 (逗号分隔)"),
  make_option("--npc",         type = "integer",  default = 30, help = "使用的 PCA / harmony 维度数 (dims = 1:npc)"),
  make_option("--resolution",  type = "character", default = "0.5", help = "分辨率，可逗号分隔多值，如 0.4,0.8"),
  make_option("--reduction",   type = "character", default = "auto", help = "降维矩阵 (auto|harmony|pca|...)"),
  make_option("--run_markers", action = "store_true", default = FALSE, help = "执行 FindAllMarkers 并输出结果")
)
opts <- parse_args(OptionParser(option_list = opt_list))

## --- 基本检查 -------------------------------------------------------------
for (req in c("rds_file", "markers_tsv"))
  if (is.null(opts[[req]])) stop("--", req, " is required")
if (!file.exists(opts$rds_file))    stop("rds_file not found: ", opts$rds_file)
if (!file.exists(opts$markers_tsv)) stop("markers_tsv not found: ", opts$markers_tsv)
dir_create(opts$out_dir, recurse = TRUE)

## --- 解析 resolution 列表 --------------------------------------------------
res_vec <- as.numeric(strsplit(opts$resolution, ",")[[1]])
if (any(is.na(res_vec))) stop("--resolution 需为数字或逗号分隔的数字，例如 0.4,0.8")

## --- 读对象 & assay -------------------------------------------------------
seu <- readRDS(opts$rds_file)
assay_use <- opts$assay %||% DefaultAssay(seu)
if (!assay_use %in% Assays(seu)) stop("Assay ", assay_use, " 不存在于对象中")

## --- 选择降维矩阵 ---------------------------------------------------------
reduction_use <- if (opts$reduction != "auto") {
  opts$reduction
} else if ("harmony" %in% Reductions(seu)) {
  "harmony"
} else {
  "pca"
}
if (!reduction_use %in% Reductions(seu)) stop("Reduction ", reduction_use, " 不存在于对象中")
message("Using reduction = ", reduction_use)

## --- FindNeighbors 一次即可 ----------------------------------------------
dims_use <- 1:opts$npc
seu <- FindNeighbors(seu, reduction = reduction_use, dims = dims_use)

## --- Marker 基因准备 ------------------------------------------------------
mk_tbl <- read.table(opts$markers_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!"gene" %in% colnames(mk_tbl)) stop("markers_tsv 必须含 'gene' 列")

genes_all <- unique(mk_tbl$gene)
pal <- strsplit(opts$palette, ",")[[1]][1:2]

## --- 主循环：多分辨率 ------------------------------------------------------
for (res in res_vec) {
  message("\n==========  Processing resolution = ", res, "  ==========")
  seu <- FindClusters(seu, resolution = res)

  # 重新计算 UMAP 仅在首次循环或 npc 变化才发生 — 已在 Read 步骤保持不变
  if (!("umap" %in% Reductions(seu))) {
    message("Running RunUMAP dims=1:", opts$npc)
    seu <- RunUMAP(seu, reduction = reduction_use, dims = dims_use)
  }

  prefix <- ifelse(nzchar(opts$prefix), paste0(opts$prefix, "_"), "")
  res_tag <- sprintf("n%d_r%.2f", opts$npc, res)

  ## --- UMAP 图 -----------------------------------------------------------
  plot_title <- sprintf("UMAP (%s dims=1-%d, res=%.2f)", reduction_use, opts$npc, res)
  p_umap <- DimPlot(seu, label = TRUE, reduction = "umap") + NoLegend() + ggtitle(plot_title)
  umap_path <- file.path(opts$out_dir, paste0(prefix, "umap_", reduction_use, "_", res_tag, ".pdf"))
  save_plot(p_umap, umap_path, w = 8, h = 6)

  ## --- DotPlot ----------------------------------------------------------
  genes_present <- intersect(genes_all, rownames(seu))
  w_dot <- max(14, 14 + (length(genes_present) - 5) * 0.2)
  dot_title <- sprintf("All marker genes (%s dims=1-%d, res=%.2f)", reduction_use, opts$npc, res)
  p_dot <- DotPlot(seu, features = genes_present, group.by = opts$group_by, assay = assay_use, dot.scale = opts$dot_scale, cols = pal) +
           RotatedAxis() + labs(title = dot_title, x = opts$group_by, y = NULL)
  dot_path <- file.path(opts$out_dir, paste0(prefix, "dotplot_", reduction_use, "_", res_tag, ".pdf"))
  save_plot(p_dot, dot_path, w = w_dot, h = 14)

  ## --- 可选 FindAllMarkers ---------------------------------------------
  if (opts$run_markers) {
    message("Running FindAllMarkers (res=", res, ") ...")
    Idents(seu) <- opts$group_by
    if (assay_use == "SCT") seu <- PrepSCTFindMarkers(seu, assay = assay_use)
    all_markers <- FindAllMarkers(seu, assay = assay_use, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", max.cells.per.ident = 3000)
    marker_path <- file.path(opts$out_dir, paste0(prefix, "FindAllMarkers_", reduction_use, "_", res_tag, ".tsv"))
    write.table(all_markers, file = marker_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Markers saved: ", marker_path)
  }

  ## --- 保存中间 Seurat 对象 --------------------------------------------
  rds_path <- file.path(opts$out_dir, paste0(prefix, "seu_", reduction_use, "_", res_tag, ".rds"))
  saveRDS(seu, rds_path)
  message("Saved Seurat object: ", rds_path)
}

message("All resolutions processed ✓")
