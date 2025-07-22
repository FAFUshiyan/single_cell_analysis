#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# sc_marker_dotplot_all.R  ‑‑  全流程：重新聚类/UMAP + DotPlot + 可选 FindAllMarkers
# ---------------------------------------------------------------------------
# 调用示例：
#   Rscript sc_marker_dotplot_all.R \
#       --rds_file my_data.rds \
#       --markers_tsv marker.tsv \
#       --out_dir results \
#       --prefix PBMC \
#       --npc 40 \
#       --resolution 0.8 \
#       --run_markers   # ← 若要计算 FindAllMarkers
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
  make_option("--npc",         type = "integer",  default = 30, help = "使用的 PCA 维度数 (dims = 1:npc)"),
  make_option("--resolution",  type = "double",   default = 0.5, help = "FindClusters 分辨率"),
  make_option("--run_markers", action = "store_true", default = FALSE, help = "执行 FindAllMarkers 并输出结果")
)
opts <- parse_args(OptionParser(option_list = opt_list))

## --- 基本检查 -------------------------------------------------------------
for (req in c("rds_file", "markers_tsv"))
  if (is.null(opts[[req]])) stop("--", req, " is required")
if (!file.exists(opts$rds_file))    stop("rds_file not found: ", opts$rds_file)
if (!file.exists(opts$markers_tsv)) stop("markers_tsv not found: ", opts$markers_tsv)
dir_create(opts$out_dir, recurse = TRUE)

## --- 读对象 & assay -------------------------------------------------------
seu <- readRDS(opts$rds_file)
assay_use <- opts$assay %||% DefaultAssay(seu)
if (!assay_use %in% Assays(seu))
  stop("Assay ", assay_use, " 不存在于对象中")

## --- 重新聚类 & UMAP ------------------------------------------------------
dims_use <- 1:opts$npc
message("Running clustering / UMAP with dims = 1:", opts$npc, ", resolution = ", opts$resolution)

seu <- FindNeighbors(seu, dims = dims_use)
seu <- FindClusters(seu, resolution = opts$resolution)
seu <- RunUMAP(seu, dims = dims_use)

## --- 保存 UMAP 图 ---------------------------------------------------------
prefix <- ifelse(nzchar(opts$prefix), paste0(opts$prefix, "_"), "")
p_umap <- DimPlot(seu, label = TRUE, reduction = "umap") + NoLegend() +
  labs(title = paste0("UMAP clusters (dims=1-", opts$npc, ")"), subtitle = paste0("resolution = ", opts$resolution))

umap_path <- file.path(opts$out_dir, paste0(prefix, "umap_n", opts$npc, "_r", opts$resolution, ".pdf"))
save_plot(p_umap, umap_path, w = 8, h = 6)

## --- DotPlot ----------------------------------------------------------------
mk_tbl <- read.table(opts$markers_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!"gene" %in% colnames(mk_tbl)) stop("markers_tsv 必须含 'gene' 列")

genes <- unique(mk_tbl$gene)
genes_present  <- genes[genes %in% rownames(seu)]
genes_missing  <- setdiff(genes, genes_present)
if (length(genes_missing) > 0) message("ℹ 缺失基因: ", paste(genes_missing, collapse = ", "))

pal <- strsplit(opts$palette, ",")[[1]][1:2]
w_dot <- max(14, 14 + (length(genes_present) - 5) * 0.2)

p_dot <- DotPlot(seu, features = genes_present, group.by = opts$group_by, assay = assay_use, dot.scale = opts$dot_scale, cols = pal) +
  RotatedAxis() +
  labs(title = paste0("All marker genes (dims=1-", opts$npc, ", res=", opts$resolution, ")"), x = opts$group_by, y = NULL)

dot_path <- file.path(opts$out_dir, paste0(prefix, "dotplot_n", opts$npc, "_r", opts$resolution, ".pdf"))
save_plot(p_dot, dot_path, w = w_dot, h = 14)

## --- 可选 FindAllMarkers ---------------------------------------------------
if (opts$run_markers) {
  message("Running FindAllMarkers ...")
  Idents(seu) <- opts$group_by  # 使用 group_by 列
  if (assay_use == "SCT") seu <- PrepSCTFindMarkers(seu, assay = assay_use)

  all_markers <- FindAllMarkers(
    object          = seu,
    assay           = assay_use,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    test.use        = "wilcox",
    max.cells.per.ident = 3000
  )

  marker_path <- file.path(opts$out_dir, paste0(prefix, "FindAllMarkers_n", opts$npc, "_r", opts$resolution, ".tsv"))
  message(" → Saving marker table: ", marker_path)
  write.table(all_markers, file = marker_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("FindAllMarkers finished ✓")
} else {
  message("Skip FindAllMarkers (use --run_markers to enable)")
}

## --- 保存 Seurat 对象 ------------------------------------------------------
rds_path <- file.path(opts$out_dir, paste0(prefix, "seu_n", opts$npc, "_r", opts$resolution, ".rds"))
message(" → Saving updated Seurat object: ", rds_path)
saveRDS(seu, rds_path)

message("All tasks finished ✓")
