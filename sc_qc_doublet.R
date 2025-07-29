#!/usr/bin/env Rscript
# sc_qc_doublet.R
# 可配置的 Seurat MAD‑QC + DoubletFinder 流水线

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(purrr)
  library(dplyr)
  library(fs)
  library(future.apply)
  library(stringr)
  library(fields)
})


plan(multisession, workers = 4)
options(future.globals.maxSize = +Inf, future.seed = TRUE)

# ---------- 1. 命令行参数 ----------
opt_list <- list(
  make_option(c("-i", "--in_dir"),        type="character", help="10X filtered_feature_bc_matrix 目录"),
  make_option(c("-o", "--out_dir"),       type="character", default="qc_mad_plots",
              help="输出目录 [default %default]"),
  make_option(c("--nmads_feat"),          type="character", default="5",
              help="特征 MAD 阈值，可逗号分隔多个，如 3,5,7"),
  make_option(c("--nmads_mt"),            type="integer",  default=3,
              help="percent.mt 的 MAD 阈值 [default %default]"),
  make_option(c("--abs_mt_cut"),          type="double",   default=8,
              help="percent.mt 绝对上限 [default %default]"),
  make_option(c("--seed"),                type="integer",  default=1234,
              help="随机种子 [default %default]"),
  make_option(c("--max_pc"),              type="integer",  default=30,
              help="后续 PCA/UMAP 使用的 PC 数量 [default %default]"),
  make_option(c("--doublet_pN"),          type="double",   default=0.25,
              help="DoubletFinder 的 pN 参数 [default %default]"),
  make_option(c("--exp_doublet_rate"),    type="double",   default=0.008,
              help="预期双细胞率 (×每千细胞) [default %default]"),
  make_option(c("--cluster_col"),         type="character", default=NULL,
              help="已有的细胞亚群列名（可选）")
)
opt <- parse_args(OptionParser(option_list=opt_list))
if (is.null(opt$in_dir)) stop("必须提供 --in_dir")

nmads_feat_vec <- as.numeric(strsplit(opt$nmads_feat, ",")[[1]])
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

# ---------- 2. 读取数据 ----------
seu <- CreateSeuratObject(Read10X(opt$in_dir), project = "SC")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, "^MT-")

# ---------- 3. QC 工具函数 ----------
qc_meta <- function(obj, nmads_feat, nmads_mt, abs_mt_cut){
  obj@meta.data |>
    rownames_to_column("cell") |>
    mutate(log_nCount   = log1p(nCount_RNA),
           log_nFeature = log1p(nFeature_RNA)) |>
    group_by(orig.ident) |>
    mutate(
      med_cnt  = median(log_nCount),   mad_cnt  = mad(log_nCount),
      med_feat = median(log_nFeature), mad_feat = mad(log_nFeature),
      med_mt   = median(percent.mt),   mad_mt   = mad(percent.mt),
      out_cnt  = log_nCount   < med_cnt  - nmads_feat*mad_cnt |
                 log_nCount   > med_cnt  + nmads_feat*mad_cnt,
      out_feat = log_nFeature < med_feat - nmads_feat*mad_feat |
                 log_nFeature > med_feat + nmads_feat*mad_feat,
      out_mt   = percent.mt   > med_mt + nmads_mt*mad_mt |
                 percent.mt   > abs_mt_cut,
      keep = !(out_cnt | out_feat | out_mt)
    ) |>
    ungroup() |>
    mutate(keep = replace_na(keep, FALSE),
           keep = factor(keep, levels = c(FALSE, TRUE),
                               labels  = c("discard", "keep")))
}

qc_stats <- function(obj, nmads_feat, nmads_mt, abs_mt_cut){
  meta <- qc_meta(obj, nmads_feat, nmads_mt, abs_mt_cut)
  tibble(nmads_feat, nmads_mt, abs_mt_cut,
         kept        = sum(meta$keep == "keep"),
         removed     = sum(meta$keep == "discard"),
         pct_removed = 100*removed/(kept+removed))
}

# ---------- 4. 遍历不同 MAD 阈值 ----------
grid <- expand.grid(
  nmads_feat = nmads_feat_vec,
  nmads_mt   = opt$nmads_mt,
  abs_mt_cut = opt$abs_mt_cut,
  KEEP.OUT.ATTRS = FALSE)

plot_once <- function(params){
  nmf <- params$nmads_feat; nmm <- params$nmads_mt; abs <- params$abs_mt_cut
  meta <- qc_meta(seu, nmf, nmm, abs)
  seu_tmp <- seu; seu_tmp$keep <- meta$keep[match(Cells(seu_tmp), meta$cell)]

  p1 <- VlnPlot(seu_tmp,
                features = c("nCount_RNA","nFeature_RNA","percent.mt"),
                group.by = "keep", pt.size = .1, ncol = 3) +
        ggtitle(sprintf("feat=%s MAD | mt<=%s MAD / %s%%", nmf, nmm, abs))
  p2 <- FeatureScatter(seu_tmp, "nCount_RNA", "percent.mt", group.by="keep") +
        scale_color_manual(values = c("red","forestgreen")) +
        ggtitle("nCount vs percent.mt")

  ggsave(file.path(opt$out_dir,
                   sprintf("qc_feat%02d_mt%02d.pdf", nmf, abs)),
         p1/p2, width = 10, height = 8)

  qc_stats(seu, nmf, nmm, abs)
}

results <- map_dfr(split(grid, seq(nrow(grid))), plot_once)

# 总体去除比例柱形图
p_bar <- ggplot(results, aes(factor(nmads_feat), pct_removed))+
  geom_col(fill = "steelblue")+
  geom_text(aes(label = sprintf("%.1f%%", pct_removed)), vjust = -0.3)+
  labs(x = "nmads_feat", y = "% removed",
       title = "Removal rate under different MAD thresholds")+
  theme_minimal()
ggsave(file.path(opt$out_dir, "removal_bar.pdf"), p_bar, width = 5, height = 4)

# ---------- 5. 选择最佳 MAD 并过滤 ----------
best_feat <- if(length(nmads_feat_vec) > 1) median(nmads_feat_vec) else nmads_feat_vec
meta_final <- qc_meta(seu, best_feat, opt$nmads_mt, opt$abs_mt_cut)
seu <- AddMetaData(seu,
                   metadata = meta_final |>
                              select(cell, keep) |>
                              column_to_rownames("cell"))
seu_filtered <- subset(seu, subset = keep == "keep")
message(sprintf("Before: %d cells | After: %d cells | Removed: %.1f%%",
        ncol(seu), ncol(seu_filtered),
        100 * (ncol(seu) - ncol(seu_filtered)) / ncol(seu)))

# ---------- 6. SCT → PCA → UMAP ----------
seu_filtered <- SCTransform(seu_filtered, vars.to.regress="percent.mt", verbose=FALSE)
seu_filtered <- RunPCA(seu_filtered, npcs = opt$max_pc, verbose=FALSE)
seu_filtered <- RunUMAP(seu_filtered, dims = 1:opt$max_pc, verbose=FALSE)

# ---------- 7. DoubletFinder ----------
if (is.null(opt$cluster_col)) {
  seu_filtered <- FindNeighbors(seu_filtered, dims=1:opt$max_pc, verbose=FALSE) |>
                  FindClusters(resolution=0.5, verbose=FALSE)
  annotations <- Idents(seu_filtered)
} else {
  if (!opt$cluster_col %in% colnames(seu_filtered@meta.data))
    stop(sprintf("找不到 cluster_col '%s'", opt$cluster_col))
  annotations <- seu_filtered[[opt$cluster_col]][,1]
}
homotypic.prop <- modelHomotypic(annotations)
exp_poi      <- opt$exp_doublet_rate * (ncol(seu_filtered) / 1000)
nExp_poi     <- round(exp_poi * ncol(seu_filtered))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

set.seed(opt$seed)  # 保证 pANN 可重复
sweep.res   <- paramSweep(seu_filtered, PCs = 1:opt$max_pc, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn       <- find.pK(sweep.stats)
mpK         <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])

seu_filtered <- doubletFinder(
  seu_filtered,
  PCs  = 1:opt$max_pc,
  pN   = opt$doublet_pN,
  pK   = mpK,
  nExp = nExp_poi.adj,
  sct  = TRUE
)

class_col <- grep("^DF\\.classifications_", colnames(seu_filtered@meta.data),
                  value = TRUE) |> tail(1)

du <- DimPlot(seu_filtered, group.by = class_col, cols = c("blue", "red"))
ggsave(file.path(opt$out_dir, "doubletFinder.pdf"), du, width = 10, height = 8)

sc_clean <- seu_filtered[, seu_filtered@meta.data[[class_col]] == "Singlet"]

# ---------- 8. 保存结果 ----------
saveRDS(seu_filtered, file = file.path(opt$out_dir, "PBMC_qc_filtered.rds"))
saveRDS(sc_clean,     file = file.path(opt$out_dir, "PBMC_qc_filtered.doubletFinder.rds"))
message("✔ 流水线完成 — 所有文件已保存到 ", opt$out_dir)

