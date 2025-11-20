#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---------------------------
# 0) 解析命令行参数
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

usage <- "
Usage:
  Rscript run_weighted_vln_celltype.R \
    --rds Seurat_obj.rds \
    --genes gene_list.tsv \
    --group_by cell_type \
    --outdir weighted_vln_out \
    [--ymax 3] [--min_expr 1] [--width 4] [--height 2] \
    [--cluster_order \"A,B,C\"] [--cluster_order_file order.txt] \
    [--palette_file palette.tsv] [--export_weights FALSE]

gene_list.tsv formats:
  (1) one column, no header:  GENE_ID
  (2) two columns, header optional:
        gene_id   gene_label
      e.g. AT4G21750 ATML1

palette.tsv format (optional):
  group   color
  CD8 T   #4E79A7
  B cells #59A14F
"

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste("No value for", flag))
  args[hit + 1]
}

rds_file            <- get_arg("--rds")
genes_file          <- get_arg("--genes")
group_by            <- get_arg("--group_by", "cell_type")
outdir              <- get_arg("--outdir", "weighted_vln_out")
ymax                <- as.numeric(get_arg("--ymax", "3"))
min_expr            <- as.numeric(get_arg("--min_expr", "1"))
width               <- as.numeric(get_arg("--width", "4"))
height              <- as.numeric(get_arg("--height", "2"))
cluster_order_str   <- get_arg("--cluster_order", NULL)
cluster_order_file  <- get_arg("--cluster_order_file", NULL)
palette_file        <- get_arg("--palette_file", NULL)
export_weights      <- tolower(get_arg("--export_weights", "false")) %in% c("true","t","1","yes")

if (is.null(rds_file) || is.null(genes_file)) {
  cat(usage)
  quit(save="no", status=1)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("rds_file   = ", rds_file)
message("genes_file = ", genes_file)
message("group_by   = ", group_by)
message("outdir     = ", outdir)

# ---------------------------
# 1) 读入 Seurat 对象
# ---------------------------
seu <- readRDS(rds_file)

if (group_by != "ident" && !group_by %in% colnames(seu@meta.data)) {
  stop("group_by = ", group_by, " 不在 seu@meta.data 中。")
}

# ---------------------------
# 2) 读入 gene list
# ---------------------------
# 允许 1 列或 2 列；允许有无表头
genes_df <- read_tsv(genes_file, col_names = FALSE, comment = "#", show_col_types = FALSE)

if (ncol(genes_df) == 1) {
  colnames(genes_df) <- c("gene_id")
  genes_df$gene_label <- genes_df$gene_id
} else {
  # 前两列分别当 gene_id / gene_label
  genes_df <- genes_df[, 1:2]
  colnames(genes_df) <- c("gene_id","gene_label")
}

genes_df <- genes_df %>% distinct()

# ---------------------------
# 3) cluster 顺序（可选）
# ---------------------------
cluster_order <- NULL
if (!is.null(cluster_order_file)) {
  cluster_order <- read_lines(cluster_order_file)
  cluster_order <- cluster_order[cluster_order != ""]
} else if (!is.null(cluster_order_str)) {
  cluster_order <- strsplit(cluster_order_str, ",")[[1]]
  cluster_order <- trimws(cluster_order)
}

# ---------------------------
# 4) palette（可选）
# ---------------------------
cluster_palette <- NULL
if (!is.null(palette_file)) {
  pal_df <- read_tsv(palette_file, show_col_types = FALSE)
  if (!all(c("group","color") %in% colnames(pal_df))) {
    stop("palette_file 需要两列：group, color")
  }
  cluster_palette <- pal_df$color
  names(cluster_palette) <- pal_df$group
}

# ---------------------------
# 5) 你的加权 violin 函数（按 group_by 支持 meta.data）
# ---------------------------
WeightedVln_clip3 <- function(
  seu,
  gene_at,
  group.by        = "ident",
  file            = "gene_weighted_vln_clip3.pdf",
  gene_label      = NULL,
  ymax            = 3,
  width           = 4,
  height          = 0.8,
  cluster_palette = NULL,
  min_expr        = 1,
  cluster_order   = NULL,
  export_weights  = FALSE
) {

  if (!gene_at %in% rownames(seu)) {
    stop(paste("基因", gene_at, "不在 Seurat features 中"))
  }
  if (is.null(gene_label)) gene_label <- gene_at

  ## 1) 取表达 + cluster 信息（全量）
  df <- FetchData(seu, vars = c(gene_at, group.by))
  colnames(df) <- c("expr", "cluster")

  ## 2) 阈值过滤：小于 min_expr 视为 0（但不再 stop）
  df$expr <- ifelse(df$expr >= min_expr, df$expr, 0)
  if (all(df$expr == 0)) {
    message(sprintf(
      "[%s] 没有细胞 expr >= %.2f，将输出空图（保留所有 cluster 位置）",
      gene_label, min_expr
    ))
  }

  ## 3) 确定“全量 cluster levels”
  cluster_chr <- as.character(df$cluster)

  if (!is.null(cluster_order)) {
    lev_all <- cluster_order[cluster_order %in% unique(cluster_chr)]
    if (length(lev_all) == 0) stop("cluster_order 与数据无匹配。")
  } else {
    # 如果 meta 是 factor，用其 levels；否则按字母顺序
    if (group.by != "ident" && is.factor(seu@meta.data[[group.by]])) {
      lev_all <- levels(seu@meta.data[[group.by]])
      lev_all <- lev_all[lev_all %in% unique(cluster_chr)]
    } else {
      lev_all <- sort(unique(cluster_chr))
    }
  }

  df$cluster <- factor(cluster_chr, levels = lev_all)

  ## 4) 每个 cluster 的阳性比例 pct_pos（>= min_expr 视为阳性）
  cluster_stats <- df %>%
    group_by(cluster) %>%
    summarise(
      n_cells = n(),
      n_pos   = sum(expr > 0),
      pct_pos = n_pos / n_cells,
      .groups = "drop"
    )

  ## 5) 权重：如果 max pct_pos = 0（全无表达），weight 全为 0
  max_pct <- max(cluster_stats$pct_pos)
  if (max_pct == 0) {
    cluster_stats$weight <- 0
  } else {
    cluster_stats <- cluster_stats %>%
      mutate(weight = pct_pos / max_pct)
  }

  ## 可选导出权重
  if (export_weights) {
    write_tsv(cluster_stats, sub("\\.pdf$", ".weights.tsv", file))
  }

  ## 6) 合并权重，计算加权表达
  df_w <- df %>%
    left_join(cluster_stats[, c("cluster", "weight")], by = "cluster") %>%
    mutate(expr_weighted = expr * weight)

  ## 7) 只用阳性细胞画“点云”
  df_pos <- df_w %>% filter(expr > 0)

  ## 8) 用阳性细胞算 violin；若某些 cluster 无阳性，补 dummy 0
  df_plot <- df_pos
  miss_lev <- setdiff(lev_all, unique(as.character(df_pos$cluster)))
  if (length(miss_lev) > 0) {
    dummy <- data.frame(
      expr          = 0,
      cluster       = factor(miss_lev, levels = lev_all),
      weight        = cluster_stats$weight[match(miss_lev, as.character(cluster_stats$cluster))],
      expr_weighted = 0
    )
    df_plot <- bind_rows(df_plot, dummy)
  }

  ## 9) 裁剪到 [0, ymax]
  df_plot <- df_plot %>%
    mutate(expr_plot = pmin(expr_weighted, ymax))

  ## 10) 中位数（只对阳性算）
  df_med <- df_pos %>%
    group_by(cluster) %>%
    summarise(med = median(pmin(expr_weighted, ymax), na.rm = TRUE), .groups="drop")

  if (nrow(df_med) > 0) {
    df_med$x     <- as.numeric(df_med$cluster)
    df_med$x_min <- df_med$x - 0.18
    df_med$x_max <- df_med$x + 0.18
  }

  ## 11) 颜色映射
  scale_col_layer <- NULL
  if (!is.null(cluster_palette)) {
    if (!is.null(names(cluster_palette)) && all(lev_all %in% names(cluster_palette))) {
      cols <- cluster_palette[lev_all]
    } else {
      cols <- cluster_palette[seq_len(min(length(cluster_palette), length(lev_all)))]
      names(cols) <- lev_all[seq_along(cols)]
    }
    scale_col_layer <- scale_colour_manual(values = cols, drop = FALSE)
  }

  ## 12) 作图
  pdf(file, width = width, height = height)

  p <- ggplot(df_plot, aes(x = cluster, y = expr_plot)) +
    geom_violin(
      aes(group = cluster),
      scale="width", width=0.5, trim=TRUE,
      fill="grey90", colour=NA
    ) +
    # 点云只画阳性细胞
    geom_point(
      data = df_pos,
      aes(colour = cluster, y = pmin(expr_weighted, ymax)),
      position = position_jitter(width=0.3, height=0),
      size=0.28, alpha=0.35
    ) +
    geom_violin(
      aes(group = cluster),
      scale="width", width=0.5, trim=TRUE,
      fill=NA, colour="black", linewidth=0.22
    ) +
    # 中位数线：仅阳性 cluster 才有
    {if (nrow(df_med) > 0)
      geom_segment(
        data=df_med,
        aes(x=x_min, xend=x_max, y=med, yend=med),
        inherit.aes=FALSE, colour="black", linewidth=0.25, lineend="round"
      )
    } +
    scale_y_continuous(
      limits=c(0, ymax),
      breaks=seq(0, ymax, by=1),
      labels=as.character(seq(0, ymax, by=1))
    ) +
    scale_x_discrete(drop = FALSE) +   # ✅ 核心：保留空 cluster 位置
    labs(x = group.by, y = "weighted expression (clipped)") +
    theme_classic(base_size=7) +
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      legend.position="none",
      axis.line.x = element_line(colour="grey75", linetype="dashed", linewidth=0.5),
      axis.line.y = element_line(colour="black", linewidth=0.5),
      axis.ticks  = element_line(linewidth=0.5)
    )

  if (!is.null(scale_col_layer)) p <- p + scale_col_layer

  print(p)
  dev.off()
}

# ---------------------------
# 6) 循环每个基因输出
# ---------------------------
sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

for (i in seq_len(nrow(genes_df))) {
  gene_id    <- genes_df$gene_id[i]
  gene_label <- genes_df$gene_label[i]

  out_pdf <- file.path(outdir, paste0(sanitize_filename(gene_label),
                                     "_weighted_vln_clip3.pdf"))

  message("Plotting gene: ", gene_label, " (", gene_id, ") -> ", out_pdf)

  tryCatch({
    WeightedVln_clip3(
      seu             = seu,
      gene_at         = gene_id,
      group.by        = group_by,      # ✅ 这里用 cell_type 或任意 meta 列
      file            = out_pdf,
      gene_label      = gene_label,
      ymax            = ymax,
      width           = width,
      height          = height,
      cluster_palette = cluster_palette,
      min_expr        = min_expr,
      cluster_order   = cluster_order,
      export_weights  = export_weights
    )
  }, error = function(e) {
    message("  [skip] ", gene_label, " : ", e$message)
  })
}

message("Done.")

