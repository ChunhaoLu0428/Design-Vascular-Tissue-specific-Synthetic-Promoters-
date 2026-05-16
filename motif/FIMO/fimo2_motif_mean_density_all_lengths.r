#!/usr/bin/env Rscript

# 依赖：如果还没装，先运行一次下面两行（可手动在 R 里执行）
# install.packages("pheatmap")
# install.packages("data.table")

library(pheatmap)
library(data.table)

input_tsv <- "/home/share_data1/luchh/vascular24ts/STREME/pto_potrihomo/enrichment/results_presence_density/motif_mean_density_all_lengths.tsv"
output_png <- "motif_density_heatmap_allmotifs_r.png"
output_pdf <- "motif_density_heatmap_allmotifs_r.pdf"

# ---- 读数据 ----
df <- fread(input_tsv, data.table = FALSE)
if (!"motif" %in% colnames(df)) {
  rownames(df) <- df[, 1]
  df <- df[, -1]
} else {
  rownames(df) <- df$motif
  df$motif <- NULL
}

# ---- 筛选 top N 变异最大的 motif ----
top_n <- 500
vars <- apply(df, 1, function(x) var(x, na.rm = TRUE))
if (length(vars) > top_n) {
  sel <- names(sort(vars, decreasing = TRUE))[1:top_n]
  mat <- df[sel, , drop = FALSE]
} else {
  mat <- df
}

# ---- log10 变换 ----
eps <- 1e-9
mat_log <- log10(mat + eps)
mat_plot <- mat_log

# ---- 绘制 PNG ----
png(output_png, width = 1200, height = 2000, res = 150)
pheatmap(
  mat_plot,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 10,
  main = "Motif mean density heatmap (log10)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  border_color = NA,
  angle_col = 45
)
dev.off()

# ---- 绘制 PDF ----
pdf(output_pdf, width = 10, height = 20)
pheatmap(
  mat_plot,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 10,
  main = "Motif mean density heatmap (log10)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  border_color = NA,
  angle_col = 45
)
dev.off()

message("Wrote heatmap to ", output_png, " and ", output_pdf)



input_tsv <- "/home/share_data1/luchh/vascular24ts/STREME/pto_potrihomo/enrichment/results_presence_density/qv_motif_mean_density_all_lengths.tsv"
output_png <- "qv_motif_density_heatmap_allmotifs_r.png"
output_pdf <- "qv_motif_density_heatmap_allmotifs_r.pdf"

# ---- 读数据 ----
df <- fread(input_tsv, data.table = FALSE)
if (!"motif" %in% colnames(df)) {
  rownames(df) <- df[, 1]
  df <- df[, -1]
} else {
  rownames(df) <- df$motif
  df$motif <- NULL
}

df[is.na(df)] <- 0


# ---- 筛选 top N 变异最大的 motif ----
top_n <- 500
vars <- apply(df, 1, function(x) var(x, na.rm = TRUE))
if (length(vars) > top_n) {
  sel <- names(sort(vars, decreasing = TRUE))[1:top_n]
  mat <- df[sel, , drop = FALSE]
} else {
  mat <- df
}

# ---- log10 变换 ----
eps <- 1e-9
mat_log <- log10(mat + eps)
mat_plot <- mat_log

# ---- 绘制 PNG ----
png(output_png, width = 1200, height = 2000, res = 150)
pheatmap(
  mat_plot,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 10,
  main = "qv Motif mean density heatmap (log10)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  border_color = NA,
  angle_col = 45
)
dev.off()

# ---- 绘制 PDF ----
pdf(output_pdf, width = 10, height = 20)
pheatmap(
  mat_plot,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 10,
  main = "qv Motif mean density heatmap (log10)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  border_color = NA,
  angle_col = 45
)
dev.off()

message("Wrote heatmap to ", output_png, " and ", output_pdf)