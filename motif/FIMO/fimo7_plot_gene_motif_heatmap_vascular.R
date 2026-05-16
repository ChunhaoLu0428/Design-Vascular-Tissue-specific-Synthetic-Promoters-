#!/usr/bin/env Rscript

# Heatmap for vascular gene × motif mean combined_score.
# Input is gene_motif_mean_combined_vascular.tsv (rows=vascular genes, cols=motifs).

suppressPackageStartupMessages({
  options(repos = "https://cloud.r-project.org")
  needed <- c("pheatmap", "RColorBrewer", "matrixStats")
  for (pkg in needed) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
})

# ---- Defaults (can be overridden via CLI) ----
input_file <- "/home/share_data1/luchh/vascular24ts/STREME/enrichment/100new/vascula_mean_combined_output/gene_motif_mean_combined_vascular.tsv"
out_prefix <- "results/100newvascular_gene_motif_heatmap"
scale_rows <- TRUE    # z-score per gene
scale_cols <- FALSE   # z-score per motif
top_genes <- 40000      # keep top N genes by variance, set NA to keep all
top_motifs <- 3000     # keep top N motifs by variance, set NA to keep all
cluster_rows <- TRUE
cluster_cols <- TRUE
# ---------------------------------------------

# parse simple CLI overrides like --input=..., --top-genes=100 etc.
args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (grepl("^--input=", arg)) input_file <- sub("^--input=", "", arg)
  if (grepl("^--out-prefix=", arg)) out_prefix <- sub("^--out-prefix=", "", arg)
  if (grepl("^--no-scale-rows$", arg)) scale_rows <- FALSE
  if (grepl("^--scale-cols$", arg)) scale_cols <- TRUE
  if (grepl("^--top-genes=", arg)) top_genes <- as.integer(sub("^--top-genes=", "", arg))
  if (grepl("^--top-motifs=", arg)) top_motifs <- as.integer(sub("^--top-motifs=", "", arg))
}

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

# Read matrix
mat <- read.table(input_file, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)
mat <- as.matrix(mat)
mode(mat) <- "numeric"
mat[is.infinite(mat)] <- NA

# Filter top genes/motifs by variance if requested
if (!is.null(top_genes) && !is.na(top_genes) && top_genes < nrow(mat)) {
  row_var <- rowVars(mat, na.rm = TRUE)
  keep_genes <- order(row_var, decreasing = TRUE)[1:top_genes]
  mat <- mat[keep_genes, , drop = FALSE]
}
if (!is.null(top_motifs) && !is.na(top_motifs) && top_motifs < ncol(mat)) {
  col_var <- colVars(mat, na.rm = TRUE)
  keep_motifs <- order(col_var, decreasing = TRUE)[1:top_motifs]
  mat <- mat[, keep_motifs, drop = FALSE]
}

# Scaling
if (scale_rows) {
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, sd, na.rm = TRUE)
  mat <- sweep(mat, 1, row_means, FUN = "-")
  mat <- sweep(mat, 1, row_sds, FUN = "/")
}
if (scale_cols) {
  col_means <- colMeans(mat, na.rm = TRUE)
  col_sds <- apply(mat, 2, sd, na.rm = TRUE)
  mat <- sweep(mat, 2, col_means, FUN = "-")
  mat <- sweep(mat, 2, col_sds, FUN = "/")
}
mat[is.infinite(mat)] <- NA

# Color scheme
if (scale_rows || scale_cols) {
  maxabs <- max(abs(mat), na.rm = TRUE)
  breaks <- seq(-maxabs, maxabs, length.out = 101)
  colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
} else {
  rng <- range(mat, na.rm = TRUE)
  if (diff(rng) == 0) rng <- c(rng[1] - 1e-3, rng[2] + 1e-3)
  breaks <- seq(rng[1], rng[2], length.out = 101)
  colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
}

# Ensure output dir
outdir <- dirname(out_prefix)
if (nchar(outdir) && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Draw and save heatmap
png_file <- paste0(out_prefix, ".png")
pdf_file <- paste0(out_prefix, ".pdf")

pheatmap(
  mat,
  color = colors,
  breaks = breaks,
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 6,
  main = "Vascular gene × motif mean combined_score",
  na_col = "grey90",
  filename = png_file,
  width = 50,
  height = 100
)

# also save PDF (re-draw to capture same clustering/layout)
library(grid)
ph <- pheatmap(
  mat,
  color = colors,
  breaks = breaks,
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 6,
  main = "Vascular gene × motif mean combined_score",
  na_col = "grey90",
  silent = TRUE
)
pdf(pdf_file, width = 50, height = 100)
grid::grid.draw(ph$gtable)
dev.off()

message("Saved heatmap to:", png_file, "and", pdf_file)
