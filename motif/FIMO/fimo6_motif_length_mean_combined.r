#!/usr/bin/env Rscript

# Robust heatmap of motif × promoter length mean combined_score.
# Usage: edit input_file / out_prefix below, then run: Rscript plot_mean_combined_heatmap_fixed.R

suppressWarnings({
  options(repos = "https://cloud.r-project.org")
  needed <- c("pheatmap", "RColorBrewer", "grid")
  for (pkg in needed) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
})
library(pheatmap)
library(RColorBrewer)
library(grid)

# ---- USER CONFIGURE ----
input_file <- "/home/share_data1/luchh/vascular24ts/STREME/enrichment/100new/vascula_mean_combined_output/motif_length_mean_combined.tsv"  # <- 改成你的实际路径
out_prefix <- "results/100newmean_combined_heatmap"                        # <- 输出前缀（不带扩展名）
scale_rows <- TRUE  # 如果想每个 motif 做 z-score 标准化设 TRUE；否则 FALSE
cluster_rows <- TRUE
cluster_cols <- TRUE
# ------------------------

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

# Create output dir if needed
outdir <- dirname(out_prefix)
if (!dir.exists(outdir) && nzchar(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Read matrix
cat("Reading matrix from", input_file, "\n")
mat <- tryCatch({
  as.data.frame(read.table(input_file, header = TRUE, sep = "\t", row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE))
}, error = function(e) {
  stop("Failed to read input matrix: ", e$message)
})

mat_num <- as.matrix(mat)
mode(mat_num) <- "numeric"
mat_num[is.infinite(mat_num)] <- NA  # replace Inf with NA for plotting

# Row scaling if requested (z-score per motif)
if (scale_rows) {
  cat("Scaling rows (motif-wise z-score)\n")
  row_means <- rowMeans(mat_num, na.rm = TRUE)
  row_sds <- apply(mat_num, 1, sd, na.rm = TRUE)
  # avoid division by zero: if sd==0, result will be NA
  mat_num <- t((t(mat_num) - row_means) / row_sds)
}

# Setup color and breaks
if (scale_rows) {
  # diverging centered at 0
  max_abs <- max(abs(mat_num), na.rm = TRUE)
  lim <- ceiling(max_abs * 10) / 10  # round up a bit
  breaks <- seq(-lim, lim, length.out = 101)
  colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
} else {
  rng <- range(mat_num, na.rm = TRUE)
  if (diff(rng) == 0) {
    # flat, make small range to avoid error
    rng <- c(rng[1] - 1e-3, rng[2] + 1e-3)
  }
  breaks <- seq(rng[1], rng[2], length.out = 101)
  colors <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
}

# Plotting
cat("Generating heatmap\n")
ph <- tryCatch({
  pheatmap(
    mat_num,
    color = colors,
    breaks = breaks,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 4,
    border_color = NA,
    main = "Motif × Promoter Length mean combined_score",
    na_col = "grey90",
    silent = TRUE
  )
}, error = function(e) {
  stop("Failed to create heatmap: ", e$message)
})

# Save PNG
png_file <- paste0(out_prefix, ".png")
cat("Saving PNG to", png_file, "\n")
png(filename = png_file, width = 1600, height = 2000, res = 200)
grid::grid.draw(ph$gtable)
dev.off()

# Save PDF
pdf_file <- paste0(out_prefix, ".pdf")
cat("Saving PDF to", pdf_file, "\n")
# Use cairo_pdf if available for better compatibility
if (capabilities("cairo")) {
  cairo_pdf(filename = pdf_file, width = 20, height = 40)
} else {
  pdf(filename = pdf_file, width = 20, height = 40)
}
grid::grid.draw(ph$gtable)
dev.off()

cat("Done. Outputs:\n", png_file, "\n", pdf_file, "\n")
