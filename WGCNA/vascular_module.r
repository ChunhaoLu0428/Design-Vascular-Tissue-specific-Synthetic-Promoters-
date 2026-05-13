library(WGCNA)
library(dplyr)
library(pheatmap)

# --------- 1. 加载 WGCNA 结果和原始表达 ----------
load("wgcna-network.Rdata")  # 期望包含 Expr3, net/moduleColors

# 恢复 moduleColors（如果是 numeric label 先转换成颜色）
if (!exists("moduleColors")) {
  if (exists("net") && !is.null(net$colors)) {
    moduleColors <- labels2colors(net$colors)
  } else {
    stop("找不到 moduleColors，也不能从 net 推断。")
  }
}
# genes actually used in WGCNA
genesUsed <- colnames(Expr3)  # Expr3: samples x genes
# 保证 moduleColors 以基因命名并对齐
if (length(moduleColors) != length(genesUsed)) stop("moduleColors 长度和 Expr3 基因数不一致。")
names(moduleColors) <- genesUsed

# 读原始表达（若之前没加载）
Expr1 <- read.table("/home/share_data1/luchh/vascular24ts/1_tau_Ptogenes_expression.txt",
                    header = TRUE, row.names = 1, check.names = FALSE)
Expr2 <- log2(Expr1 + 1)  # genes x tissues

# 取出 WGCNA 实际用过并在原始表达中都有的基因
common_genes <- intersect(genesUsed, rownames(Expr2))
if (length(common_genes) == 0) stop("没有交集基因，请检查基因命名一致性。")
# 子集化并对齐
Expr2_filtered <- Expr2[common_genes, , drop = FALSE]              # genes x tissues
moduleColors_filtered <- moduleColors[common_genes]                # named by gene

# 确保顺序一致
moduleColors_filtered <- moduleColors_filtered[rownames(Expr2_filtered)]
stopifnot(all(names(moduleColors_filtered) == rownames(Expr2_filtered)))

# --------- 2. 定义维管组织（根据你的列名） ----------
vascular_tissues <- c("Mature_xylem", "Cambium", "Immature_xylem", "Phloem")
all_tissues <- colnames(Expr2_filtered)
# 校正拼写差（如果原始列名里是 Cabium 而非 Cambium）
if (!"Cambium" %in% all_tissues && "Cabium" %in% all_tissues) {
  vascular_tissues[vascular_tissues == "Cambium"] <- "Cabium"
}

nonvascular_tissues <- setdiff(all_tissues, vascular_tissues)
if (length(nonvascular_tissues) == 0) {
  warning("没有非维管组织可对比，后续判定会基于相对表达。")
}

# --------- 3. 基因内部标准化（z-score across tissues） ----------
expr_mat <- as.matrix(Expr2_filtered)  # genes x tissues
gene_means <- rowMeans(expr_mat, na.rm = TRUE)
gene_sds <- apply(expr_mat, 1, sd, na.rm = TRUE)
gene_sds[gene_sds == 0] <- 1e-6  # 防止除 0
zscore_mat <- sweep(expr_mat, 1, gene_means, "-")
zscore_mat <- sweep(zscore_mat, 1, gene_sds, "/")  # genes x tissues

# --------- 4. 判定每个基因在每个 vascular tissue 上“特异高表达” ----------
gene_flag_list <- list()
for (t in vascular_tissues) {
  if (!t %in% all_tissues) next
  expr_t <- expr_mat[, t]  # 当前 vascular tissue 表达 (vector length = genes)

  # 非维管中每基因最大表达（fallback：若没有非维管就用除去该 tissue 的最大）
  if (length(nonvascular_tissues) > 0) {
    max_nonv <- apply(expr_mat[, nonvascular_tissues, drop = FALSE], 1, max, na.rm = TRUE)
  } else {
    # 用除本组织之外的最大
    others <- setdiff(colnames(expr_mat), t)
    max_nonv <- apply(expr_mat[, others, drop = FALSE], 1, max, na.rm = TRUE)
  }

  logFC <- expr_t - max_nonv  # log2 scale 差值
  z_t <- zscore_mat[, t]

  # 定义上调/特异：logFC >= 1 (>=2-fold vs 非维管最高) 且 zscore >=1.5
  up_flag <- (logFC >= 0) & (z_t >= 1.5)

  gene_flag_list[[t]] <- data.frame(
    Gene = rownames(expr_mat),
    VascularTissue = t,
    Expr_in_vascular = expr_t,
    Max_NonV = max_nonv,
    logFC = logFC,
    Zscore = z_t,
    UpInThisTissue = up_flag,
    stringsAsFactors = FALSE
  )
}
gene_DE_df <- bind_rows(gene_flag_list)

# --------- 5. 归到 module 级别汇总 ----------
gene_module_map <- data.frame(
  Gene = rownames(expr_mat),
  Module = moduleColors_filtered,
  stringsAsFactors = FALSE
)
merged <- gene_DE_df %>%
  inner_join(gene_module_map, by = "Gene")

module_summary <- merged %>%
  group_by(Module, VascularTissue) %>%
  summarise(
    ModuleSize = n(),
    UpGenes = sum(UpInThisTissue, na.rm = TRUE),
    PropUp = UpGenes / ModuleSize,
    .groups = "drop"
  )

# --------- 6. 判定 vascular module（可调阈值） ----------
prop_threshold <- 0.2   # 至少 20% 的基因在该维管组织“特异上调”
min_up_genes <- 5      # 至少 5 个这样的基因
module_summary <- module_summary %>%
  mutate(IsVascularInThatTissue = (PropUp >= prop_threshold & UpGenes >= min_up_genes))

vascular_module_call <- module_summary %>%
  group_by(Module) %>%
  summarise(
    PassedTissues = paste(VascularTissue[IsVascularInThatTissue], collapse = ";"),
    IsVascularModule = any(IsVascularInThatTissue),
    .groups = "drop"
  )

# --------- 7. 保存与输出 ----------
write.csv(module_summary, "module_gene_level_vascular_summary.csv", row.names = FALSE)
write.csv(vascular_module_call, "vascular_module_calls.csv", row.names = FALSE)

# 打印被判为 vascular module 的
print(subset(vascular_module_call, IsVascularModule))

# --------- 8. （可选）可视化某个 module 在所有组织上的热图 ---------
# 例如看 red module
target_module <- "red"
genes_in_mod <- gene_module_map$Gene[gene_module_map$Module == target_module]
library(ggplot2)
library(pheatmap)
library(grid)

# 画图但不用 filename，让它画在当前设备
p <- pheatmap(
  mat,
  scale = "row",
  main = paste0("Module ", target_module, ": gene expression across tissues"),
  fontsize_row = 4,
  fontsize_col = 8,
  silent = TRUE  # 抑制直接绘制
)

# 把 pheatmap 的结果画到一个新图形中再保存
pdf_file <- paste0("Module_", target_module, "_expression_heatmap.pdf")
pdf(pdf_file, width = 6, height = 8)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

#——————————————————————————————————————————#
library(pheatmap)

pheatmap(
  mat,
  scale = "row",
  main = paste0("Module ", target_module, ": gene expression across tissues"),
  fontsize_row = 4,
  fontsize_col = 8,
  filename = paste0("Module_", target_module, "_expression_heatmap.pdf"),  # 直接存
  width = 6,
  height = 8
)
