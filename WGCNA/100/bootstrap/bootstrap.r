load("/home/share_data1/luchh/vascular24ts/trans/WGCNA/WGCNAPto/099/wgcna-network.Rdata")
setwd("/home/share_data1/luchh/vascular24ts/trans/WGCNA/WGCNAPto/099/bootstrap/")

# -------------- 前提：你已经有以下对象在环境里 ---------------
# Expr3           : 预处理后的样本×基因表达矩阵（行是样本, 列是基因）
# moduleColors    : 原始模块的颜色标签向量，顺序对应 colnames(Expr3)
# sft             : pickSoftThreshold 结果（含 fitIndices 和 powerEstimate）
# 如果没有，可以 load("wgcna-network.Rdata") 恢复： contains Expr3, net, sft, moduleColors, etc.

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()  # 如果需要并行（可选）

# ----------------- 1. 选 optimal_power（若是 NA 时） -----------------
# 尝试从 sft$powerEstimate 取；如果是 NA 或不满足，用手动策略挑一个
if (!exists("optimal_power") || is.na(optimal_power)) {
  if (exists("sft") && !is.null(sft$powerEstimate) && !is.na(sft$powerEstimate)) {
    optimal_power <- sft$powerEstimate
    message("使用 sft$powerEstimate 作为 optimal_power: ", optimal_power)
  } else if (exists("sft") && !is.null(sft$fitIndices)) {
    fit <- sft$fitIndices
    signedR2 <- -sign(fit[,3]) * fit[,2]
    powers <- fit[,1]
    threshold <- 0.8  # 可以下调到 0.8/0.85 如果没有 power 达标
    candidates <- which(signedR2 >= threshold)
    if (length(candidates) > 0) {
      optimal_power <- powers[candidates[1]]
      message("找到满足阈值的最小 power: ", optimal_power, " (signedR2=", round(signedR2[candidates[1]],3), ")")
    } else {
      best_idx <- which.max(signedR2)
      optimal_power <- powers[best_idx]
      warning(sprintf("没有 power 达到 signed R² >= %.2f，退而求其次选 signedR2 最大的 power = %s (signedR2=%.3f)", 
                      threshold, optimal_power, signedR2[best_idx]))
    }
  } else {
    stop("缺失 sft 结果，无法确定 optimal_power。请先运行 pickSoftThreshold 再继续。")
  }
}

# 检查
if (is.na(optimal_power) || length(optimal_power) != 1) {
  stop("optimal_power 不是单一有效值，请检查 soft-threshold 结果。")
}
cat("最终使用的 optimal_power:", optimal_power, "\n")

# ----------------- 2. 定义辅助函数 -----------------
# 映射 bootstrap 模块标签到原始模块（基于最大 overlap）
matchBootToOriginal <- function(origColors, bootColors) {
  tab <- table(origColors, bootColors)
  boot2orig <- apply(tab, 2, function(col) {
    if (all(col == 0)) NA else names(which.max(col))
  })
  remapped <- bootColors
  for (bootCol in names(boot2orig)) {
    mapped <- boot2orig[[bootCol]]
    if (!is.na(mapped)) {
      remapped[bootColors == bootCol] <- mapped
    }
  }
  return(remapped)
}

# bootstrap WGCNA 稳定性计算函数（返回 gene/module stability 等）
runBootstrapWGCNA <- function(Expr3, moduleColors, power,
                              networkType = c("unsigned", "signed"),
                              corType = c("pearson", "bicor"),
                              minModuleSize = 30, mergeCutHeight = 0.25,
                              B = 100, seed = 1234, verbose = TRUE) {
  networkType <- match.arg(networkType)
  corType <- match.arg(corType)
  set.seed(seed)
  geneNames <- colnames(Expr3)
  origModules <- moduleColors
  uniqueModules <- setdiff(unique(origModules), "grey")

  bootAssigned <- matrix(NA, nrow = length(geneNames), ncol = B,
                         dimnames = list(geneNames, paste0("b", 1:B)))
  moduleJaccard <- matrix(NA, nrow = length(uniqueModules), ncol = B,
                          dimnames = list(uniqueModules, paste0("b", 1:B)))

  for (b in 1:B) {
    if (verbose) cat("Bootstrap iteration:", b, "\n")
    bootIdx <- sample(nrow(Expr3), replace = TRUE)
    Expr_boot <- Expr3[bootIdx, ]

    net_boot <- blockwiseModules(
      Expr_boot,
      power = power,
      networkType = networkType,
      corType = corType,
      TOMType = networkType,
      minModuleSize = minModuleSize,
      mergeCutHeight = mergeCutHeight,
      numericLabels = FALSE,
      verbose = 0,
      saveTOMs = FALSE
    )

    bootColors <- net_boot$colors
    remapped <- matchBootToOriginal(origColors = origModules, bootColors = bootColors)
    bootAssigned[, b] <- remapped

    for (mod in uniqueModules) {
      origGenes <- geneNames[which(origModules == mod)]
      bootGenes <- geneNames[which(remapped == mod)]
      if (length(origGenes) == 0 || length(bootGenes) == 0) {
        moduleJaccard[mod, b] <- 0
      } else {
        moduleJaccard[mod, b] <- length(intersect(origGenes, bootGenes)) / 
                                 length(union(origGenes, bootGenes))
      }
    }
  }

  geneStability <- sapply(seq_along(geneNames), function(i) {
    orig <- origModules[i]
    mean(bootAssigned[i, ] == orig, na.rm = TRUE)
  })
  names(geneStability) <- geneNames

  moduleStability <- apply(moduleJaccard, 1, median, na.rm = TRUE)

  list(
    geneStability = geneStability,
    moduleStability = moduleStability,
    bootAssigned = bootAssigned,
    moduleJaccard = moduleJaccard
  )
}

# ----------------- 3. 运行 bootstrap 估计稳定性 -----------------
# 确保 moduleColors 和 Expr3 一致
stopifnot(length(moduleColors) == ncol(Expr3))

# 设定参数
type <- "unsigned"    # networkType
corType <- "pearson"  # correlation
B <- 100              # bootstrap 轮数，可调 200~500 提高稳定性估计

bootRes <- runBootstrapWGCNA(
  Expr3 = Expr3,
  moduleColors = moduleColors,
  power = optimal_power,
  networkType = type,
  corType = corType,
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  B = B,
  seed = 2025,
  verbose = TRUE
)

# ----------------- 4. 生成 moduleStability 表格并保存 -----------------
moduleStability_df <- data.frame(
  Module = names(bootRes$moduleStability),
  Median_Jaccard = as.numeric(bootRes$moduleStability),
  stringsAsFactors = FALSE
)
# 加上每个模块基因数（含 grey 可以按需要过滤）
moduleStability_df$Size <- sapply(moduleStability_df$Module, function(m) {
  sum(moduleColors == m)
})
# 按稳定性降序
moduleStability_df <- moduleStability_df[order(-moduleStability_df$Median_Jaccard), ]

# 保存到 CSV
write.csv(moduleStability_df, "moduleStability_table.csv", row.names = FALSE, quote = FALSE)

# 简要输出预览
print(moduleStability_df)

# 可选：画条形图直观展示
barplot(moduleStability_df$Median_Jaccard, names.arg = moduleStability_df$Module,
        las=2, main="Bootstrap 模块稳定性（median Jaccard）",
        ylab="Median Jaccard", cex.names=0.8)
abline(h=0.5, col="red", lty=2)

# ----------------- 5. 生成 moduleStability gene表格并保存 -----------------

# 1. 保证 moduleColors 有名字（基因名），否则按 Expr3 列名赋
if (is.null(names(moduleColors))) {
  names(moduleColors) <- colnames(Expr3)
}

# 2. 基因级 stability 向量
geneStability <- bootRes$geneStability  # names 是基因

# 3. 构造基础表：Gene, Module, Stability
core_df <- data.frame(
  Gene = names(geneStability),
  Module = moduleColors[names(geneStability)],
  Stability = as.numeric(geneStability),
  stringsAsFactors = FALSE
)

# 4. （可选）加上 kME：基因与所属模块 eigengene 的相关性
# 先计算 module eigengenes（如果还没算）
MEs <- moduleEigengenes(Expr3, moduleColors)$eigengenes
geneMM <- as.data.frame(cor(Expr3, MEs, use = "p"))  # gene-module correlations

# 取每个基因在其自己的模块的 kME
core_df$kME <- mapply(function(g, m) {
  eigName <- paste0("ME", m)
  if (eigName %in% colnames(geneMM)) {
    geneMM[g, eigName]
  } else {
    NA
  }
}, core_df$Gene, core_df$Module)

# 5. 设稳定性阈值并筛选 core genes
stability_threshold <- 0
core_df_filtered <- subset(core_df, Stability >= stability_threshold)

# 6. 结果排序：先模块再 stability 降序
core_df_filtered <- core_df_filtered[order(core_df_filtered$Module, -core_df_filtered$Stability), ]

# 7. 写出 CSV（包含 gene, module, stability, kME）
write.csv(core_df_filtered, "0Bootstrap_Core_Genes_by_Stability.csv", row.names = FALSE, quote = FALSE)

# 8. 简要输出统计
cat("每个模块 high-stability gene 数量：\n")
print(table(core_df_filtered$Module))
cat("前几行预览：\n")
print(head(core_df_filtered))
