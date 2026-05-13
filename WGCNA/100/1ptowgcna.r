# 验证安装
library(WGCNA) 
library(flashClust)
library(reshape2)
library(stringr)
library(BiocParallel)
library(pheatmap) 

setwd("/home/share_data1/luchh/vascular24ts/trans/WGCNA/WGCNAPto/100/")
rm(list = ls())
Expr1 <- read.table("/home/share_data1/luchh/vascular24ts/1_tau_Ptogenes_expression.txt", header = TRUE, row.names = 1)
Expr2 <- log(Expr1+1, 2)
dim(Expr2)
head(rownames(Expr2))
head(colnames(Expr2))

corType = "pearson" # correlation method,another correlation type is “bicor”
maxPOutliers = 0.05
robustY = FALSE  # Dealing with binary data
Expr_t <- as.data.frame(t(Expr2))

gsg = goodSamplesGenes(Expr_t, verbose = 5)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(Expr2)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(Expr2)[!gsg$goodSamples],collapse = ",")));
  # Remove the offending genes and samples from the data:
  Expr_t = Expr_t[gsg$goodSamples, gsg$goodGenes]
}

#----------将删除的基因输出-----------#

if (!gsg$allOK) {
  # 输出被移除的基因到CSV
  if (sum(!gsg$goodGenes) > 0) {
    removed_genes <- data.frame(Gene = names(Expr_t)[!gsg$goodGenes])
    write.csv(removed_genes, "removed_genes.csv", row.names = FALSE)
  }
  
  # 输出被移除的样本到CSV
  if (sum(!gsg$goodSamples) > 0) {
    removed_samples <- data.frame(Sample = rownames(Expr_t)[!gsg$goodSamples])
    write.csv(removed_samples, "removed_samples.csv", row.names = FALSE)
  }
}
#-------------------------------------#

nGenes = ncol(Expr_t)
nSamples = nrow(Expr_t)
dim(Expr_t)

gsg$allOK

#-----------------------#
gsg = goodSamplesGenes(Expr_t, verbose = 5)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(Expr2)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(Expr2)[!gsg$goodSamples],collapse = ",")));
  # Remove the offending genes and samples from the data:
  Expr_t = Expr_t[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(Expr_t)
nSamples = nrow(Expr_t)
dim(Expr_t)

gsg$allOK
#----------------------#

sampleTree <- hclust(dist(Expr_t), method = "average")
par(mar = c(0,4,2,0))
pdf(file = "sampleTree.pdf", width = 12, height = 9)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
# call outlier samples
clust <- cutreeStatic(sampleTree, cutHeight = 260, minSize = 8)
rownames(Expr_t)[clust==0]
keepSamples <- (clust != 0)
Expr3 <- Expr_t[keepSamples, ]
dim(Expr3)
geneNames <- colnames(Expr3)


powers <- c(seq(1, 10, by=1), seq(12, 30, by=2))
type ="unsigned"    
sft = pickSoftThreshold(Expr3, powerVector=powers, networkType=type, verbose=5)
sizeGrWindow (9,5)
par(mfrow= c(1,2))   
cex1=0.9        
pdf("wgcna_soft.thresholding.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.9,col="red")
dev.off()
power = sft$powerEstimate 
power

#----------------检查网络拓扑-----------#
# 初始化结果存储
results <- data.frame(
  power = 2:30,
  mean.k = numeric(29),
  R.squared = numeric(29)
)

# 循环计算不同 power 下的指标
for (i in 2:30) {
  # 计算邻接矩阵
  adjacency <- adjacency(
    datExpr = Expr3,
    power = i,
    type = "unsigned",
    corFnc = "cor",
    corOptions = list(use = 'p')
  )
  
  # 计算节点连接度
  k <- colSums(adjacency) - 1
  
  # 计算 R²（无标度拟合指数）
  logk <- log10(k + 1)
  logPk <- log10(rank(-k) / length(k))
  fit <- lm(logPk ~ logk)
  R.squared <- summary(fit)$r.squared
  
  # 计算平均连接度
  mean.k <- mean(k)
  
  # 存储结果
  results[i-1, "mean.k"] <- mean.k
  results[i-1, "R.squared"] <- R.squared
}
# 保存结果
write.csv(results, file = "WGCNA_power_analysis_results.csv", row.names = FALSE)

# 找到 R.squared 最大值对应的行
optimal_row <- results[which.max(results$R.squared), ]

# 提取对应的 Power 值
optimal_power <- optimal_row$power

# 输出结果（注意使用英文逗号）
cat("Selected Power:", optimal_power, 
    "(R² =", round(optimal_row$R.squared, 3), 
    ", mean.k =", round(optimal_row$mean.k, 1), ")\n")
    

net <- blockwiseModules(
  Expr3, 
  maxBlockSize = dim(Expr3)[2],
  corType = corType,
  power = optimal_power,
  networkType = type,
  TOMType = type,
  saveTOMs = TRUE, 
  saveTOMFileBase = "blockwiseTOM",
  minModuleSize = 30,
  mergeCutHeight =0.25,    
  numericLabels = F, # modudule named in number
  nThreads = 0, 
  verbose = 3)
#——————————————————————————#

table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
pdf("plotDendroAndColors222.pdf")
plotDendroAndColors(net$dendrograms[[1]],moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
save(Expr3, sft, net, moduleColors, file = "wgcna-network.Rdata")##如果断了，可以直接加载这个Rdata，不用再跑一遍net和sft，但需要重新加载数据之类的。

MEs = net$MEs
MEs = moduleEigengenes(Expr3, moduleColors)$eigengenes
MET = orderMEs(MEs)
sizeGrWindow(7, 6)
pdf("module_correlation333.pdf")
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

module_colors <- setdiff(unique(moduleColors), "grey")
for (color in module_colors){module <- geneNames[which(moduleColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE) }
# Export the network into edge and node list files Cytoscape can read
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dimnames(TOM) <- list(geneNames, geneNames)


for(i in module_colors)
{
  modules = i
  probes = colnames(Expr3)
  inModule = is.finite(match(moduleColors, modules))
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("cyt_edges_", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile=paste("cyt_nodes_", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,threshold = 0.1, nodeNames = modProbes, nodeAttr = moduleColors[inModule]) 
}

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(Expr3, MET, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
             as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = pearsonAndPvalue(Expr3, MET, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$pearson
  MMPvalue   = geneModuleMembershipA$p
}
write.table(geneModuleMembership,file="geneModuleMembership")
write.table(MMPvalue,file="MMPvalue")


#-------------筛选枢纽基因---------##
# 计算模块内基因的连接度（kME）
geneModuleMembership <- as.data.frame(cor(Expr3, MET, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# 筛选每个模块中kME值最高的基因
hubGenes <- chooseTopHubInEachModule(Expr3, moduleColors)
write.csv(hubGenes, "Hub_Genes_Each_Module.csv")

# 假设 hubGenes 是之前通过 chooseTopHubInEachModule 得到的结果
hubGenes <- chooseTopHubInEachModule(Expr3, moduleColors)

# 提取所有枢纽基因的表达量（从转置前的原始数据 Expr2 中提取，因为 Expr3 是转置后的）
hub_expr <- Expr2[rownames(Expr2) %in% hubGenes, ]

# 检查提取结果
head(hub_expr)
dim(hub_expr)  # 应为 (枢纽基因数量 × 样本数)

# 添加基因名列并保存
hub_expr_with_gene <- data.frame(
  Gene = rownames(hub_expr),
  hub_expr,
  row.names = NULL  # 移除行名
)

write.csv(
  hub_expr_with_gene,
  file = "Hub_Genes_Expression.csv",
  row.names = FALSE,
  quote = FALSE
)

# 对枢纽基因表达矩阵按行（基因）标准化
hub_expr_scaled <- t(scale(t(hub_expr)))

# 检查标准化后的数据
summary(hub_expr_scaled)

library(pheatmap)

# 设置输出PDF
pdf("Hub_Genes_Expression_Heatmap.pdf", width = 20, height = 40)

# 绘制热图
pheatmap(
  hub_expr_scaled,
  color = colorRampPalette(c("blue", "white", "red"))(50),  # 颜色梯度
  cluster_rows = TRUE,     # 对基因聚类
  cluster_cols = TRUE,     # 对样本聚类
  show_rownames = TRUE,    # 显示基因名
  show_colnames = TRUE,    # 显示样本名
  fontsize_row = 4,        # 基因名字体大小
  fontsize_col = 8,        # 样本名字体大小
  main = "Expression of Hub Genes (Z-score normalized)"
)

# 关闭图形设备
dev.off()
##--------------------------##


#------------获得所有模块热图和表达量————————
module_colors <- unique(moduleColors)  # 包含 "grey"

for (color in module_colors) {
  # 提取当前模块的基因名
  module_genes <- colnames(Expr3)[moduleColors == color]
  
  # 从原始表达矩阵 Expr2 中提取基因表达量
  module_expr <- Expr2[rownames(Expr2) %in% module_genes, ]
  
  # 保存为CSV文件
  write.csv(
    module_expr,
    file = paste0("Module_", color, "_Expression.csv"),
    quote = FALSE
  )
}

library(pheatmap)

for (color in module_colors) {
  # 提取当前模块的基因和表达量
  module_genes <- colnames(Expr3)[moduleColors == color]
  module_expr <- Expr2[rownames(Expr2) %in% module_genes, ]
  
  # 跳过空模块或基因数过少的模块（可选）
  if (nrow(module_expr) < 2) next  # 至少2个基因才能聚类
  
  # 对表达矩阵按行标准化（Z-score）
  expr_scaled <- t(scale(t(module_expr)))
  
  # 设置输出文件名
  pdf_file <- paste0("Module_", color, "_Heatmap.pdf")
  
  # 绘制热图
  pdf(pdf_file, width = 10, height = 8)
  pheatmap(
    expr_scaled,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    clustering_method = "average",
    show_rownames = ifelse(nrow(module_expr) <= 30, TRUE, FALSE),  # 基因少时显示名称
    main = paste("Module:", color, "-", length(module_genes), "genes")
  )
  dev.off()
}

library(ggplot2)
library(reshape2)

# 准备数据：合并所有模块（含grey）的表达矩阵
all_modules_expr <- do.call(rbind, lapply(module_colors, function(color) {
  genes <- colnames(Expr3)[moduleColors == color]
  expr <- Expr2[rownames(Expr2) %in% genes, ]
  if (nrow(expr) > 0) {
    data.frame(
      Gene = rownames(expr),
      Module = color,
      expr
    )
  }
}))

# 转换为长格式并标准化（按基因）
expr_long <- melt(all_modules_expr, id.vars = c("Gene", "Module"))
expr_long$value_scaled <- ave(expr_long$value, expr_long$Gene, FUN = scale)

# 绘制分面热图（按模块分组）
ggplot(expr_long, aes(x = variable, y = Gene, fill = value_scaled)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  facet_wrap(~Module, scales = "free_y", ncol = 3) +  # 自由y轴缩放
  labs(title = "Gene Expression Heatmaps (All Modules, Z-score by Gene)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_blank(),  # 隐藏基因名
    strip.text = element_text(size = 8)  # 调整分面标题字体
  )

ggsave("All_Modules_Including_Grey_Heatmap.pdf", width = 16, height = 20)






