1. fimo1_epresenc_density.py - 核心分析脚本
功能：读取 FIMO 输出，计算每个基因-基序组合的：presence：是否有至少1个结合位点（0/1），density：结合位点数 / 启动子长度
输入：
--fimo-base：FIMO输出目录，结构如 {base}/{length}/fimo.tsv
--special-template：各长度的真实启动子长度文件模板
--lengths：要分析的启动子长度列表
输出（写入 --outdir）：
motif_mean_density_all_lengths.tsv	motif × 长度 的平均密度矩阵
qv_motif_mean_density_all_lengths.tsv	同上，但只包含q-value筛选后的motif
abundance_{length}.tsv	基因×motif密度矩阵（每个长度）
cochran_none/qvalue.tsv	Cochran's Q检验（存在性差异）
density_friedman_allmotifs.tsv	Friedman检验（密度差异）
combined_presence_density_ranking.tsv	综合排名

2. fimo2_motif_mean_density_all_lengths.r - 密度热图
功能：将 motif_mean_density_all_lengths.tsv 绘制为热图
输入：motif_mean_density_all_lengths.tsv
输出：PNG/PDF热图

3. fimo3_enrichment.py - 维管基因富集分析
功能：检测每个motif是否在维管基因中富集
输入：abundance_{length}.tsv（来自 fimo1）
vascular-ids：维管基因ID列表
输出：
enrichment_zscore_{length}.tsv	每个基因的Z-score
enrichment_foldchange_{length}.tsv	每个基因的倍数变化
motif_level_enrichment_all_lengths.tsv	motif级别汇总
统计方法：
Z-score：(vascular_mean - background_mean) / background_std
permutation test：随机抽样背景基因构建零分布
Mann-Whitney U test：vascular vs background

4. fimo4_combined_score.py - 综合得分
功能：将 Z-score 和 log2FC 合并为单一综合得分
公式：
不标准化：combined = zscore + log2(FC)
标准化：combined = z_scaled + log2fc_scaled（范围[-1,1]）
输入：enrichment_zscore_{length}.tsv + enrichment_foldchange_{length}.tsv
输出：combined_score_{length}.tsv

5. fimo5_compute_mean_combined_vascular.py - 汇总平均
功能：跨长度计算平均综合得分
输出：
gene_motif_mean_combined.tsv	基因×motif平均得分（所有基因）
gene_motif_mean_combined_vascular.tsv	同上，仅维管基因
motif_length_mean_combined.tsv	motif×长度平均得分

6. fimo6_motif_length_mean_combined.r - motif×长度热图
输入：motif_length_mean_combined.tsv
输出：热图

7. fimo7_plot_gene_motif_heatmap_vascular.R - 维管基因热图
输入：gene_motif_mean_combined_vascular.tsv
输出：维管基因×motif热图
