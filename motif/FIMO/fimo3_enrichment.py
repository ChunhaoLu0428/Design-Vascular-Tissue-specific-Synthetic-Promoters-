#!/usr/bin/env python3
#nohup python "/home/share_data1/luchh/vascular24ts/STREME/100/fimo_out/enrichment/motif_enrichment_vascular.py"  --abundance-dir /home/share_data1/luchh/vascular24ts/STREME/100/fimo_out/results_presence_density/results_presence_density --vascular-ids "/home/share_data1/luchh/vascular24ts/STREME/100/vasculargeneID.txt" --lengths 500 1000 1500 2000 2500 3000 3200 500intron --permutations 1000 --out motif_enrichment_vascular.tsv &
#nohup python "/home/share_data1/luchh/vascular24ts/STREME/enrichment/fimo3_enrichment.py" --abundance-dir /home/share_data1/luchh/vascular24ts/STREME/enrichment/099/results_presence_density/ --vascular-ids "/home/share_data1/luchh/vascular24ts/STREME/enrichment/099/099vasculargeneID.txt" --lengths 500 1000 1500 2000 2500 3000 3200 500intron --permutations 1000 --out 099motif_enrichment_vascular &


import argparse
from pathlib import Path
import logging
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# ------------------ Utilities ------------------ #

def setup_logging():
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")

def read_vascular_ids(path: Path) -> set:
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())

def gene_motif_enrichment_tables(abundance_df: pd.DataFrame, vascular_ids: set, pseudocount: float = 1e-9):
    """
    输入某 length 的 abundance gene x motif matrix（NA 已处理为 0）。
    Background = 非 vascular genes。
    返回两个 DataFrame：z-score 和 fold-change (gene x motif)。
    """
    genes = abundance_df.index.astype(str)
    motifs = abundance_df.columns.tolist()

    vascular_genes = [g for g in genes if g in vascular_ids]
    background_genes = [g for g in genes if g not in vascular_ids]

    if len(background_genes) == 0:
        raise ValueError("No background (non-vascular) genes available to build null model.")

    # ensure no NaN, treat absence as zero
    abundance_df = abundance_df.fillna(0.0)
    bg_mat = abundance_df.loc[background_genes].astype(float)

    mean_bg = bg_mat.mean(axis=0, skipna=True)
    std_bg = bg_mat.std(axis=0, ddof=1, skipna=True)
    std_bg_safe = std_bg.replace(0, np.nan)

    zscore_mat = pd.DataFrame(index=genes, columns=motifs, dtype=float)
    foldchange_mat = pd.DataFrame(index=genes, columns=motifs, dtype=float)

    for motif in motifs:
        m_bg = mean_bg.get(motif, np.nan)
        s_bg = std_bg_safe.get(motif, np.nan)
        vals = abundance_df[motif].astype(float)

        foldchange_mat[motif] = (vals + pseudocount) / (m_bg + pseudocount)

        if pd.isna(s_bg):
            zscore_mat[motif] = np.nan
        else:
            zscore_mat[motif] = (vals - m_bg) / s_bg

    return zscore_mat, foldchange_mat

def analyze_motif_level(abundance_df: pd.DataFrame, vascular_ids: set, perms: int = 1000, pseudocount: float = 1e-9):
    """
    motif-level enrichment: vascular genes vs background per motif.
    返回 DataFrame 每 motif 一行，含 mean_vascular, mean_background, fold_change,
    empirical z, empirical p, Mann-Whitney p.
    """
    genes = abundance_df.index.astype(str)
    motifs = abundance_df.columns.tolist()

    vascular_genes = [g for g in genes if g in vascular_ids]
    background_genes = [g for g in genes if g not in vascular_ids]

    if len(vascular_genes) == 0:
        logging.warning("No vascular genes present for this length.")
        return pd.DataFrame()

    # treat missing as zero
    abundance_df = abundance_df.fillna(0.0)
    bg_mat = abundance_df.loc[background_genes].astype(float) if len(background_genes) > 0 else pd.DataFrame()

    results = []

    for motif in motifs:
        v_vals = abundance_df.loc[[g for g in vascular_genes if g in abundance_df.index], motif].astype(float).to_numpy()
        b_vals = bg_mat[motif].astype(float).to_numpy() if len(background_genes) > 0 else np.array([])

        mean_v = np.nanmean(v_vals) if v_vals.size > 0 else np.nan
        mean_b = np.nanmean(b_vals) if b_vals.size > 0 else np.nan

        fc = (mean_v + pseudocount) / (mean_b + pseudocount)

        # empirical null via permutation sampling from background
        z_emp = np.nan
        p_emp = np.nan
        n_v = len(v_vals)
        if b_vals.size > 0 and n_v > 0:
            replace = False
            if len(background_genes) < n_v:
                replace = True
                logging.warning(
                    "Background size (%d) < vascular size (%d); doing permutation with replacement for motif-level null.",
                    len(background_genes),
                    n_v,
                )
            null_means = []
            for _ in range(perms):
                sampled = np.random.choice(b_vals, size=n_v, replace=replace)
                null_means.append(np.nanmean(sampled))
            null_means = np.array(null_means)
            null_mean = np.nanmean(null_means)
            null_std = np.nanstd(null_means, ddof=1)
            if null_std > 0:
                z_emp = (mean_v - null_mean) / null_std
            p_emp = (np.sum(null_means >= mean_v) + 1) / (len(null_means) + 1)
        # Mann-Whitney U test (vascular > background)
        try:
            if b_vals.size > 0 and len(v_vals) > 0:
                _, p_mwu = mannwhitneyu(v_vals, b_vals, alternative="greater")
            else:
                p_mwu = np.nan
        except Exception:
            p_mwu = np.nan

        results.append({
            "motif": motif,
            "mean_vascular": mean_v,
            "mean_background": mean_b,
            "fold_change": fc,
            "z_empirical": z_emp,
            "p_empirical": p_emp,
            "p_mwu": p_mwu,
        })

    return pd.DataFrame(results)

# ------------------ Main ------------------ #

def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="Vascular motif enrichment: gene×motif and motif-level (fixed NaN handling).")
    parser.add_argument("--abundance-dir", required=True, help="Directory containing abundance_{length}.tsv")
    parser.add_argument("--vascular-ids", required=True, help="File with vascular gene IDs, one per line")
    parser.add_argument("--lengths", nargs="+", required=True, help="Promoter length labels, e.g., 500 1000 1500 500intron")
    parser.add_argument("--permutations", type=int, default=1000, help="Permutations for empirical null")
    parser.add_argument("--pseudocount", type=float, default=1e-9, help="Small constant to avoid zeros in fold change")
    parser.add_argument("--outdir", default="motif_enrichment_vascular_fixed", help="Output directory")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()

    np.random.seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    vascular_ids = read_vascular_ids(Path(args.vascular_ids))
    all_motif_level = []

    for length in args.lengths:
        abundance_path = Path(args.abundance_dir) / f"abundance_{length}.tsv"
        if not abundance_path.exists():
            logging.warning("Missing abundance file for length %s: %s", length, abundance_path)
            continue

        logging.info("Processing promoter length: %s", length)
        df_abun = pd.read_csv(abundance_path, sep="\t", index_col=0)
        df_abun.index = df_abun.index.astype(str)
        df_abun = df_abun.fillna(0.0)  # 修复：缺失当作 0

        # gene×motif enrichment (z-score and fold-change)
        try:
            zscore_df, fc_df = gene_motif_enrichment_tables(df_abun, vascular_ids, pseudocount=args.pseudocount)
        except ValueError as e:
            logging.error("Skipping gene×motif enrichment for %s: %s", length, e)
            continue

        zscore_out = outdir / f"enrichment_zscore_{length}.tsv"
        fc_out = outdir / f"enrichment_foldchange_{length}.tsv"
        zscore_df.to_csv(zscore_out, sep="\t", na_rep="NA")
        fc_df.to_csv(fc_out, sep="\t", na_rep="NA")
        logging.info("Wrote gene×motif enrichment (z-score) to %s", zscore_out)
        logging.info("Wrote gene×motif enrichment (fold-change) to %s", fc_out)

        # motif-level enrichment
        motif_df = analyze_motif_level(df_abun, vascular_ids, perms=args.permutations, pseudocount=args.pseudocount)
        if motif_df.empty:
            continue
        motif_df["length"] = length
        all_motif_level.append(motif_df)

    if not all_motif_level:
        logging.error("No motif-level data was generated; exiting.")
        return

    combined = pd.concat(all_motif_level, ignore_index=True)

    # multiple testing correction across all motif×length for empirical and mwu p-values
    for col in ["p_empirical", "p_mwu"]:
        if col in combined.columns:
            mask = combined[col].notna()
            if mask.sum() > 0:
                adj = multipletests(combined.loc[mask, col], method="fdr_bh")[1]
                combined.loc[mask, f"{col}_adj"] = adj
            else:
                combined[f"{col}_adj"] = np.nan
        else:
            combined[f"{col}_adj"] = np.nan

    combined["log2FC"] = np.log2(combined["fold_change"].replace(0, np.nan))
    combined = combined.sort_values(
        by=["p_empirical_adj", "p_mwu_adj", "log2FC"],
        ascending=[True, True, False],
        na_position="last",
    )

    motif_out = outdir / "motif_level_enrichment_all_lengths.tsv"
    combined.to_csv(motif_out, sep="\t", index=False)
    logging.info("Wrote motif-level enrichment to %s", motif_out)

if __name__ == "__main__":
    main()
