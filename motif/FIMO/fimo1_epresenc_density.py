#!/usr/bin/env python3
#python script.py --lengths 500 500intron 1000 1500 2000 2500 3000 3200 --fimo-base /path/to/fimo_root --special-template "/path/to/special_lengths_{length}.tsv" --qvalue 0.05 --outdir results_presence_density --min-genes 3
#nohup python "/home/share_data1/luchh/vascular24ts/STREME/enrichment/fimo1_epresenc_density.py" --lengths 500 500intron 1000 1500 2000 2500 3000 3200 --fimo-base "/home/share_data1/luchh/vascular24ts/STREME/099/099fimo_out/" --special-template "/home/share_data1/luchh/promoter_seq/LM50/{length}/promoter_output/1_5_pro_len.txt" --qvalue 0.05 --outdir results_presence_density --min-genes 3 &
#--lengths 500 500intron 1000 1500 2000 2500 3000 3200
#--fimo-base /home/share_data1/luchh/vascular24ts/STREME/pto_potrihomo/
#--special-template /home/share_data1/luchh/promoter_seq/Potri/{length}_promoter_output/1_5_pro_len.txt
#--qvalue 0.05 --outdir results_presence_density --min-genes 3 &



import argparse
import logging
import re
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2, friedmanchisquare, spearmanr, wilcoxon, norm
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

# ------------------ Config & Utilities ------------------ #

@dataclass
class PipelineConfig:
    lengths: List[str]
    fimo_base: Path
    special_template: str  # expects .format(length=...)
    qvalue: float = 0.05
    outdir: Path = Path("results_presence_density")
    min_genes: int = 3
    log_level: int = logging.INFO

def setup_logging(level: int) -> None:
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )

@lru_cache(maxsize=None)
def parse_default_length(label: str) -> Optional[int]:
    m = re.match(r"^(\d+)", label)
    return int(m.group(1)) if m else None

# ------------------ Data Loading ------------------ #

def load_special_lengths(label: str, template: str) -> Dict[str, float]:
    path = Path(template.format(length=label))
    if not path.exists():
        logging.debug("No special length file for %s at %s", label, path)
        return {}
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, comment="#")
    if df.shape[1] < 2:
        logging.warning("Special lengths file %s format error", path)
        return {}
    df = df.iloc[:, :2]
    df.columns = ["gene", "actual_length"]
    df["actual_length"] = pd.to_numeric(df["actual_length"], errors="coerce")
    mapping = dict(zip(df["gene"], df["actual_length"]))
    logging.info("Loaded %d special lengths from %s", len(mapping), path)
    return mapping

def load_fimo(file_path: Path, filter_mode: str = "qvalue", qvalue_thresh: float = 0.05) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing FIMO file: {file_path}")
    df = pd.read_csv(file_path, sep="\t", comment="#", dtype=str)
    required = {"motif_id", "sequence_name", "q-value", "p-value"}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise ValueError(f"FIMO file {file_path} missing columns: {missing}")
    df["q-value"] = df["q-value"].astype(float)
    df["p-value"] = df["p-value"].astype(float)
    if filter_mode == "qvalue":
        df = df[df["q-value"] <= qvalue_thresh].copy()
    elif filter_mode != "none":
        raise ValueError("filter_mode must be 'none' or 'qvalue'")
    df = df.rename(columns={"motif_id": "motif", "sequence_name": "gene"})
    return df[["gene", "motif"]]

def build_long_format(
    length_labels: List[str],
    fimo_base: Path,
    special_template: str,
    filter_mode: str = "qvalue",
    qvalue_thresh: float = 0.05,
) -> pd.DataFrame:
    all_rows: List[pd.DataFrame] = []
    for label in length_labels:
        fimo_path = fimo_base / label / "fimo.tsv"
        try:
            logging.info("Loading FIMO %s (%s) from %s", label, filter_mode, fimo_path)
            df = load_fimo(fimo_path, filter_mode=filter_mode, qvalue_thresh=qvalue_thresh)
        except Exception as e:
            logging.warning("Skipped %s for %s: %s", label, filter_mode, e)
            continue

        default_len = parse_default_length(label)
        special_map = load_special_lengths(label, special_template)

        summary = df.groupby(["gene", "motif"]).size().reset_index(name="hits")
        summary["length"] = label

        def effective_length(gene: str) -> float:
            if gene in special_map:
                return special_map[gene]
            if default_len is not None:
                return float(default_len)
            return float("nan")

        summary["effective_length"] = summary["gene"].apply(effective_length)
        summary["presence"] = (summary["hits"] >= 1).astype(int)
        summary["density"] = summary.apply(
            lambda row: row["hits"] / row["effective_length"]
            if pd.notna(row["effective_length"]) and row["effective_length"] > 0
            else float("nan"),
            axis=1,
        )
        all_rows.append(summary)

    if not all_rows:
        raise RuntimeError(f"No data loaded for filter_mode={filter_mode}")
    long_df = pd.concat(all_rows, ignore_index=True)
    logging.info(
        "Built long-format for %s: %d motifs, %d genes, lengths: %s",
        filter_mode,
        long_df["motif"].nunique(),
        long_df["gene"].nunique(),
        sorted(long_df["length"].unique()),
    )
    return long_df

# ------------------ Presence Analysis ------------------ #

def presence_matrix_for_motif(long_df: pd.DataFrame, motif: str) -> pd.DataFrame:
    sub = long_df[long_df["motif"] == motif]
    mat = sub.pivot(index="gene", columns="length", values="presence")
    return mat.fillna(0).astype(int)

def cochran_q_test(presence_matrix: pd.DataFrame) -> Tuple[float, float]:
    M = presence_matrix.values
    k = M.shape[1]
    T_i = M.sum(axis=1)
    R_j = M.sum(axis=0)
    T = R_j.sum()
    numerator = (k * (R_j ** 2).sum() - T ** 2)
    denominator = (k * T - (T_i ** 2).sum())
    if denominator == 0:
        return float("nan"), float("nan")
    Q = (k - 1) * numerator / denominator
    p = chi2.sf(Q, df=k - 1)
    return Q, p

def pairwise_mcnemar(presence_matrix: pd.DataFrame) -> pd.DataFrame:
    lengths = list(presence_matrix.columns)
    recs = []
    for a, b in combinations(lengths, 2):
        pa = presence_matrix[a]
        pb = presence_matrix[b]
        b_only = ((pa == 1) & (pb == 0)).sum()
        c_only = ((pa == 0) & (pb == 1)).sum()
        table = [
            [((pa == 1) & (pb == 1)).sum(), b_only],
            [c_only, ((pa == 0) & (pb == 0)).sum()],
        ]
        try:
            res = mcnemar(table, exact=True)
        except Exception:
            res = mcnemar(table, exact=False)
        oratio = (b_only + 0.5) / (c_only + 0.5)
        log2_or = np.log2(oratio)
        recs.append({
            "len_a": a, "len_b": b,
            "log2_oddsratio": log2_or,
            "pvalue": res.pvalue,
            "b_only": b_only, "c_only": c_only,
        })
    df = pd.DataFrame(recs)
    if not df.empty:
        df["p_adj"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    return df

def summarize_presence(long_df: pd.DataFrame, min_genes: int = 3) -> Tuple[pd.DataFrame, pd.DataFrame]:
    motifs = sorted(long_df["motif"].unique())
    cochran_records, pairwise_frames = [], []
    tested = skipped = 0
    for motif in motifs:
        mat = presence_matrix_for_motif(long_df, motif)
        if mat.shape[0] < min_genes or mat.shape[1] < 2:
            skipped += 1
            continue
        tested += 1
        Q, p = cochran_q_test(mat)
        cochran_records.append({"motif": motif, "Q": Q, "pvalue": p})
        pw = pairwise_mcnemar(mat)
        if not pw.empty:
            pw["motif"] = motif
            pairwise_frames.append(pw)
    logging.info("Presence: tested %d motifs, skipped %d", tested, skipped)
    cochran_df = pd.DataFrame(cochran_records)
    if not cochran_df.empty:
        cochran_df["p_adj"] = multipletests(cochran_df["pvalue"].fillna(1), method="fdr_bh")[1]
    pairwise_df = pd.concat(pairwise_frames, ignore_index=True) if pairwise_frames else pd.DataFrame()
    return cochran_df, pairwise_df

# ------------------ Density Analysis ------------------ #

def density_matrix_for_motif(
    long_df: pd.DataFrame, motif: str, lengths_ordered: List[str]
) -> pd.DataFrame:
    sub = long_df[long_df["motif"] == motif]
    pdens = sub.pivot(index="gene", columns="length", values="density")
    peff = sub.pivot(index="gene", columns="length", values="effective_length")
    genes = pdens.index.union(peff.index)
    pdens = pdens.reindex(index=genes, columns=lengths_ordered)
    peff = peff.reindex(index=genes, columns=lengths_ordered)
    invalid = peff.isna() | (peff <= 0)
    filled = pdens.fillna(0.0).where(~invalid, other=np.nan)
    return filled.astype(float)

def friedman_density_test(dmat: pd.DataFrame, min_genes: int = 3) -> Tuple[float, float]:
    df = dmat.dropna(axis=0, how="any")
    if df.shape[0] < min_genes or df.shape[1] < 2:
        return float("nan"), float("nan")
    try:
        args = [df[col].values for col in df.columns]
        stat, p = friedmanchisquare(*args)
        return stat, p
    except Exception as e:
        logging.warning("Friedman failed: %s", e)
        return float("nan"), float("nan")

def pairwise_wilcoxon_density(dmat: pd.DataFrame) -> pd.DataFrame:
    lengths = list(dmat.columns)
    recs = []
    for a, b in combinations(lengths, 2):
        da, db = dmat[a], dmat[b]
        mask = (~da.isna()) & (~db.isna())
        if mask.sum() < 3:
            continue
        try:
            _, pval = wilcoxon(da[mask], db[mask], alternative="two-sided", zero_method="wilcox")
        except:
            continue
        eps = 1e-9
        log2fc = np.log2((da[mask].mean() + eps)/(db[mask].mean() + eps))
        recs.append({"len_a": a, "len_b": b, "log2_mean_ratio": log2fc, "pvalue": pval})
    df = pd.DataFrame(recs)
    if not df.empty:
        df["p_adj"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    return df

# ------------------ Mean density table for R heatmap ------------------ #

def build_mean_density_matrix(long_df: pd.DataFrame, lengths_ordered: List[str]) -> pd.DataFrame:
    """
    motif x promoter-type mean density (average over genes).
    """
    motifs = sorted(long_df["motif"].unique())
    records = []
    for mot in motifs:
        dmat = density_matrix_for_motif(long_df, mot, lengths_ordered)  # gene x length
        mean_per_length = dmat.mean(axis=0, skipna=True)
        row = mean_per_length.to_dict()
        row["motif"] = mot
        records.append(row)
    df = pd.DataFrame(records).set_index("motif")
    df = df.reindex(columns=lengths_ordered)
    return df

# ------------------ Combined ranking ------------------ #

def build_combined_ranking(cochran_df: pd.DataFrame, friedman_df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    """
    Combine presence (Cochran’s Q adjusted p) and density (Friedman adjusted p) signals.
    Uses Stouffer's method and rank product to rank motifs.
    """
    df = pd.DataFrame({"motif": sorted(set(cochran_df["motif"]).union(friedman_df["motif"]))})
    df = df.merge(cochran_df[["motif", "p_adj"]].rename(columns={"p_adj": "presence_padj"}), on="motif", how="left")
    df = df.merge(friedman_df[["motif", "p_adj"]].rename(columns={"p_adj": "density_padj"}), on="motif", how="left")
    df["presence_padj"] = df["presence_padj"].fillna(1.0)
    df["density_padj"] = df["density_padj"].fillna(1.0)

    eps = 1e-300
    df["z_presence"] = norm.ppf(1 - df["presence_padj"].clip(eps, 1 - eps))
    df["z_density"] = norm.ppf(1 - df["density_padj"].clip(eps, 1 - eps))

    df["combined_z"] = (df["z_presence"] + df["z_density"]) / np.sqrt(2)
    df["combined_p_stouffer"] = 1 - norm.cdf(df["combined_z"])

    df["rank_presence"] = df["presence_padj"].rank(method="average", ascending=True)
    df["rank_density"] = df["density_padj"].rank(method="average", ascending=True)
    df["rank_product"] = np.sqrt(df["rank_presence"] * df["rank_density"])

    df = df.sort_values(["combined_p_stouffer", "rank_product"])
    return df.head(top_n)

# ------------------ Plotting (trend + cochran compare) ------------------ #

def plot_density_trend(
    motif: str,
    density_matrix: pd.DataFrame,
    friedman_p_adj: float,
    lengths_ordered: List[str],
    outdir: Path,
) -> None:
    mean = density_matrix.mean(axis=0, skipna=True)
    sem = density_matrix.sem(axis=0, skipna=True)
    xs = list(range(len(lengths_ordered)))
    plt.figure(figsize=(6,4))
    plt.errorbar(xs, mean[lengths_ordered], yerr=sem[lengths_ordered], marker="o", linestyle="-")
    plt.xticks(xs, lengths_ordered, rotation=45)
    plt.xlabel("Promoter type / length")
    plt.ylabel("Mean density")
    title = f"{motif} density trend"
    if not np.isnan(friedman_p_adj):
        title += f" (Friedman adj p={friedman_p_adj:.2e})"
    plt.title(title)
    plt.tight_layout()
    fname = outdir / f"density_trend_{motif.replace('/', '_')}.png"
    plt.savefig(fname)
    plt.close()
    logging.debug("Saved density trend plot for %s to %s", motif, fname)

def compare_cochran_scatter(
    merged_df: pd.DataFrame, strategy_a: str, strategy_b: str, outdir: Path
) -> None:
    def neglog10(s: pd.Series) -> np.ndarray:
        arr = np.array(s, dtype=float)
        with np.errstate(divide="ignore"):
            nl = -np.log10(arr)
        nl[~np.isfinite(nl)] = np.nan
        return nl
    xa = neglog10(merged_df[f"pvalue_{strategy_a}"])
    yb = neglog10(merged_df[f"pvalue_{strategy_b}"])
    mask = ~np.isnan(xa) & ~np.isnan(yb)
    rho = np.nan
    if mask.sum() >= 3:
        rho, _ = spearmanr(xa[mask], yb[mask])
    plt.figure()
    plt.scatter(xa, yb, alpha=0.5)
    lim = np.nanmax([np.nanmax(xa), np.nanmax(yb)])
    if not np.isnan(lim):
        plt.plot([0, lim], [0, lim], "--", color="gray")
    plt.xlabel(f"-log10 pvalue {strategy_a}")
    plt.ylabel(f"-log10 pvalue {strategy_b}")
    plt.title(f"Cochran Q: {strategy_a} vs {strategy_b} (rho={rho:.2f})")
    plt.tight_layout()
    outpath = outdir / f"cochran_compare_{strategy_a}_vs_{strategy_b}.png"
    plt.savefig(outpath)
    plt.close()
    logging.info("Saved Cochran comparison plot to %s", outpath)

# ------------------ Orchestration ------------------ #

def parse_args() -> PipelineConfig:
    parser = argparse.ArgumentParser(description="Compare motif presence/density across promoter lengths.")
    parser.add_argument("--lengths", nargs="+", required=True)
    parser.add_argument("--fimo-base", required=True)
    parser.add_argument("--special-template", required=True)
    parser.add_argument("--qvalue", type=float, default=0.05)
    parser.add_argument("--outdir", default="results_presence_density")
    parser.add_argument("--min-genes", type=int, default=3)
    args = parser.parse_args()
    return PipelineConfig(
        lengths=args.lengths,
        fimo_base=Path(args.fimo_base),
        special_template=args.special_template,
        qvalue=args.qvalue,
        outdir=Path(args.outdir),
        min_genes=args.min_genes,
    )

def main() -> None:
    config = parse_args()
    setup_logging(config.log_level)
    logging.info("Pipeline configuration: %s", config)

    config.outdir.mkdir(parents=True, exist_ok=True)

    # Build long-format data
    long_none = build_long_format(
        config.lengths, config.fimo_base, config.special_template, filter_mode="none", qvalue_thresh=config.qvalue
    )
    long_qv = build_long_format(
        config.lengths, config.fimo_base, config.special_template, filter_mode="qvalue", qvalue_thresh=config.qvalue
    )

    # === 输出 motif x promoter type 的 mean density 表（所有 motif，用 long_none） ===
    mean_density_all = build_mean_density_matrix(long_none, config.lengths)
    mean_density_path = config.outdir / "motif_mean_density_all_lengths.tsv"
    mean_density_all.to_csv(mean_density_path, sep="\t", na_rep="NA")
    logging.info("Wrote mean density table for all motifs to %s", mean_density_path)
    
    # mean density for qvalue-filtered motifs
    mean_density_qv = build_mean_density_matrix(long_qv, config.lengths)
    mean_density_qv.to_csv(config.outdir / "qv_motif_mean_density_all_lengths.tsv", sep="\t", na_rep="NA")
    logging.info("Wrote mean density tables to %s", config.outdir)

    # === abundance per length (all motifs and qvalue filtered) ===
    for length in config.lengths:
        sub = long_none[long_none["length"] == length]
        if sub.empty:
            logging.warning("No data for length %s when writing abundance matrix", length)
            continue
        abundance_mat = sub.pivot(index="gene", columns="motif", values="density")
        outpath = config.outdir / f"abundance_{length}.tsv"
        abundance_mat.to_csv(outpath, sep="\t", na_rep="NA")
        logging.info("Wrote abundance matrix (all motifs) for length %s to %s", length, outpath)

    for length in config.lengths:
        sub_qv = long_qv[long_qv["length"] == length]
        if sub_qv.empty:
            continue
        abundance_qv = sub_qv.pivot(index="gene", columns="motif", values="density")
        outpath_qv = config.outdir / f"abundance_qvalue_{length}.tsv"
        abundance_qv.to_csv(outpath_qv, sep="\t", na_rep="NA")
        logging.info("Wrote abundance matrix (qvalue-filtered motifs) for length %s to %s", length, outpath_qv)

    # Presence analysis
    coch_none, pw_none = summarize_presence(long_none, min_genes=config.min_genes)
    coch_qv, pw_qv = summarize_presence(long_qv, min_genes=config.min_genes)

    coch_none.to_csv(config.outdir / "cochran_none.tsv", sep="\t", index=False)
    coch_qv.to_csv(config.outdir / "cochran_qvalue.tsv", sep="\t", index=False)
    pw_none.to_csv(config.outdir / "pairwise_none.tsv", sep="\t", index=False)
    pw_qv.to_csv(config.outdir / "pairwise_qvalue.tsv", sep="\t", index=False)
    logging.info("Wrote presence tables.")

    # Merge Cochran comparison
    all_motifs = sorted(set(coch_none["motif"]).union(coch_qv["motif"]))
    merged = []
    for mot in all_motifs:
        row = {"motif": mot}
        for label, df in [("none", coch_none), ("qvalue", coch_qv)]:
            sub = df[df["motif"] == mot]
            if not sub.empty:
                row[f"Q_{label}"] = sub.iloc[0]["Q"]
                row[f"pvalue_{label}"] = sub.iloc[0]["pvalue"]
                row[f"p_adj_{label}"] = sub.iloc[0]["p_adj"]
            else:
                row[f"Q_{label}"] = row[f"pvalue_{label}"] = row[f"p_adj_{label}"] = np.nan
        merged.append(row)
    merged_df = pd.DataFrame(merged)
    merged_df.to_csv(config.outdir / "merged_cochran_comparison.tsv", sep="\t", index=False)
    logging.info("Wrote merged Cochran comparison.")

    compare_cochran_scatter(merged_df, "none", "qvalue", config.outdir)

    # Density analysis using all motifs (long_none)
    df_fr_all, density_pw_all = [], []
    motifs_all = sorted(long_none["motif"].unique())
    logging.info("Starting density analysis (ALL motifs) on %d motifs", len(motifs_all))
    for mot in motifs_all:
        dmat = density_matrix_for_motif(long_none, mot, config.lengths)
        stat, p = friedman_density_test(dmat, min_genes=config.min_genes)
        df_fr_all.append({"motif": mot, "stat": stat, "pvalue": p})
        pw = pairwise_wilcoxon_density(dmat)
        if not pw.empty:
            pw["motif"] = mot
            density_pw_all.append(pw)
    df_fr_all = pd.DataFrame(df_fr_all)
    if not df_fr_all.empty:
        df_fr_all["p_adj"] = multipletests(df_fr_all["pvalue"].fillna(1), method="fdr_bh")[1]
    df_fr_all.to_csv(config.outdir / "density_friedman_allmotifs.tsv", sep="\t", index=False)
    if density_pw_all:
        pd.concat(density_pw_all, ignore_index=True).to_csv(
            config.outdir / "density_pairwise_wilcoxon_allmotifs.tsv", sep="\t", index=False
        )
    logging.info("Wrote density tables for all motifs.")

    # (Optional) density analysis on qvalue-filtered motifs for comparison
    df_fr_qv, density_pw_qv = [], []
    motifs_qv = sorted(long_qv["motif"].unique())
    logging.info("Starting density analysis (qvalue-filtered motifs) on %d motifs", len(motifs_qv))
    for mot in motifs_qv:
        dmat = density_matrix_for_motif(long_qv, mot, config.lengths)
        stat, p = friedman_density_test(dmat, min_genes=config.min_genes)
        df_fr_qv.append({"motif": mot, "stat": stat, "pvalue": p})
        pw = pairwise_wilcoxon_density(dmat)
        if not pw.empty:
            pw["motif"] = mot
            density_pw_qv.append(pw)
    df_fr_qv = pd.DataFrame(df_fr_qv)
    if not df_fr_qv.empty:
        df_fr_qv["p_adj"] = multipletests(df_fr_qv["pvalue"].fillna(1), method="fdr_bh")[1]
    df_fr_qv.to_csv(config.outdir / "density_friedman_qvalue.tsv", sep="\t", index=False)
    if density_pw_qv:
        pd.concat(density_pw_qv, ignore_index=True).to_csv(
            config.outdir / "density_pairwise_wilcoxon_qvalue.tsv", sep="\t", index=False
        )
    logging.info("Wrote density tables for qvalue-filtered motifs.")

    # === combined ranking presence + density ===
    combined_df = build_combined_ranking(coch_qv, df_fr_all, top_n=100)
    combined_out = config.outdir / "combined_presence_density_ranking.tsv"
    combined_df.to_csv(combined_out, sep="\t", index=False)
    logging.info("Wrote combined presence+density ranking to %s", combined_out)

    # Plot trends for all-motifs density (主做 all motifs)
    for mot in motifs_all:
        dmat = density_matrix_for_motif(long_none, mot, config.lengths)
        p_adj = df_fr_all.loc[df_fr_all["motif"] == mot, "p_adj"].squeeze() if "p_adj" in df_fr_all else np.nan
        plot_density_trend(mot, dmat, p_adj, config.lengths, config.outdir)

    logging.info("Pipeline complete. Results in %s", config.outdir)

if __name__ == "__main__":
    main()
