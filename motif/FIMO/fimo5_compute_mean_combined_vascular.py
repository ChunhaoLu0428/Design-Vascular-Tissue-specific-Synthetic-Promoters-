#!/usr/bin/env python3
# nohup python "/home/share_data1/luchh/vascular24ts/STREME/100/fimo_out/enrichment/motif_enrichment_vascula/combined_scores/compute_mean_combined_vascular.py" --dir /home/share_data1/luchh/vascular24ts/STREME/100/fimo_out/enrichment/motif_enrichment_vascula/combined_scores/ --lengths 500 1000 1500 2000 2500 3000 3200 500intron --vascular-genes "/home/share_data1/luchh/vascular24ts/STREME/100/vasculargeneID.txt" --outdir vascula_mean_combined_output &
"""
From combined_score_{length}.tsv produce:
  1. gene×motif mean combined_score across lengths (all genes)
  1b. gene×motif mean combined_score across lengths restricted to vascular genes (if provided and overlapping)
  2. motif×length mean combined_score across genes (optionally restricted to vascular genes)

Includes diagnostics and robust matching of vascular genes (stripping whitespace).
Writes overlapping vascular genes to a file for inspection.
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import logging
import sys

def load_gene_list(path: Path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )

def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="Compute mean combined_score across lengths and motif×length summary.")
    parser.add_argument("--dir", required=True, help="Directory containing combined_score_{length}.tsv files")
    parser.add_argument("--lengths", nargs="+", required=True, help="Length labels, e.g., 500 1000 1500 2000 2500 3000 3200 500intron")
    parser.add_argument("--vascular-genes", help="Optional file of vascular gene IDs to restrict motif×length averaging and produce vascular gene×motif mean")
    parser.add_argument("--outdir", default="mean_combined_results", help="Output directory")
    args = parser.parse_args()

    base = Path(args.dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    vascular_set = None
    if args.vascular_genes:
        vascular_set = load_gene_list(Path(args.vascular_genes))
        vascular_set = {g.strip() for g in vascular_set}
        sample = sorted(list(vascular_set))[:5]
        logging.info("Loaded %d vascular genes (examples: %s ...)", len(vascular_set), sample)

    records = []
    for length in args.lengths:
        f = base / f"combined_score_{length}.tsv"
        if not f.exists():
            logging.warning("Missing %s, skipping length %s", f, length)
            continue
        try:
            df = pd.read_csv(f, sep="\t", index_col=0)
        except Exception as e:
            logging.error("Failed to read %s: %s", f, e)
            continue
        df.index = df.index.astype(str)
        # melt into long form: gene, motif, combined_score
        df_long = df.reset_index().rename(columns={df.index.name or "index": "gene"})
        melted = df_long.melt(id_vars="gene", var_name="motif", value_name="combined_score")
        melted["length"] = length
        records.append(melted)

    if not records:
        logging.error("No combined_score files loaded; check paths/lengths.")
        sys.exit(1)

    long_df = pd.concat(records, ignore_index=True)
    # clean gene names (strip whitespace)
    long_df["gene"] = long_df["gene"].astype(str).str.strip()

    # diagnostic intersection
    all_genes = set(long_df["gene"])
    logging.info("Total distinct genes in combined data: %d", len(all_genes))
    use_vascular_filter = False
    if vascular_set is not None:
        overlap = sorted(all_genes & vascular_set)
        logging.info("Vascular gene list size: %d, intersection with data: %d", len(vascular_set), len(overlap))
        if len(overlap) == 0:
            logging.warning("No vascular genes overlapped; proceeding without vascular restriction for motif×length mean.")
        else:
            use_vascular_filter = True
            # write intersecting vascular genes for manual inspection
            overlap_path = outdir / "vascular_genes_in_data.txt"
            with open(overlap_path, "w") as fo:
                for g in overlap:
                    fo.write(g + "\n")
            logging.info("Wrote %d overlapping vascular genes to %s", len(overlap), overlap_path)

    # 1. gene–motif mean across lengths (all lengths, all genes)
    gm_mean_all = (
        long_df.groupby(["gene", "motif"])["combined_score"]
        .mean()
        .reset_index()
        .pivot(index="gene", columns="motif", values="combined_score")
    )
    gene_motif_out = outdir / "gene_motif_mean_combined.tsv"
    gm_mean_all.to_csv(gene_motif_out, sep="\t", na_rep="NA")
    logging.info("Wrote gene×motif mean combined_score across lengths (ALL genes) to %s", gene_motif_out)

    # 1b. gene–motif mean across lengths restricted to vascular genes (if applicable)
    if use_vascular_filter:
        long_df_vasc = long_df[long_df["gene"].isin(vascular_set)]
        gm_mean_vasc = (
            long_df_vasc.groupby(["gene", "motif"])["combined_score"]
            .mean()
            .reset_index()
            .pivot(index="gene", columns="motif", values="combined_score")
        )
        gene_motif_vasc_out = outdir / "gene_motif_mean_combined_vascular.tsv"
        gm_mean_vasc.to_csv(gene_motif_vasc_out, sep="\t", na_rep="NA")
        logging.info("Wrote gene×motif mean combined_score across lengths (vascular genes only) to %s", gene_motif_vasc_out)

    # 2. motif×length mean across genes (restricted to vascular if available and overlapping)
    df_for_motif_length = long_df
    if use_vascular_filter:
        df_for_motif_length = long_df[long_df["gene"].isin(vascular_set)]
    ml_mean = (
        df_for_motif_length.groupby(["motif", "length"])["combined_score"]
        .mean()
        .reset_index()
        .pivot(index="motif", columns="length", values="combined_score")
    )
    motif_length_out = outdir / "motif_length_mean_combined.tsv"
    ml_mean.to_csv(motif_length_out, sep="\t", na_rep="NA")
    logging.info("Wrote motif×length mean combined_score to %s", motif_length_out)

if __name__ == "__main__":
    main()
