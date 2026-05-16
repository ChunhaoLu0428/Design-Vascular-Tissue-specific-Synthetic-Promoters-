#!/usr/bin/env python3

#nohup python "/home/share_data1/luchh/vascular24ts/STREME/100/fimo_out/enrichment/motif_enrichment_vascula/combined_score.py" --dir /home/share_data1/luchh/vascular24ts/STREME/100/fimo_out/enrichment/motif_enrichment_vascula/ --lengths 500 1000 1500 2000 2500 3000 3200 500intron --two-sided-z --outdir /home/share_data1/luchh/vascular24ts/STREME/100/fimo_out/enrichment/motif_enrichment_vascula/combined_scores &
#nohup python "/home/share_data1/luchh/vascular24ts/STREME/enrichment/fimo4_combined_score.py" --dir /home/share_data1/luchh/vascular24ts/STREME/enrichment/099/motif_enrichment_vascular/ --lengths 500 1000 1500 2000 2500 3000 3200 500intron --normalize --outdir /home/share_data1/luchh/vascular24ts/STREME/enrichment/099/combined_scores &


#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import argparse

def compute_combined(z_df: pd.DataFrame, fc_df: pd.DataFrame, two_sided_z: bool, normalize: bool):
    # 转成数值型
    zdf = z_df.apply(pd.to_numeric, errors="coerce")
    fcdf = fc_df.apply(pd.to_numeric, errors="coerce")

    genes = zdf.index
    motifs = [c for c in zdf.columns if c in fcdf.columns]
    combined = pd.DataFrame(index=genes, columns=motifs, dtype=float)

    for motif in motifs:
        zcol = zdf[motif]
        fccol = fcdf[motif]

        # log2FC，保留方向性：>1 正，<1 负
        with np.errstate(divide='ignore', invalid='ignore'):
            log2fc = np.where((fccol > 0) & (~fccol.isna()), np.log2(fccol), np.nan)

        # z component
        if two_sided_z:
            zcomp = zcol.abs()
        else:
            zcomp = zcol  # 保留方向

        if normalize:
            # min-max 标准化到 [-1, 1]
            z_max = np.nanmax(np.abs(zcomp))
            fc_max = np.nanmax(np.abs(log2fc))
            z_scaled = zcomp / z_max if z_max and not np.isnan(z_max) else zcomp
            log2fc_scaled = log2fc / fc_max if fc_max and not np.isnan(fc_max) else log2fc
            combined_score = z_scaled + log2fc_scaled
        else:
            # 直接相加（NaN 的 log2fc 当作 0）
            log2fc_filled = pd.Series(np.where(np.isnan(log2fc), 0.0, log2fc), index=genes)
            combined_score = zcomp + log2fc_filled

        combined.loc[:, motif] = combined_score

    return combined

def main():
    parser = argparse.ArgumentParser(description="Compute gene×motif combined_score = zscore + log2FC per promoter length.")
    parser.add_argument("--dir", required=True, help="Directory containing enrichment_zscore_{length}.tsv and enrichment_foldchange_{length}.tsv")
    parser.add_argument("--lengths", nargs="+", required=True, help="Length labels, e.g., 500 1000 1500 500intron")
    parser.add_argument("--two-sided-z", action="store_true", help="Use absolute z-score when computing combined_score")
    parser.add_argument("--normalize", action="store_true", help="Normalize z-score and log2FC to the same range [-1, 1] before combining (keeps z-score sign unless two-sided-z)")
    parser.add_argument("--outdir", default=None, help="Where to write combined_score_{length}.tsv (defaults to same dir)")
    args = parser.parse_args()

    base = Path(args.dir)
    outbase = Path(args.outdir) if args.outdir else base
    outbase.mkdir(parents=True, exist_ok=True)

    for length in args.lengths:
        zfile = base / f"enrichment_zscore_{length}.tsv"
        fcfile = base / f"enrichment_foldchange_{length}.tsv"
        if not zfile.exists() or not fcfile.exists():
            print(f"[WARN] missing for {length}: {zfile.name if not zfile.exists() else ''} {fcfile.name if not fcfile.exists() else ''}")
            continue

        zdf = pd.read_csv(zfile, sep="\t", index_col=0)
        fcdf = pd.read_csv(fcfile, sep="\t", index_col=0)

        combined_df = compute_combined(zdf, fcdf, two_sided_z=args.two_sided_z, normalize=args.normalize)
        outpath = outbase / f"combined_score_{length}.tsv"
        combined_df.to_csv(outpath, sep="\t", na_rep="NA")
        print(f"Wrote combined_score for {length} to {outpath}")

if __name__ == "__main__":
    main()

