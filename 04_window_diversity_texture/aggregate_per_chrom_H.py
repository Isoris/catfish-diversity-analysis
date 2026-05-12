#!/usr/bin/env python3
"""
aggregate_per_chrom_H.py

Consolidates per-sample per-chromosome theta_pi sums (one TSV per sample,
produced by `thetaStat print | awk`) into:

  - long-format per-sample table with H_chr = sum_theta_pi / n_callable_chr
  - wide cohort matrix (rows = sample, cols = chrom)
  - per-chromosome cohort summary (median, IQR, range, n_outliers @ 1.5*IQR)

Inputs are TSVs with columns: sample, chrom, n_sites_obs, sum_theta_pi.
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-sample-dir", required=True,
                   help="Directory of <sample>.per_chrom_tP.tsv files")
    p.add_argument("--callable-cache", required=True,
                   help="TSV (chrom, n_callable) for the cohort callable mask")
    p.add_argument("--out-long", required=True)
    p.add_argument("--out-wide", required=True)
    p.add_argument("--out-summary", required=True)
    return p.parse_args()


def natural_lg_sort_key(chrom: str):
    # Sort C_gar_LG01 < ... < C_gar_LG28 < anything else (alphabetical)
    import re
    m = re.search(r"LG(\d+)", str(chrom))
    if m:
        return (0, int(m.group(1)))
    return (1, str(chrom))


def main():
    args = parse_args()

    callable_df = pd.read_csv(args.callable_cache, sep="\t")
    if not {"chrom", "n_callable"}.issubset(callable_df.columns):
        sys.exit(f"{args.callable_cache}: expected columns chrom, n_callable")
    callable_map = dict(zip(callable_df["chrom"], callable_df["n_callable"]))

    files = sorted(glob.glob(os.path.join(args.per_sample_dir,
                                          "*.per_chrom_tP.tsv")))
    if not files:
        sys.exit(f"No per-sample files in {args.per_sample_dir}")

    rows = []
    for fp in files:
        df = pd.read_csv(fp, sep="\t")
        if df.empty:
            continue
        df["n_callable_chr"] = df["chrom"].map(callable_map)
        missing = df[df["n_callable_chr"].isna()]
        if not missing.empty:
            for c in sorted(missing["chrom"].unique()):
                print(f"  WARNING: {os.path.basename(fp)}: chrom {c} not in "
                      f"callable cache; H_chr will be NaN", file=sys.stderr)
        df["H_chr"] = df["sum_theta_pi"] / df["n_callable_chr"]
        rows.append(df)

    if not rows:
        sys.exit("No non-empty per-sample tables found.")

    long_df = pd.concat(rows, ignore_index=True)
    long_df = long_df[["sample", "chrom", "n_sites_obs",
                       "n_callable_chr", "sum_theta_pi", "H_chr"]]

    # Sort long output by sample then chrom (natural LG order)
    long_df["_key"] = long_df["chrom"].map(natural_lg_sort_key)
    long_df = long_df.sort_values(["sample", "_key"]).drop(columns="_key")
    long_df.to_csv(args.out_long, sep="\t", index=False)
    print(f"wrote {args.out_long}  ({len(long_df):,} rows)", file=sys.stderr)

    # Wide matrix: sample x chrom
    wide = long_df.pivot(index="sample", columns="chrom", values="H_chr")
    wide = wide[sorted(wide.columns, key=natural_lg_sort_key)]
    wide.to_csv(args.out_wide, sep="\t")
    print(f"wrote {args.out_wide}  ({wide.shape[0]} samples x "
          f"{wide.shape[1]} chroms)", file=sys.stderr)

    # Per-chromosome cohort summary
    summary_rows = []
    for chrom in wide.columns:
        v = wide[chrom].dropna().values
        if v.size == 0:
            continue
        q25, med, q75 = np.percentile(v, [25, 50, 75])
        iqr = q75 - q25
        lo = q25 - 1.5 * iqr
        hi = q75 + 1.5 * iqr
        n_out_low = int((v < lo).sum())
        n_out_high = int((v > hi).sum())
        summary_rows.append({
            "chrom": chrom,
            "n_samples": int(v.size),
            "median_H_chr": float(med),
            "q25_H_chr": float(q25),
            "q75_H_chr": float(q75),
            "iqr_H_chr": float(iqr),
            "min_H_chr": float(v.min()),
            "max_H_chr": float(v.max()),
            "n_outliers_low": n_out_low,
            "n_outliers_high": n_out_high,
        })
    summ = pd.DataFrame(summary_rows)
    summ["_key"] = summ["chrom"].map(natural_lg_sort_key)
    summ = summ.sort_values("_key").drop(columns="_key")
    summ.to_csv(args.out_summary, sep="\t", index=False)
    print(f"wrote {args.out_summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
