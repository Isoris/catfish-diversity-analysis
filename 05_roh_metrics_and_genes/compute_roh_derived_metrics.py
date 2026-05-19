#!/usr/bin/env python3
"""
compute_roh_derived_metrics.py

Per-sample ROH-derived metrics for the catfish hatchery cohort:

  N_ROH      total / short / medium / long tract counts
  S_ROH      total / short / medium / long ROH bp (also reported as Mb)
  F_ROH      S_ROH / callable_nonrepeat_bp   (joined from STEP_A04 summary)
  F_HOM      cohort-relative inbreeding deficit:
               F_HOM_mean   = 1 - H_obs / mean(H_obs cohort)
               F_HOM_median = 1 - H_obs / median(H_obs cohort)

Inputs
  --roh-bed       chr<TAB>start<TAB>end<TAB>sample<TAB>length
  --roh-summary   STEP_A04 catfish_roh.per_sample_roh.tsv (for F_ROH + denom)
  --het-summary   STEP_A02 genomewide_heterozygosity.tsv (for H_obs)

The H_obs centring used by F_HOM is documented in the methods file and
deliberately exposed in two flavours (mean / median). The HW-expected
denominator that would yield the classical PLINK F_HOM is a follow-up;
note this is a cohort-relative estimator, not the absolute F_HW.
"""

import argparse
import sys

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--roh-bed", required=True)
    p.add_argument("--roh-summary", required=True)
    p.add_argument("--het-summary", required=True)
    p.add_argument("--min-roh-bp", type=int, default=100_000,
                   help="Minimum tract length to count (default 100 kb)")
    p.add_argument("--bin-medium-bp", type=int, default=1_000_000,
                   help="short<medium boundary (default 1 Mb)")
    p.add_argument("--bin-long-bp", type=int, default=5_000_000,
                   help="medium<long boundary (default 5 Mb)")
    p.add_argument("--out-per-sample", required=True)
    p.add_argument("--out-summary", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    # ── Per-sample tract aggregation from the BED ────────────────────────────
    bed = pd.read_csv(args.roh_bed, sep="\t", header=None,
                      names=["chrom", "start", "end", "sample", "length"])
    if bed.empty:
        sys.exit(f"{args.roh_bed} is empty")
    bed["length"] = pd.to_numeric(bed["length"], errors="coerce")
    bed = bed.dropna(subset=["length"])
    bed["length"] = bed["length"].astype(int)
    bed = bed[bed["length"] >= args.min_roh_bp].copy()

    def bin_of(L: int) -> str:
        if L < args.bin_medium_bp:
            return "short"
        if L < args.bin_long_bp:
            return "medium"
        return "long"

    bed["bin"] = bed["length"].map(bin_of)

    grp = bed.groupby("sample")
    counts_all = grp["length"].count().rename("N_ROH_total")
    sums_all = grp["length"].sum().rename("S_ROH_total_bp")

    bin_counts = (bed.groupby(["sample", "bin"])
                  .size().unstack(fill_value=0))
    bin_sums = (bed.groupby(["sample", "bin"])["length"]
                .sum().unstack(fill_value=0))

    for col in ("short", "medium", "long"):
        if col not in bin_counts.columns:
            bin_counts[col] = 0
        if col not in bin_sums.columns:
            bin_sums[col] = 0

    metrics = pd.DataFrame({
        "N_ROH_total": counts_all,
        "N_ROH_short": bin_counts["short"],
        "N_ROH_medium": bin_counts["medium"],
        "N_ROH_long": bin_counts["long"],
        "S_ROH_total_bp": sums_all,
        "S_ROH_short_bp": bin_sums["short"],
        "S_ROH_medium_bp": bin_sums["medium"],
        "S_ROH_long_bp": bin_sums["long"],
    }).fillna(0).astype({
        "N_ROH_total": int, "N_ROH_short": int,
        "N_ROH_medium": int, "N_ROH_long": int,
        "S_ROH_total_bp": int, "S_ROH_short_bp": int,
        "S_ROH_medium_bp": int, "S_ROH_long_bp": int,
    })
    metrics["S_ROH_total_Mb"] = metrics["S_ROH_total_bp"] / 1e6
    metrics["S_ROH_long_Mb"] = metrics["S_ROH_long_bp"] / 1e6
    metrics = metrics.reset_index()

    # ── Join F_ROH + callable denominator from STEP_A04 summary ─────────────
    roh_sum = pd.read_csv(args.roh_summary, sep="\t")
    keep = [c for c in ("sample", "FROH", "callable_nonrepeat_bp",
                        "ancestry_label") if c in roh_sum.columns]
    roh_sum = roh_sum[keep].rename(columns={"FROH": "F_ROH"})
    metrics = metrics.merge(roh_sum, on="sample", how="left")

    # ── Join H_obs from STEP_A02 summary and compute F_HOM ──────────────────
    het = pd.read_csv(args.het_summary, sep="\t")
    if "het_genomewide" not in het.columns:
        sys.exit(f"{args.het_summary}: missing 'het_genomewide' column")
    het = het[["sample", "het_genomewide"]].rename(
        columns={"het_genomewide": "H_obs"})
    het["H_obs"] = pd.to_numeric(het["H_obs"], errors="coerce")
    metrics = metrics.merge(het, on="sample", how="left")

    H = metrics["H_obs"].dropna().values
    H_mean = float(np.mean(H)) if H.size else float("nan")
    H_median = float(np.median(H)) if H.size else float("nan")
    metrics["F_HOM_mean"] = 1.0 - metrics["H_obs"] / H_mean
    metrics["F_HOM_median"] = 1.0 - metrics["H_obs"] / H_median

    # Stable, readable column order
    col_order = [
        "sample", "ancestry_label",
        "H_obs", "F_HOM_mean", "F_HOM_median",
        "N_ROH_total", "N_ROH_short", "N_ROH_medium", "N_ROH_long",
        "S_ROH_total_bp", "S_ROH_short_bp", "S_ROH_medium_bp",
        "S_ROH_long_bp", "S_ROH_total_Mb", "S_ROH_long_Mb",
        "F_ROH", "callable_nonrepeat_bp",
    ]
    col_order = [c for c in col_order if c in metrics.columns]
    metrics = metrics[col_order].sort_values("sample")
    metrics.to_csv(args.out_per_sample, sep="\t", index=False)
    print(f"wrote {args.out_per_sample}  ({len(metrics)} samples)",
          file=sys.stderr)
    print(f"  cohort H_obs mean={H_mean:.4e}  median={H_median:.4e}",
          file=sys.stderr)

    # ── Cohort summary of each metric ───────────────────────────────────────
    summary_cols = ["H_obs", "F_HOM_mean", "F_HOM_median",
                    "N_ROH_total", "N_ROH_long",
                    "S_ROH_total_Mb", "S_ROH_long_Mb", "F_ROH"]
    rows = []
    for c in summary_cols:
        if c not in metrics.columns:
            continue
        v = pd.to_numeric(metrics[c], errors="coerce").dropna().values
        if v.size == 0:
            continue
        q25, med, q75 = np.percentile(v, [25, 50, 75])
        rows.append({
            "metric": c,
            "n": int(v.size),
            "mean": float(np.mean(v)),
            "median": float(med),
            "q25": float(q25),
            "q75": float(q75),
            "iqr": float(q75 - q25),
            "min": float(np.min(v)),
            "max": float(np.max(v)),
        })
    pd.DataFrame(rows).to_csv(args.out_summary, sep="\t", index=False)
    print(f"wrote {args.out_summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
