#!/usr/bin/env python3
"""
compute_window_metrics.py

Two-pass aggregator over per-sample pestPG files:

  Pass 1: read every sample's windowed pestPG, build a long (sample, chrom,
          win_start, n_sites, H_w) table where H_w = tP / nSites.
          Drop windows with nSites < min_callable_frac * win.

  Pass 2: cohort-median H_w per (chrom, win_start) → smoothed median for
          chi_min denominator → per-sample DDI and chi_min.

Outputs:
  --out-per-sample      sample-level summary (DDI, chi_min, ...)
  --out-matrix          gzip-compressed long matrix of kept windows
  --out-cohort-median   per-window cohort median (raw + smoothed)
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd


PESTPG_COLS = [
    "WinInfo", "Chr", "WinCenter",
    "tW", "tP", "tF", "tH", "tL",
    "Tajima", "fuf", "fud", "fayh", "zeng", "nSites",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pestpg-manifest", required=True,
                   help="TSV header sample\\tpestPG_path")
    p.add_argument("--win", type=int, required=True,
                   help="Nominal window size (bp) — used for min-callable cutoff")
    p.add_argument("--min-callable-frac", type=float, default=0.3)
    p.add_argument("--smooth-window", type=int, default=11,
                   help="Odd integer; cohort-median smoothing for chi_min "
                        "denominator (1 disables smoothing)")
    p.add_argument("--mad-constant", type=float, default=1.4826,
                   help="MAD scale constant; the ratio DDI = MAD/median is "
                        "invariant in the constant but reported here for "
                        "transparency")
    p.add_argument("--out-per-sample", required=True)
    p.add_argument("--out-matrix", required=True,
                   help="Long sample x window matrix, gzip TSV")
    p.add_argument("--out-cohort-median", required=True)
    return p.parse_args()


def parse_win_start(wininfo: str):
    """
    pestPG column 1 is "(indexStart,indexStop)(firstPos,lastPos)(WinStart,WinStop)".
    Extract WinStart as int. Returns NaN if parsing fails.
    """
    try:
        last = wininfo.rsplit("(", 1)[-1]
        win_start = last.split(",")[0]
        return int(win_start)
    except Exception:
        return np.nan


def natural_lg_sort_key(chrom: str):
    import re
    m = re.search(r"LG(\d+)", str(chrom))
    if m:
        return (0, int(m.group(1)))
    return (1, str(chrom))


def read_one(sample: str, path: str, min_nsites: int) -> pd.DataFrame:
    if not os.path.exists(path):
        print(f"  WARN: {sample}: pestPG missing — {path}", file=sys.stderr)
        return None
    with open(path) as fh:
        header = fh.readline()
    if not header.startswith("#"):
        print(f"  WARN: {sample}: unexpected pestPG header: {header[:80]!r}",
              file=sys.stderr)
        return None
    df = pd.read_csv(path, sep="\t", skiprows=1, names=PESTPG_COLS)
    df["win_start"] = df["WinInfo"].map(parse_win_start)
    for c in ["tP", "nSites", "win_start"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["tP", "nSites", "win_start"])
    df = df[df["nSites"] >= min_nsites].copy()
    if df.empty:
        return None
    df["sample"] = sample
    df["H_w"] = df["tP"] / df["nSites"]
    df["win_start"] = df["win_start"].astype(int)
    df["nSites"] = df["nSites"].astype(int)
    df = df.rename(columns={"Chr": "chrom"})
    return df[["sample", "chrom", "win_start", "nSites", "H_w"]]


def smoothed_median(series: pd.Series, win: int) -> pd.Series:
    if win <= 1:
        return series
    if win % 2 == 0:
        win += 1
    # min_periods=1 lets edges keep partial windows rather than producing NaN.
    return series.rolling(window=win, center=True, min_periods=1).median()


def main():
    args = parse_args()
    min_nsites = int(round(args.min_callable_frac * args.win))
    print(f"  min nSites per window = {min_nsites} "
          f"({args.min_callable_frac} × {args.win})", file=sys.stderr)

    manifest = pd.read_csv(args.pestpg_manifest, sep="\t")
    if not {"sample", "pestPG_path"}.issubset(manifest.columns):
        sys.exit(f"{args.pestpg_manifest}: expected columns sample, pestPG_path")

    # ── Pass 1: read all samples ────────────────────────────────────────────
    parts = []
    for _, row in manifest.iterrows():
        d = read_one(row["sample"], row["pestPG_path"], min_nsites)
        if d is not None and not d.empty:
            parts.append(d)
    if not parts:
        sys.exit("No usable pestPG data after filtering.")
    long_df = pd.concat(parts, ignore_index=True)
    print(f"  long matrix: {len(long_df):,} (sample, window) rows; "
          f"{long_df['sample'].nunique()} samples, "
          f"{long_df.groupby(['chrom','win_start']).ngroups:,} unique windows",
          file=sys.stderr)

    # Write the gzip-compressed long matrix
    long_df_sorted = long_df.copy()
    long_df_sorted["_key"] = long_df_sorted["chrom"].map(natural_lg_sort_key)
    long_df_sorted = long_df_sorted.sort_values(
        ["sample", "_key", "win_start"]
    ).drop(columns="_key")
    with gzip.open(args.out_matrix, "wt") as fh:
        long_df_sorted.to_csv(fh, sep="\t", index=False)
    print(f"  wrote {args.out_matrix}", file=sys.stderr)

    # ── Pass 2: cohort medians per (chrom, win_start) + smoothing ──────────
    cohort = (long_df
              .groupby(["chrom", "win_start"])["H_w"]
              .agg(["median", "size"])
              .reset_index()
              .rename(columns={"median": "median_H_w_cohort",
                               "size": "n_samples"}))
    cohort["_key"] = cohort["chrom"].map(natural_lg_sort_key)
    cohort = cohort.sort_values(["_key", "win_start"]).drop(columns="_key")

    cohort["median_H_w_smoothed"] = (
        cohort
        .groupby("chrom", group_keys=False)["median_H_w_cohort"]
        .apply(lambda s: smoothed_median(s, args.smooth_window))
        .reset_index(drop=True)
    )
    cohort.to_csv(args.out_cohort_median, sep="\t", index=False)
    print(f"  wrote {args.out_cohort_median}", file=sys.stderr)

    # ── Per-sample DDI + chi_min ───────────────────────────────────────────
    cohort_idx = cohort.set_index(["chrom", "win_start"])[
        ["median_H_w_cohort", "median_H_w_smoothed"]
    ]
    long_df = long_df.join(cohort_idx, on=["chrom", "win_start"])

    # chi_min uses the smoothed cohort median as the denominator
    long_df["chi"] = long_df["H_w"] / long_df["median_H_w_smoothed"]

    out_rows = []
    for sample, g in long_df.groupby("sample"):
        h = g["H_w"].values
        med = float(np.median(h))
        mad_raw = float(np.median(np.abs(h - med)))
        mad_scaled = mad_raw * args.mad_constant
        ddi = mad_scaled / med if med > 0 else np.nan

        # chi_min over windows with a finite smoothed-median denominator
        g_ok = g[np.isfinite(g["chi"])]
        if g_ok.empty:
            chi_min = np.nan
            chi_min_chrom = ""
            chi_min_pos = -1
        else:
            i = g_ok["chi"].idxmin()
            chi_min = float(g_ok.loc[i, "chi"])
            chi_min_chrom = str(g_ok.loc[i, "chrom"])
            chi_min_pos = int(g_ok.loc[i, "win_start"])

        out_rows.append({
            "sample": sample,
            "n_windows_kept": int(len(g)),
            "median_H_w": med,
            "mad_H_w": mad_scaled,
            "DDI": ddi,
            "chi_min": chi_min,
            "chi_min_chrom": chi_min_chrom,
            "chi_min_pos": chi_min_pos,
        })

    summary = pd.DataFrame(out_rows).sort_values("sample")
    summary.to_csv(args.out_per_sample, sep="\t", index=False)
    print(f"  wrote {args.out_per_sample}", file=sys.stderr)


if __name__ == "__main__":
    main()
