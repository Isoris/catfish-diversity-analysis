#!/usr/bin/env python3
"""
adapt_texture_metrics.py

Convert the outputs of 04_window_diversity_texture/ into the
texture_metrics_v1.schema.json payload consumed by the diversity-atlas
`texture` page.

Upstream inputs (from STEP_A06_window_H_and_DDI.sh in module 04):
  --per-sample        per_sample_summary.tsv
                      cols: sample, n_windows_kept, median_H_w, mad_H_w,
                            DDI, chi_min, chi_min_chrom, chi_min_pos
  --cohort-median     cohort_median_per_window.tsv
                      cols: chrom, win_start, median_H_w_cohort,
                            n_samples, median_H_w_smoothed
  --long-matrix       long_window_matrix.tsv.gz
                      cols: sample, chrom, win_start, nSites, H_w
  --het-summary       (optional) genomewide_heterozygosity.tsv
                      cols: sample, H_genomewide, ...   (joined as h_gw)
  --win-bp            nominal window size (bp), default 50000
  --step-bp           step size (bp), default 50000
  --min-callable      min callable sites per window, default 15000
  --out               output JSON path (default stdout)

Emits texture_metrics_v1 with:
  params.{win_bp, step_bp, min_callable}
  cohort_summary.{n_samples, median_ddi, median_h_gw, median_chi_min,
                  mad_h_w_cohort}
  per_sample[]    (one entry per sample)
  windows.{chroms, cohort_median_H_w, cohort_q25_H_w, cohort_q75_H_w}
  per_sample_H_w  keyed by sample_id

When inputs are absent, emits the v0-stub shape so the page renders
'data pending' rather than failing.
"""

import argparse
import gzip
import json
import os
import sys
from datetime import datetime, timezone


def _now_iso():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--per-sample", help="per_sample_summary.tsv from module 04")
    p.add_argument("--cohort-median", help="cohort_median_per_window.tsv")
    p.add_argument("--long-matrix", help="long_window_matrix.tsv(.gz)")
    p.add_argument("--het-summary", help="genomewide_heterozygosity.tsv (optional)")
    p.add_argument("--win-bp", type=int, default=50000)
    p.add_argument("--step-bp", type=int, default=50000)
    p.add_argument("--min-callable", type=int, default=15000)
    p.add_argument("--out", help="output JSON path (default stdout)")
    return p.parse_args()


def _emit_stub(out_path):
    payload = {
        "_doc": ("texture_metrics v0-stub — upstream pipeline outputs not "
                 "available at adapter run time. Re-run "
                 "STEP_A06_window_H_and_DDI.sh, then re-run this adapter."),
        "version": "v0-stub",
        "generated_at": _now_iso(),
        "params": {"win_bp": 50000, "step_bp": 50000, "min_callable": None},
        "cohort_summary": {
            "n_samples": 0,
            "median_ddi": None, "median_h_gw": None,
            "median_chi_min": None, "mad_h_w_cohort": None,
        },
        "per_sample": [],
        "windows": {
            "chroms": [], "cohort_median_H_w": [],
            "cohort_q25_H_w": [], "cohort_q75_H_w": [],
        },
        "per_sample_H_w": {},
    }
    _write(out_path, payload)


def _write(out_path, payload):
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(payload, fh, indent=2)
        print(f"wrote {out_path}", file=sys.stderr)
    else:
        json.dump(payload, sys.stdout, indent=2)


def _open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main():
    import numpy as np
    import pandas as pd

    args = _parse_args()
    have_all = (args.per_sample and os.path.exists(args.per_sample)
                and args.cohort_median and os.path.exists(args.cohort_median)
                and args.long_matrix and os.path.exists(args.long_matrix))
    if not have_all:
        print("adapt_texture_metrics: inputs incomplete — emitting v0-stub",
              file=sys.stderr)
        _emit_stub(args.out)
        return

    per_sample_df = pd.read_csv(args.per_sample, sep="\t")
    cohort_df = pd.read_csv(args.cohort_median, sep="\t")
    with _open_maybe_gz(args.long_matrix) as fh:
        long_df = pd.read_csv(fh, sep="\t")

    # Optional genome-wide H join
    h_gw_by_sample = {}
    if args.het_summary and os.path.exists(args.het_summary):
        het = pd.read_csv(args.het_summary, sep="\t")
        sample_col = next((c for c in ("sample", "sample_id", "ind") if c in het.columns), None)
        h_col = next((c for c in ("H_genomewide", "H_GW", "het", "H") if c in het.columns), None)
        if sample_col and h_col:
            h_gw_by_sample = dict(zip(het[sample_col].astype(str), het[h_col].astype(float)))

    # ---- cohort windows (parallel arrays in chrom order) -------------------
    cohort_df = cohort_df.sort_values(["chrom", "win_start"]).reset_index(drop=True)
    chroms = cohort_df["chrom"].astype(str).tolist()
    cohort_med = cohort_df["median_H_w_cohort"].astype(float).tolist()
    n_win = len(cohort_df)
    # Per-window q25/q75 across samples — compute from the long matrix.
    pivot = (long_df.assign(_win=long_df["chrom"].astype(str) + ":"
                                  + long_df["win_start"].astype(str))
                    .pivot_table(index="_win", columns="sample",
                                 values="H_w", aggfunc="mean"))
    cohort_df["_win"] = cohort_df["chrom"].astype(str) + ":" + cohort_df["win_start"].astype(str)
    pivot = pivot.reindex(cohort_df["_win"].tolist())
    q25 = pivot.quantile(0.25, axis=1, numeric_only=True).fillna(np.nan).astype(float).tolist()
    q75 = pivot.quantile(0.75, axis=1, numeric_only=True).fillna(np.nan).astype(float).tolist()

    # ---- per-sample H_w (parallel to chroms/cohort_med) --------------------
    samples = sorted(long_df["sample"].astype(str).unique().tolist())
    per_sample_H_w = {}
    for s in samples:
        row = pivot.get(s)
        if row is None:
            per_sample_H_w[s] = [None] * n_win
        else:
            per_sample_H_w[s] = [
                None if (v is None or (isinstance(v, float) and not np.isfinite(v))) else float(v)
                for v in row.tolist()
            ]

    # ---- per-sample summary -------------------------------------------------
    per_sample_list = []
    for _, r in per_sample_df.iterrows():
        sid = str(r["sample"])
        per_sample_list.append({
            "sample": sid,
            "h_gw": h_gw_by_sample.get(sid),
            "ddi": _nan_to_none(r.get("DDI")),
            "chi_min": _nan_to_none(r.get("chi_min")),
            "chi_min_chr": _nan_to_none(r.get("chi_min_chrom")),
            "chi_min_pos": _int_or_none(r.get("chi_min_pos")),
            "median_H_w": _nan_to_none(r.get("median_H_w")),
            "mad_H_w": _nan_to_none(r.get("mad_H_w")),
        })

    # ---- cohort summary -----------------------------------------------------
    h_gw_vals = [v for v in h_gw_by_sample.values() if v is not None and np.isfinite(v)]
    ddi_vals = per_sample_df["DDI"].dropna().astype(float)
    chi_vals = per_sample_df["chi_min"].dropna().astype(float)
    mad_vals = per_sample_df["mad_H_w"].dropna().astype(float)
    cohort_summary = {
        "n_samples": int(per_sample_df["sample"].nunique()),
        "median_ddi": float(ddi_vals.median()) if not ddi_vals.empty else None,
        "median_h_gw": float(np.median(h_gw_vals)) if h_gw_vals else None,
        "median_chi_min": float(chi_vals.median()) if not chi_vals.empty else None,
        "mad_h_w_cohort": float(mad_vals.median()) if not mad_vals.empty else None,
    }

    payload = {
        "_doc": ("texture_metrics derived from "
                 "catfish-diversity-analysis/04_window_diversity_texture "
                 "via adapt_texture_metrics.py."),
        "version": "v1-derived",
        "generated_at": _now_iso(),
        "params": {
            "win_bp": int(args.win_bp),
            "step_bp": int(args.step_bp),
            "min_callable": int(args.min_callable),
        },
        "cohort_summary": cohort_summary,
        "per_sample": per_sample_list,
        "windows": {
            "chroms": chroms,
            "cohort_median_H_w": cohort_med,
            "cohort_q25_H_w": q25,
            "cohort_q75_H_w": q75,
        },
        "per_sample_H_w": per_sample_H_w,
    }
    _write(args.out, payload)


def _nan_to_none(v):
    try:
        import numpy as np
        if v is None: return None
        if isinstance(v, str): return v if v else None
        f = float(v)
        return None if not np.isfinite(f) else f
    except Exception:
        return None


def _int_or_none(v):
    try:
        if v is None: return None
        i = int(v)
        return None if i < 0 else i
    except Exception:
        return None


if __name__ == "__main__":
    main()
