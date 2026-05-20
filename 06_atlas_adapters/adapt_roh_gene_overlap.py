#!/usr/bin/env python3
"""
adapt_roh_gene_overlap.py

Convert the outputs of 05_roh_metrics_and_genes/ into the
roh_gene_overlap_v1.schema.json payload consumed by the diversity-atlas
`roh` page extension (Plot A cumulative gene-burden + Plot B biotype ×
peak heatmap-table).

Upstream inputs:
  --roh-bed                 per-sample ROH BED
                            cols: chrom, start, end, sample, length
                            (from STEP_A04_parse_roh_and_het.sh)
  --genes-bed               gene BED produced by `genes_in_roh.py gff-to-bed`
                            cols: chrom, start, end, gene_id, gene_name, strand
                            (optional 7th col: biotype — if present, used for
                            peaks.biotype_counts; otherwise everything counts as
                            'protein_coding'.)
  --cohort-recurrence       genes_in_roh_cohort_recurrence.tsv
                            cols: gene_id, gene_name, gene_chr, gene_start,
                                  gene_end, strand, n_samples_with_overlap, ...
                            (from `genes_in_roh.py aggregate`)
  --k8-assign               (optional) per-sample K=8 cluster assignment TSV
                            cols: sample, k8_cluster
  --family-assign           (optional) per-sample family-id TSV
                            cols: sample, family_id
  --froh-summary            (optional) per-sample F_ROH summary for quartile
                            stratification
                            cols: sample, F_ROH
  --constraint-proxy        (optional) per-gene constraint-positive flag TSV
                            cols: gene_id, constraint_pos   (0/1)
                            Without this, n_genes_constraint_pos = null and
                            per_group_cumulative falls back to n_genes_total.
  --n-top-peaks             how many recurrent ROH peaks to surface for Plot B
                            (default 24)
  --out                     output JSON path (default stdout)

Emits roh_gene_overlap_v1 with:
  params.{rank_criterion, constraint_proxy, gene_model_release}
  blocks[]        one entry per kept ROH tract (carries gene_ids + biotype_counts)
  peaks[]         top-N recurrence-ranked ROH peaks with biotype × cohort counts
  per_group_cumulative['K=8'][]   cumulative curves per cluster
  per_family_cumulative[]         per family-id
  per_quartile_cumulative[]       per F_ROH quartile (Q1..Q4)
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone


def _now_iso():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--roh-bed")
    p.add_argument("--genes-bed")
    p.add_argument("--cohort-recurrence")
    p.add_argument("--k8-assign")
    p.add_argument("--family-assign")
    p.add_argument("--froh-summary")
    p.add_argument("--constraint-proxy")
    p.add_argument("--n-top-peaks", type=int, default=24)
    p.add_argument("--gene-model-release", default=None)
    p.add_argument("--out")
    return p.parse_args()


def _write(out_path, payload):
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(payload, fh, indent=2)
        print(f"wrote {out_path}", file=sys.stderr)
    else:
        json.dump(payload, sys.stdout, indent=2)


def _emit_stub(out_path, gene_model_release=None):
    payload = {
        "_doc": ("roh_gene_overlap v0-stub — upstream pipeline outputs not "
                 "available at adapter run time. Run "
                 "STEP_A07_roh_derived_metrics.sh + STEP_A08_genes_in_roh.sh, "
                 "then re-run this adapter."),
        "version": "v0-stub",
        "generated_at": _now_iso(),
        "params": {
            "rank_criterion": None,
            "constraint_proxy": None,
            "gene_model_release": gene_model_release,
        },
        "blocks": [],
        "peaks": [],
        "per_group_cumulative": {"K=8": []},
        "per_family_cumulative": [],
        "per_quartile_cumulative": [],
    }
    _write(out_path, payload)


def main():
    import numpy as np
    import pandas as pd

    args = _parse_args()
    needed = [args.roh_bed, args.genes_bed, args.cohort_recurrence]
    if not all(p and os.path.exists(p) for p in needed):
        print("adapt_roh_gene_overlap: core inputs incomplete — emitting v0-stub",
              file=sys.stderr)
        _emit_stub(args.out, args.gene_model_release)
        return

    roh = pd.read_csv(args.roh_bed, sep="\t", header=None,
                      names=["chrom", "start", "end", "sample", "length"])
    genes = pd.read_csv(args.genes_bed, sep="\t", header=None,
                        names=["chrom", "start", "end", "gene_id",
                               "gene_name", "strand", "biotype"][:7],
                        usecols=lambda c: True).iloc[:, :7]
    if genes.shape[1] == 6:
        genes["biotype"] = "protein_coding"
    genes.columns = ["chrom", "start", "end", "gene_id",
                     "gene_name", "strand", "biotype"]
    rec = pd.read_csv(args.cohort_recurrence, sep="\t")

    # Optional joins
    k8 = _read_assign(args.k8_assign, "k8_cluster")
    fam = _read_assign(args.family_assign, "family_id")
    constraint = None
    if args.constraint_proxy and os.path.exists(args.constraint_proxy):
        constraint = pd.read_csv(args.constraint_proxy, sep="\t")
        constraint["constraint_pos"] = constraint["constraint_pos"].astype(int)

    # ---- blocks: one per kept ROH tract -------------------------------------
    roh = roh.copy()
    roh["block_id"] = (roh["sample"].astype(str) + ":" + roh["chrom"].astype(str)
                       + ":" + roh["start"].astype(str)
                       + "-" + roh["end"].astype(str))
    roh["len_mb"] = (roh["end"] - roh["start"]) / 1e6

    # Per-block gene intersections via a single-pass sort+merge per chrom.
    block_genes = _intersect(roh, genes)
    biotype_by_gene = dict(zip(genes["gene_id"], genes["biotype"]))
    cons_set = set(constraint.loc[constraint["constraint_pos"] == 1, "gene_id"]
                   .tolist()) if constraint is not None else None

    blocks_out = []
    for _, r in roh.iterrows():
        gids = block_genes.get(r["block_id"], [])
        biotype_counts = {}
        for g in gids:
            b = biotype_by_gene.get(g, "protein_coding")
            biotype_counts[b] = biotype_counts.get(b, 0) + 1
        n_cpos = (sum(1 for g in gids if g in cons_set)
                  if cons_set is not None else None)
        blocks_out.append({
            "block_id": r["block_id"],
            "chrom": str(r["chrom"]),
            "start": int(r["start"]),
            "end": int(r["end"]),
            "len_mb": float(r["len_mb"]),
            "rank_value": float(r["length"]),
            "sample_id": str(r["sample"]),
            "k_group": k8.get(str(r["sample"])),
            "family_id": fam.get(str(r["sample"])),
            "gene_ids": gids,
            "n_genes_total": len(gids),
            "n_genes_constraint_pos": n_cpos,
            "biotype_counts": biotype_counts,
        })

    # ---- peaks: top-N most recurrent (n_samples_with_overlap) ---------------
    # Group recurring genes by chrom + nearby start so we get peak labels
    # rather than per-gene rows; collapse with simple distance-based merge.
    peak_rows = []
    rec_sorted = rec.sort_values("n_samples_with_overlap", ascending=False)
    seen_regions = []
    PEAK_PAD_BP = 200_000
    for _, r in rec_sorted.iterrows():
        chrom = str(r["gene_chr"])
        start = int(r["gene_start"]); end = int(r["gene_end"])
        merged = False
        for region in seen_regions:
            if region["chrom"] != chrom:
                continue
            if start <= region["end"] + PEAK_PAD_BP and end >= region["start"] - PEAK_PAD_BP:
                region["start"] = min(region["start"], start)
                region["end"] = max(region["end"], end)
                region["genes"].append(r["gene_id"])
                region["cohort_count"] = max(region["cohort_count"],
                                             int(r["n_samples_with_overlap"]))
                merged = True
                break
        if not merged:
            seen_regions.append({
                "chrom": chrom, "start": start, "end": end,
                "genes": [r["gene_id"]],
                "cohort_count": int(r["n_samples_with_overlap"]),
            })
        if len(seen_regions) >= args.n_top_peaks * 3:
            break
    seen_regions.sort(key=lambda x: -x["cohort_count"])
    for i, region in enumerate(seen_regions[:args.n_top_peaks]):
        mid_mb = ((region["start"] + region["end"]) / 2) / 1e6
        label = f"ROH {region['chrom']} {mid_mb:.1f} Mb"
        biotype_counts = {}
        for g in region["genes"]:
            b = biotype_by_gene.get(g, "protein_coding")
            biotype_counts[b] = biotype_counts.get(b, 0) + 1
        peak_rows.append({
            "peak_id": f"peak_{i:03d}",
            "label": label,
            "chrom": region["chrom"],
            "start": int(region["start"]),
            "end": int(region["end"]),
            "cohort_sample_count": region["cohort_count"],
            "biotype_counts": biotype_counts,
        })

    # ---- per-group cumulative (K=8) -----------------------------------------
    blocks_df = pd.DataFrame(blocks_out)
    per_k = _per_group_cumulative(blocks_df, "k_group")
    per_family = _per_group_cumulative(blocks_df, "family_id")

    per_quartile = []
    if args.froh_summary and os.path.exists(args.froh_summary):
        froh = pd.read_csv(args.froh_summary, sep="\t")
        sample_col = next((c for c in ("sample", "sample_id") if c in froh.columns), None)
        if sample_col and "F_ROH" in froh.columns:
            q = pd.qcut(froh["F_ROH"], 4, labels=["Q1", "Q2", "Q3", "Q4"],
                        duplicates="drop")
            quartile_by_sample = dict(zip(froh[sample_col].astype(str), q.astype(str)))
            blocks_df["_qrt"] = blocks_df["sample_id"].map(quartile_by_sample)
            per_quartile = _per_group_cumulative(blocks_df.rename(columns={"_qrt": "k_group"}),
                                                 "k_group")

    payload = {
        "_doc": ("roh_gene_overlap derived from "
                 "catfish-diversity-analysis/05_roh_metrics_and_genes via "
                 "adapt_roh_gene_overlap.py."),
        "version": "v1-derived",
        "generated_at": _now_iso(),
        "params": {
            "rank_criterion": "length_bp",
            "constraint_proxy": ("external"
                                 if cons_set is not None else None),
            "gene_model_release": args.gene_model_release,
        },
        "blocks": blocks_out,
        "peaks": peak_rows,
        "per_group_cumulative": {"K=8": per_k},
        "per_family_cumulative": per_family,
        "per_quartile_cumulative": per_quartile,
    }
    _write(args.out, payload)


def _read_assign(path, value_col):
    if not path or not os.path.exists(path):
        return {}
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    sample_col = next((c for c in ("sample", "sample_id") if c in df.columns), None)
    if not sample_col or value_col not in df.columns:
        return {}
    return dict(zip(df[sample_col].astype(str), df[value_col].astype(str)))


def _intersect(roh_df, genes_df):
    """Return dict: block_id -> list of gene_ids overlapping the tract."""
    out = {}
    g_by_chrom = {c: g.sort_values("start").reset_index(drop=True)
                  for c, g in genes_df.groupby("chrom")}
    for _, r in roh_df.iterrows():
        g = g_by_chrom.get(str(r["chrom"]))
        if g is None or g.empty:
            out[r["block_id"]] = []
            continue
        mask = (g["end"] > r["start"]) & (g["start"] < r["end"])
        out[r["block_id"]] = g.loc[mask, "gene_id"].astype(str).tolist()
    return out


def _per_group_cumulative(blocks_df, group_col):
    """Sort each group's blocks by rank_value desc, then emit x=[0..N-1],
    y=cumulative n_genes_constraint_pos (or n_genes_total fallback)."""
    out = []
    if blocks_df.empty or group_col not in blocks_df.columns:
        return out
    use_cpos = blocks_df["n_genes_constraint_pos"].notna().any()
    for grp, g in blocks_df.dropna(subset=[group_col]).groupby(group_col):
        g_sorted = g.sort_values("rank_value", ascending=False).reset_index(drop=True)
        if use_cpos:
            counts = g_sorted["n_genes_constraint_pos"].fillna(0).astype(int).cumsum()
        else:
            counts = g_sorted["n_genes_total"].fillna(0).astype(int).cumsum()
        out.append({
            "group": str(grp),
            "x": list(range(len(counts))),
            "y": [int(v) for v in counts.tolist()],
            "n_samples": int(g_sorted["sample_id"].nunique()),
        })
    return out


if __name__ == "__main__":
    main()
