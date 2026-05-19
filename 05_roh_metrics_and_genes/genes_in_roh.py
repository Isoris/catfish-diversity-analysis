#!/usr/bin/env python3
"""
genes_in_roh.py — Two subcommands:

  gff-to-bed     parse a GFF/GFF3 and emit a 6-column BED of gene features:
                   chrom  start  end  gene_id  gene_name  strand

  aggregate      consume `bedtools intersect -wa -wb` output of
                   (per-sample ROH tracts × gene BED)
                 and emit per-sample + cohort gene-level tables, including
                 a "private ROH gene" set (genes overlapping ROH in exactly
                 one cohort sample — analogous to the wild-private-ROH gene
                 table referenced in the manuscript).

The intersection input must have the column layout produced by
STEP_A08_genes_in_roh.sh:

  0  roh_chr
  1  roh_start
  2  roh_end
  3  roh_sample
  4  roh_length
  5  gene_chr
  6  gene_start
  7  gene_end
  8  gene_id
  9  gene_name
 10  strand
"""

import argparse
import gzip
import os
import re
import sys

import pandas as pd


# ════════════════════════════════════════════════════════════════════════════
# gff-to-bed
# ════════════════════════════════════════════════════════════════════════════

def gff_open(path):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "r")


# GFF3 attribute parser: key=value;key=value...
ATTR_GFF3 = re.compile(r"([^=;\s]+)=([^;]+)")
# Loose GTF-style fallback: key "value"; key "value";
ATTR_GTF = re.compile(r'([^ =]+)\s+"([^"]+)"')


def parse_attrs(attr_str: str) -> dict:
    out = dict(ATTR_GFF3.findall(attr_str))
    if not out:
        out = {k: v for k, v in ATTR_GTF.findall(attr_str)}
    return out


def cmd_gff_to_bed(args):
    types_keep = set(t.strip() for t in args.feature_types.split(",") if t.strip())
    name_keys = [k.strip() for k in args.name_keys.split(",") if k.strip()]
    n_in = n_out = 0

    rows = []
    with gff_open(args.gff) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            n_in += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _src, ftype, start, end, _score, strand, _phase, attrs = f[:9]
            if ftype not in types_keep:
                continue
            attrs_d = parse_attrs(attrs)
            gene_id = (attrs_d.get("ID") or attrs_d.get("gene_id") or
                       attrs_d.get("Parent") or "")
            gene_name = ""
            for k in name_keys:
                if k in attrs_d and attrs_d[k]:
                    gene_name = attrs_d[k]
                    break
            if not gene_name:
                gene_name = gene_id
            try:
                s = int(start) - 1  # GFF is 1-based inclusive → BED 0-based
                e = int(end)
            except ValueError:
                continue
            if e <= s:
                continue
            rows.append((chrom, s, e, gene_id, gene_name, strand or "."))
            n_out += 1

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(args.out_bed, "w") as out:
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    print(f"  GFF rows scanned: {n_in:,}  features kept: {n_out:,}",
          file=sys.stderr)


# ════════════════════════════════════════════════════════════════════════════
# aggregate
# ════════════════════════════════════════════════════════════════════════════

INTERSECT_COLS = [
    "roh_chr", "roh_start", "roh_end", "roh_sample", "roh_length",
    "gene_chr", "gene_start", "gene_end", "gene_id", "gene_name", "strand",
]


def overlap_bp(row) -> int:
    return max(0, min(row["roh_end"], row["gene_end"]) -
                  max(row["roh_start"], row["gene_start"]))


def cmd_aggregate(args):
    df = pd.read_csv(args.intersection, sep="\t", header=None,
                     names=INTERSECT_COLS)
    if df.empty:
        print("  WARNING: empty intersection — no genes overlap any ROH",
              file=sys.stderr)
    for c in ("roh_start", "roh_end", "roh_length",
              "gene_start", "gene_end"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["roh_start", "roh_end",
                           "gene_start", "gene_end"]).copy()
    for c in ("roh_start", "roh_end", "roh_length",
              "gene_start", "gene_end"):
        df[c] = df[c].astype(int)

    df["overlap_bp"] = df.apply(overlap_bp, axis=1)
    df = df[df["overlap_bp"] >= args.min_overlap_bp].copy()
    df["gene_len"] = df["gene_end"] - df["gene_start"]
    df["overlap_frac_of_gene"] = df["overlap_bp"] / df["gene_len"]
    df["fully_contained"] = (
        (df["roh_start"] <= df["gene_start"]) &
        (df["roh_end"] >= df["gene_end"])
    )

    # ── Per-sample TSVs ─────────────────────────────────────────────────────
    os.makedirs(args.per_sample_dir, exist_ok=True)
    per_sample_keep = [
        "gene_id", "gene_name", "gene_chr", "gene_start", "gene_end",
        "strand", "roh_start", "roh_end", "roh_length",
        "overlap_bp", "overlap_frac_of_gene", "fully_contained",
    ]
    for sample, g in df.groupby("roh_sample"):
        g_out = g[per_sample_keep].rename(columns={"gene_chr": "chrom"}) \
                                  .sort_values(["chrom", "gene_start"])
        out = os.path.join(args.per_sample_dir,
                           f"{sample}.genes_in_roh.tsv")
        g_out.to_csv(out, sep="\t", index=False)

    # ── Per-sample summary ──────────────────────────────────────────────────
    # De-dupe to one row per (sample, gene) before counting unique-gene stats,
    # so a gene split across two ROH tracts in the same sample is counted once.
    unique_sg = (df.sort_values("overlap_bp", ascending=False)
                   .drop_duplicates(["roh_sample", "gene_id"]))
    per_sample_unique = (unique_sg.groupby("roh_sample")
                         .agg(n_genes_overlap=("gene_id", "size"),
                              n_genes_fully_contained=(
                                  "fully_contained",
                                  lambda s: int(s.sum())))
                         .reset_index())
    per_sample_overlap = (df.groupby("roh_sample")
                          .agg(total_gene_overlap_bp=("overlap_bp", "sum"))
                          .reset_index())
    per_sample = (per_sample_unique.merge(per_sample_overlap,
                                          on="roh_sample", how="left")
                  .rename(columns={"roh_sample": "sample"})
                  .sort_values("sample"))
    per_sample.to_csv(args.out_per_sample_summary, sep="\t", index=False)
    print(f"wrote {args.out_per_sample_summary}", file=sys.stderr)

    # ── Cohort long table (one row per gene × sample) ───────────────────────
    cohort_long = (df
                   .groupby(["gene_id", "gene_name", "gene_chr",
                             "gene_start", "gene_end", "strand",
                             "roh_sample"])
                   .agg(overlap_bp=("overlap_bp", "sum"),
                        fully_contained=("fully_contained", "any"))
                   .reset_index()
                   .rename(columns={"roh_sample": "sample"}))
    cohort_long.to_csv(args.out_cohort_long, sep="\t", index=False)
    print(f"wrote {args.out_cohort_long}", file=sys.stderr)

    # ── Per-gene cohort recurrence ──────────────────────────────────────────
    gene_keys = ["gene_id", "gene_name", "gene_chr", "gene_start",
                 "gene_end", "strand"]
    n_cohort = df["roh_sample"].nunique() if not df.empty else 0
    rec = (cohort_long.groupby(gene_keys)
           .agg(n_samples_with_overlap=("sample", "nunique"),
                samples=("sample", lambda s: ",".join(sorted(set(s)))),
                n_fully_contained=("fully_contained",
                                   lambda s: int(s.sum())),
                total_overlap_bp=("overlap_bp", "sum"))
           .reset_index())
    if n_cohort > 0:
        rec["frac_cohort"] = rec["n_samples_with_overlap"] / n_cohort
    else:
        rec["frac_cohort"] = 0.0
    rec = rec.sort_values(["n_samples_with_overlap", "gene_chr",
                           "gene_start"], ascending=[False, True, True])
    rec.to_csv(args.out_cohort_recurrence, sep="\t", index=False)
    print(f"wrote {args.out_cohort_recurrence}", file=sys.stderr)

    # ── Private-ROH genes (overlap in exactly 1 sample) ─────────────────────
    private = rec[rec["n_samples_with_overlap"] == 1].copy()
    private = private.rename(columns={"samples": "sample"})
    private.to_csv(args.out_private, sep="\t", index=False)
    print(f"wrote {args.out_private}  ({len(private)} private-ROH genes)",
          file=sys.stderr)


# ════════════════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("gff-to-bed")
    g.add_argument("--gff", required=True)
    g.add_argument("--feature-types", default="gene",
                   help="Comma-separated GFF type column to keep")
    g.add_argument("--name-keys", default="Name,gene_name,gene",
                   help="Comma-separated attribute keys to try for gene_name")
    g.add_argument("--out-bed", required=True)
    g.set_defaults(func=cmd_gff_to_bed)

    a = sub.add_parser("aggregate")
    a.add_argument("--intersection", required=True,
                   help="bedtools intersect -wa -wb output (ROH × genes)")
    a.add_argument("--min-overlap-bp", type=int, default=1)
    a.add_argument("--per-sample-dir", required=True)
    a.add_argument("--out-per-sample-summary", required=True)
    a.add_argument("--out-cohort-long", required=True)
    a.add_argument("--out-cohort-recurrence", required=True)
    a.add_argument("--out-private", required=True)
    a.set_defaults(func=cmd_aggregate)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
