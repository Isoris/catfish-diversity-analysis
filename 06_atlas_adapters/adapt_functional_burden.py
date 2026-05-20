#!/usr/bin/env python3
"""
adapt_functional_burden.py

Convert functional-burden / selection-efficacy pipeline outputs into
the functional_burden_v1.schema.json payload consumed by the
diversity-atlas `burden` page (five co-located layers: VESM burden,
πN/πS, π0/π4, LOF count, ROH-overlap fraction).

Upstream pipeline: NOT YET SHIPPED. Expected stack per spec §6.5:
  bcftools csq + snpEff + VESM + custom splice module + GERP track.

Until the pipeline ships, this adapter emits a schema-valid v0-stub
so the atlas page renders the 'data pending' fallback rather than
failing.

Expected inputs (when ready):

  --per-sample          per_sample_burden.tsv
                        cols: sample_id, k8, family_id, F_ROH,
                              pi, piN, piS, piN_piS,
                              piN_piS_ci_lo, piN_piS_ci_hi,
                              pi0, pi4, pi0_pi4,
                              pi0_pi4_ci_lo, pi0_pi4_ci_hi,
                              vesm_burden, vesm_n_sites,
                              lof_count, lof_count_strict, lof_count_loose,
                              splice_disruption_count,
                              lof_in_roh, vesm_in_roh_frac
  --per-group           per_group_burden.tsv  (group, n_samples, pi, ...)
                        Optional --per-group-key NAME=path (repeatable)
  --variant-inventory   variant_inventory.tsv  (type, impact, n)
  --snpeff-totals       snpeff_totals.tsv  (impact, n)
  --gerp-inventory      gerp_inventory.tsv  (bin, n, n_lof, n_missense)
  --transcripts-json    transcripts payload (passthrough)
  --msa-links-json      msa_links payload (passthrough)
  --splice-events       splice_events.tsv  (variant_id, event, ...)
  --pairwise-ks-json    pairwise_ks payload (passthrough, per spec §12)
  --pin-pis-method      e.g. 'NG86', 'fraction'
  --pi0-pi4-method      e.g. 'codon-degeneracy'
  --block-size-bp       int (default 5_000_000)
  --out                 output JSON path (default stdout)

When inputs are absent, emits the v0-stub shape.
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
    p.add_argument("--per-sample")
    p.add_argument("--per-group", help="K=8 per_group TSV (default key)")
    p.add_argument("--per-group-key", action="append", default=[],
                   help="NAME=path/to/per_group.tsv (repeatable; for "
                        "'family', 'F_ROH_quartile', etc.)")
    p.add_argument("--variant-inventory")
    p.add_argument("--snpeff-totals")
    p.add_argument("--gerp-inventory")
    p.add_argument("--transcripts-json")
    p.add_argument("--msa-links-json")
    p.add_argument("--splice-events")
    p.add_argument("--pairwise-ks-json")
    p.add_argument("--pin-pis-method", default=None)
    p.add_argument("--pi0-pi4-method", default=None)
    p.add_argument("--block-size-bp", type=int, default=5_000_000)
    p.add_argument("--out")
    return p.parse_args()


def _write(out_path, payload):
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(payload, fh, indent=2)
        print(f"wrote {out_path}", file=sys.stderr)
    else:
        json.dump(payload, sys.stdout, indent=2)


def _emit_stub(out_path, args):
    payload = {
        "_doc": ("functional_burden v0-stub — upstream bcftools csq + "
                 "snpEff + VESM + splice pipeline not shipped at adapter "
                 "run time. See 06_atlas_adapters/adapt_functional_burden.py "
                 "header for the expected input contract."),
        "version": "v0-stub",
        "generated_at": _now_iso(),
        "pipeline_steps": [],
        "estimator": {
            "piN_piS_method": args.pin_pis_method,
            "pi0_pi4_method": args.pi0_pi4_method,
            "ci_method": "block-jackknife across chromosomes",
            "block_size_bp": int(args.block_size_bp),
            "annotator_stack": ["bcftools csq", "snpEff", "VESM",
                                "custom splice module (TBD)"],
        },
        "cohort_summary": {
            "n_samples": 0, "pi": None,
            "piN_piS": None, "piN_piS_ci": [None, None],
            "pi0_pi4": None, "pi0_pi4_ci": [None, None],
            "vesm_burden": None, "lof_count": None, "splice_count": None,
        },
        "variant_inventory": [],
        "snpeff_totals": [],
        "gerp_inventory": None,
        "per_sample": [],
        "per_group": {"K=8": [], "pairwise_ks": {}},
        "alternative_stratifications": {
            "family": [], "F_ROH_quartile": [], "karyotype": {},
        },
        "top_burden_genes_by_group": {},
        "transcripts": {},
        "splice_events": [],
        "msa_links": {},
    }
    _write(out_path, payload)


def main():
    args = _parse_args()
    core_present = args.per_sample and os.path.exists(args.per_sample)
    if not core_present:
        print("adapt_functional_burden: --per-sample missing — emitting v0-stub",
              file=sys.stderr)
        _emit_stub(args.out, args)
        return

    import numpy as np
    import pandas as pd

    per_sample = _load_per_sample(args.per_sample)
    cohort = _summarise_cohort(per_sample)

    per_group_k8 = _load_per_group(args.per_group) if args.per_group else []
    alt_per_group = {}
    for spec in args.per_group_key:
        if "=" not in spec:
            print(f"adapt_functional_burden: bad --per-group-key '{spec}'",
                  file=sys.stderr)
            continue
        name, path = spec.split("=", 1)
        if os.path.exists(path):
            alt_per_group[name] = _load_per_group(path)

    variant_inventory = _load_table(args.variant_inventory,
                                    ["type", "impact", "n"])
    snpeff_totals = _load_table(args.snpeff_totals, ["impact", "n"])
    gerp_inventory = (_load_table(args.gerp_inventory,
                                  ["bin", "n", "n_lof", "n_missense"])
                      if args.gerp_inventory else None)
    splice_events = _load_table(args.splice_events,
                                ["variant_id", "event"], allow_extras=True)
    pairwise_ks = _load_json(args.pairwise_ks_json, default={})
    transcripts = _load_json(args.transcripts_json, default={})
    msa_links = _load_json(args.msa_links_json, default={})

    payload = {
        "_doc": ("functional_burden derived via adapt_functional_burden.py "
                 "from upstream bcftools csq + snpEff + VESM + splice outputs."),
        "version": "v1-derived",
        "generated_at": _now_iso(),
        "pipeline_steps": ["bcftools csq", "snpEff", "VESM",
                           "custom splice module"],
        "estimator": {
            "piN_piS_method": args.pin_pis_method,
            "pi0_pi4_method": args.pi0_pi4_method,
            "ci_method": "block-jackknife across chromosomes",
            "block_size_bp": int(args.block_size_bp),
            "annotator_stack": ["bcftools csq", "snpEff", "VESM",
                                "custom splice module"],
        },
        "cohort_summary": cohort,
        "variant_inventory": variant_inventory,
        "snpeff_totals": snpeff_totals,
        "gerp_inventory": gerp_inventory,
        "per_sample": per_sample,
        "per_group": {
            "K=8": per_group_k8,
            "pairwise_ks": pairwise_ks,
        },
        "alternative_stratifications": {
            "family": alt_per_group.get("family", []),
            "F_ROH_quartile": alt_per_group.get("F_ROH_quartile", []),
            "karyotype": {k: v for k, v in alt_per_group.items()
                          if k.startswith("karyotype:")},
        },
        "top_burden_genes_by_group": {},
        "transcripts": transcripts,
        "splice_events": splice_events,
        "msa_links": msa_links,
    }
    _write(args.out, payload)


def _load_per_sample(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    out = []
    for _, r in df.iterrows():
        out.append({
            "sample_id": str(r["sample_id"]),
            "k8": _s(r.get("k8")),
            "family_id": _s(r.get("family_id")),
            "F_ROH": _f(r.get("F_ROH")),
            "pi": _f(r.get("pi")),
            "piN": _f(r.get("piN")),
            "piS": _f(r.get("piS")),
            "piN_piS": _f(r.get("piN_piS")),
            "piN_piS_ci": [_f(r.get("piN_piS_ci_lo")), _f(r.get("piN_piS_ci_hi"))],
            "pi0": _f(r.get("pi0")),
            "pi4": _f(r.get("pi4")),
            "pi0_pi4": _f(r.get("pi0_pi4")),
            "pi0_pi4_ci": [_f(r.get("pi0_pi4_ci_lo")), _f(r.get("pi0_pi4_ci_hi"))],
            "vesm_burden": _f(r.get("vesm_burden")),
            "vesm_n_sites": _i(r.get("vesm_n_sites")),
            "lof_count": _i(r.get("lof_count")),
            "lof_count_strict": _i(r.get("lof_count_strict")),
            "lof_count_loose": _i(r.get("lof_count_loose")),
            "splice_disruption_count": _i(r.get("splice_disruption_count")),
            "lof_in_roh": _i(r.get("lof_in_roh")),
            "vesm_in_roh_frac": _f(r.get("vesm_in_roh_frac")),
        })
    return out


def _load_per_group(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    out = []
    for _, r in df.iterrows():
        out.append({
            "group": str(r["group"]),
            "n_samples": _i(r.get("n_samples")),
            "pi": _f(r.get("pi")),
            "piN_piS": _f(r.get("piN_piS")),
            "piN_piS_ci": [_f(r.get("piN_piS_ci_lo")), _f(r.get("piN_piS_ci_hi"))],
            "pi0_pi4": _f(r.get("pi0_pi4")),
            "pi0_pi4_ci": [_f(r.get("pi0_pi4_ci_lo")), _f(r.get("pi0_pi4_ci_hi"))],
            "vesm_burden_mean": _f(r.get("vesm_burden_mean")),
            "lof_count_mean": _f(r.get("lof_count_mean")),
            "roh_overlap_frac": _f(r.get("roh_overlap_frac")),
        })
    return out


def _summarise_cohort(per_sample):
    import numpy as np
    if not per_sample:
        return {
            "n_samples": 0, "pi": None,
            "piN_piS": None, "piN_piS_ci": [None, None],
            "pi0_pi4": None, "pi0_pi4_ci": [None, None],
            "vesm_burden": None, "lof_count": None, "splice_count": None,
        }

    def med(key):
        vals = [r[key] for r in per_sample if r.get(key) is not None]
        return float(np.median(vals)) if vals else None

    def total(key):
        vals = [r[key] for r in per_sample if r.get(key) is not None]
        return int(sum(vals)) if vals else None

    return {
        "n_samples": len(per_sample),
        "pi": med("pi"),
        "piN_piS": med("piN_piS"),
        "piN_piS_ci": [None, None],
        "pi0_pi4": med("pi0_pi4"),
        "pi0_pi4_ci": [None, None],
        "vesm_burden": med("vesm_burden"),
        "lof_count": total("lof_count"),
        "splice_count": total("splice_disruption_count"),
    }


def _load_table(path, cols, allow_extras=False):
    if not path or not os.path.exists(path):
        return []
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    out = []
    for _, r in df.iterrows():
        row = {c: _maybe_convert(r.get(c)) for c in cols if c in df.columns}
        if allow_extras:
            for c in df.columns:
                if c not in row:
                    row[c] = _maybe_convert(r[c])
        out.append(row)
    return out


def _load_json(path, default=None):
    if not path or not os.path.exists(path):
        return default
    with open(path) as fh:
        return json.load(fh)


def _maybe_convert(v):
    if v is None: return None
    try:
        f = float(v)
        if f.is_integer():
            return int(f)
        return f if (f == f) else None  # NaN guard
    except Exception:
        s = str(v).strip()
        return s if s and s.lower() != "nan" else None


def _f(v):
    try:
        import numpy as np
        if v is None: return None
        f = float(v)
        return None if not np.isfinite(f) else f
    except Exception:
        return None


def _i(v):
    try:
        if v is None: return None
        return int(v)
    except Exception:
        return None


def _s(v):
    try:
        if v is None: return None
        s = str(v).strip()
        return s if s and s.lower() != "nan" else None
    except Exception:
        return None


if __name__ == "__main__":
    main()
