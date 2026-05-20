#!/usr/bin/env python3
"""
adapt_divergence_network.py

Convert pairwise group-divergence outputs (FST / DXY / dA / Nei's D)
into the divergence_network_v1.schema.json payload consumed by the
diversity-atlas `divergence` page.

Upstream pipeline: NOT YET SHIPPED. The expected inputs (when ready)
are documented below. Until then, the adapter emits a schema-valid
v0-stub so the atlas page renders the 'data pending' fallback rather
than failing.

Expected inputs (when the FST pipeline ships):

  --groups-tsv          per-group node attributes
                        cols: id, n_samples, pi_within, pi_within_ci_lo,
                              pi_within_ci_hi, h_mean, f_roh_mean, theta_w,
                              colour
  --edges-tsv           pairwise group edges
                        cols: i, j, fst, fst_ci_lo, fst_ci_hi,
                              dxy, dxy_ci_lo, dxy_ci_hi,
                              da,  da_ci_lo,  da_ci_hi,
                              nei_d, n_i, n_j
  --grouping            label for params.grouping (default 'K=8')
  --fst-estimator       'Weir-Cockerham' | 'Hudson' | 'Reynolds'
  --n-bootstrap         integer
  --min-callable        min callable sites per group pair
  --alt-grouping        repeatable: NAME=path/to/groups.tsv:path/to/edges.tsv
  --out                 output JSON path (default stdout)

Each --alt-grouping bundle is parsed identically to the primary
groups+edges pair and attached under alternative_groupings[NAME].

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
    p.add_argument("--groups-tsv")
    p.add_argument("--edges-tsv")
    p.add_argument("--grouping", default="K=8")
    p.add_argument("--fst-estimator",
                   choices=["Weir-Cockerham", "Hudson", "Reynolds"],
                   default=None)
    p.add_argument("--n-bootstrap", type=int, default=None)
    p.add_argument("--min-callable", type=int, default=None)
    p.add_argument("--alt-grouping", action="append", default=[],
                   help="NAME=groups.tsv:edges.tsv (repeatable)")
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
        "_doc": ("divergence_network v0-stub — upstream pairwise FST/DXY "
                 "pipeline not shipped at adapter run time. See "
                 "06_atlas_adapters/adapt_divergence_network.py header for "
                 "the expected input contract."),
        "version": "v0-stub",
        "generated_at": _now_iso(),
        "params": {
            "grouping": args.grouping,
            "fst_estimator": args.fst_estimator,
            "n_bootstrap": args.n_bootstrap,
            "min_callable_sites_per_group_pair": args.min_callable,
        },
        "groups": [],
        "edges": [],
        "alternative_groupings": {},
    }
    _write(out_path, payload)


def main():
    args = _parse_args()
    have_core = (args.groups_tsv and os.path.exists(args.groups_tsv)
                 and args.edges_tsv and os.path.exists(args.edges_tsv))
    if not have_core:
        print("adapt_divergence_network: core inputs missing — emitting v0-stub",
              file=sys.stderr)
        _emit_stub(args.out, args)
        return

    import pandas as pd

    groups = _build_groups(args.groups_tsv)
    edges = _build_edges(args.edges_tsv)

    # Alternative groupings — each value is its own {groups, edges} bundle.
    alt = {}
    for spec in args.alt_grouping:
        if "=" not in spec or ":" not in spec.split("=", 1)[1]:
            print(f"adapt_divergence_network: bad --alt-grouping '{spec}' "
                  "(expected NAME=groups.tsv:edges.tsv), skipping",
                  file=sys.stderr)
            continue
        name, paths = spec.split("=", 1)
        gp, ep = paths.split(":", 1)
        if not (os.path.exists(gp) and os.path.exists(ep)):
            print(f"adapt_divergence_network: alt-grouping {name} files missing, "
                  "skipping", file=sys.stderr)
            continue
        alt[name] = {
            "groups": _build_groups(gp),
            "edges":  _build_edges(ep),
        }

    payload = {
        "_doc": ("divergence_network derived via adapt_divergence_network.py."),
        "version": "v1-derived",
        "generated_at": _now_iso(),
        "params": {
            "grouping": args.grouping,
            "fst_estimator": args.fst_estimator,
            "n_bootstrap": args.n_bootstrap,
            "min_callable_sites_per_group_pair": args.min_callable,
        },
        "groups": groups,
        "edges": edges,
        "alternative_groupings": alt,
    }
    _write(args.out, payload)


def _build_groups(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    out = []
    for _, r in df.iterrows():
        out.append({
            "id": str(r["id"]),
            "n_samples": int(r["n_samples"]),
            "pi_within": _f(r.get("pi_within")),
            "pi_within_ci95": [_f(r.get("pi_within_ci_lo")),
                               _f(r.get("pi_within_ci_hi"))],
            "h_mean": _f(r.get("h_mean")),
            "f_roh_mean": _f(r.get("f_roh_mean")),
            "theta_w": _f(r.get("theta_w")),
            "colour": _s(r.get("colour")),
        })
    return out


def _build_edges(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    out = []
    for _, r in df.iterrows():
        out.append({
            "i": str(r["i"]),
            "j": str(r["j"]),
            "fst": _f(r.get("fst")),
            "fst_ci95": [_f(r.get("fst_ci_lo")), _f(r.get("fst_ci_hi"))],
            "dxy": _f(r.get("dxy")),
            "dxy_ci95": [_f(r.get("dxy_ci_lo")), _f(r.get("dxy_ci_hi"))],
            "da": _f(r.get("da")),
            "da_ci95": [_f(r.get("da_ci_lo")), _f(r.get("da_ci_hi"))],
            "nei_d": _f(r.get("nei_d")),
            "n_i": _i(r.get("n_i")),
            "n_j": _i(r.get("n_j")),
        })
    return out


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
