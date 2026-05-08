# 02_roh — Runs of homozygosity + F_ROH

ngsF-HMM ROH calling on BEAGLE genotype likelihoods (from the
population-analysis sibling repo, MODULE_2A), parsed to per-sample BED,
with F_ROH genome-wide and per-LG plus H inside vs. outside ROH.

## Pipeline

```
BEAGLE GLs + .pos + samples.ind   (from population-analysis)
  + per-sample θπ tracks           (from 01_het_theta_pi, for H in/out ROH)
  │
  ├─ STEP_A03_run_ngsF_HMM.sh        10 reps × seeds 42–51, keep best by LL
  └─ STEP_A04_parse_roh_and_het.sh   .ibd → BED → F_ROH + H in/out ROH
       │
       ├─ utils/convert_ibd.pl         per-sample .ibd → BED tracts
       ├─ utils/summarize_ibd_roh.py   ROH/F_ROH summary with length bins
       └─ utils/parse_roh_and_het.py   master summary table
```

## Outputs

The third result bucket:

- **ROH / F_ROH** → `${OUT_ROH}/`
  - `per_sample/<SAMPLE>.roh.bed` — ROH intervals per sample
  - `F_ROH.tsv` — F_ROH genome-wide and per-LG
  - `het_in_out_roh.tsv` — H contrast inside vs. outside ROH
  - Length bins: short (<1 Mb, ancestral), medium (1–5 Mb, historical),
    long (>5 Mb, recent inbreeding / consanguinity)

## F_ROH|H framework

The diversity-contextualised F_ROH (named methodological contribution of
the manuscript) is computed downstream from the outputs of this module —
it conditions F_ROH on local diversity to distinguish recent inbreeding
from ancestrally low-diversity regions. See `docs/methods/MODULE_3_methods.md`.

## Helpers used (from `utils/` at repo root)

`plot_roh_core.R`, `plot_roh_by_chromosome.R`,
`plot_roh_metadata_overlays.R`, `plot_froh_heatmap.R`, `convert_ibd.pl`,
`summarize_ibd_roh.py`, `parse_roh_and_het.py`.
