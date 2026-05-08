# 02_roh — Runs of homozygosity and F_ROH

ngsF-HMM ROH calls (10 replicates, best-by-likelihood) on the BEAGLE
genotype-likelihood input from the population-analysis sibling repo
(MODULE_2A). ROH BEDs are then parsed and F_ROH computed both
genome-wide and per-chromosome, alongside H inside vs. outside ROH.

## Pipeline

```
BEAGLE GLs + .pos + samples.ind   (from population-analysis)
  + per-sample θπ tracks           (from 01_het_theta_pi, for H in/out ROH)
  │
  ├─ STEP_A03_run_ngsF_HMM.sh        10 reps × seeds 42–51, keep best
  └─ STEP_A04_parse_roh_and_het.sh   .ibd → BED → F_ROH + H in/out ROH
       │
       ├─ utils/convert_ibd.pl         per-sample .ibd → BED tracts
       ├─ utils/summarize_ibd_roh.py   ROH/F_ROH summary with length bins
       └─ utils/parse_roh_and_het.py   master summary table
```

## Outputs

This module emits the **third of the three independently-consumable
result buckets**:

- **ROH bucket** — `${OUT_ROH}/`
  - `per_sample/<SAMPLE>.roh.bed` — ROH intervals per sample
  - `F_ROH.tsv` — F_ROH genome-wide and per-LG
  - `het_in_out_roh.tsv` — H contrast inside vs. outside ROH
  - ROH length bins: short (<1 Mb, ancestral), medium (1–5 Mb, historical),
    long (>5 Mb, recent inbreeding / consanguinity)

Consumed by inbreeding/conservation summaries, founder-pack analysis,
manuscript Results section on autozygosity, and the F_ROH|H framework.

This bucket can be produced and consumed independently of the
heterozygosity and θπ buckets — but the per-sample H from
`01_het_theta_pi` is needed if you want H in/out ROH, and the ngsF-HMM
input requires the BEAGLE GLs from the population-analysis sibling.

## Plots (called from `launchers/STEP_B01_run_all_plots.sh`)

- `plot_roh_core.R` — genome-wide ROH/F_ROH distributions
- `plot_roh_by_chromosome.R` — per-chromosome ROH heatmaps
- `plot_roh_metadata_overlays.R` — ancestry/relatedness overlays
- `plot_froh_heatmap.R` — F_ROH heatmap (sample × chromosome)

## F_ROH|H framework

The diversity-contextualised F_ROH (named methodological contribution of
the manuscript) is computed downstream from the outputs of this module —
it conditions F_ROH on local diversity to distinguish recent inbreeding
from ancestrally low-diversity regions. See `docs/methods/MODULE_3_methods.md`.
