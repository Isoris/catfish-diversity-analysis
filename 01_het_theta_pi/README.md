# 01_het_theta_pi — Heterozygosity + θπ

Per-sample SAF → folded SFS → genome-wide heterozygosity, plus per-sample
windowed θπ at multiple physical scales. Heterozygosity and θπ share the
SAF→SFS pipeline upstream and split downstream into the two result
buckets.

## Pipeline

```
BAMs (from variant-analysis) + callable mask
  │
  ├─ STEP_A01_prep_inputs.sh        validate inputs, build QC-pass BAM list
  └─ STEP_A02_run_heterozygosity.sh per-sample SAF → SFS → genome-wide H
                                     → windowed θπ (main + multiscale)
       └─ SLURM_A02_*.sh             array worker, one sample per task
```

For the 226-sample cohort, run from the repo root with
`bash LAUNCH.sh --step 2 --slurm` to submit the SLURM array.

## Outputs

Two of the three result buckets:

- **Heterozygosity** → `${OUT_HETEROZYGOSITY}/`
  Per-sample, genome-wide. SAF, SFS, and one H per individual.

- **θπ** → `${OUT_THETA_PI}/` (= `${OUT_HETEROZYGOSITY}/03_theta`)
  Cohort-level windowed θπ. Main scale 500 kb non-overlapping;
  multiscale set 5 kb/1 kb, 10 kb/2 kb, 50 kb/10 kb (configurable via
  `THETA_SCALES` in `00_config.sh`).

## θπ scaling caveat

ANGSD's pestPG `tP` column is a per-window sum, not a per-site density.
The diversity-comparable quantity is `tP / nSites`. This module currently
emits raw pestPG; the per-site adapter is a known TODO. See
`docs/methods/theta_pi_scaling.md`.

## Helpers used (from `utils/` at repo root)

`plot_heterozygosity_core.R`, `plot_theta_ideogram.R`,
`plot_theta_pi_multiscale.R`, `plot_scatter_stats.R`, `run_stats.R`,
`build_theta_supp_tables.py`, `run_theta_all_scales.sh`.
