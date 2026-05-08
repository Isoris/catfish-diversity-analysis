# 01_het_theta_pi — Per-sample heterozygosity + θπ tracks

Per-sample SAF → folded SFS → genome-wide heterozygosity, plus per-sample
windowed θπ at multiple scales. These two quantities share the same
ANGSD/realSFS pipeline up to the SFS step, so they live together.

## Pipeline

```
BAMs (from variant-analysis)
  + callable mask
  + sample list (built by STEP_A01)
  │
  ├─ STEP_A01_prep_inputs.sh        validate inputs, build QC-pass BAM list
  └─ STEP_A02_run_heterozygosity.sh  per-sample SAF → SFS → genome-wide H
                                     → windowed θπ (main + multiscale)
       └─ SLURM_A02_*.sh              array worker, one sample per task
```

Sequential or SLURM-array. For the 226-sample cohort, use `--slurm` from
the top-level launcher (`launchers/LAUNCH_module3.sh --step 2 --slurm`).

## Outputs

- `${OUT_HETEROZYGOSITY}/` — per-sample SAF, SFS, genome-wide H summary table
- `${OUT_THETA_PI}/` — windowed θπ at the main scale (default 500 kb) plus
  the multiscale set (`THETA_SCALES` in `00_config.sh`)

## θπ scaling caveat

ANGSD's pestPG `tP` column is a per-window sum, not a per-site density.
The diversity-comparable quantity is `tP / nSites`. This module currently
emits raw pestPG; the per-site adapter is a known TODO (see
`docs/methods/theta_pi_scaling.md`).

## Plots / stats (called from `launchers/STEP_B01_run_all_plots.sh`)

- `plot_heterozygosity_core.R` — genome-wide H distribution
- `plot_theta_ideogram.R` — θπ ideogram tracks
- `plot_theta_pi_multiscale.R` — multi-scale θπ comparison
- `plot_scatter_stats.R` — H vs depth / ROH burden scatters
- `run_stats.R` — sample- and chromosome-level summary stats
- `build_theta_supp_tables.py` — θπ supplementary tables
