# 04_window_diversity_texture — per-sample H texture (H_chr, DDI, χ_min)

Extends the genome-wide / per-sample H from
[01_het_theta_pi](../01_het_theta_pi/) into two within-genome texture
descriptors. Operates on the per-sample `.thetas.idx` already produced
by `STEP_A02_run_heterozygosity.sh`; no SAF/SFS recomputation.

## Why

Genome-wide H is one number per individual and collapses two distinct
biological signals:

1. **How much diversity does this fish carry?** (level)
2. **How is that diversity distributed across the genome?** (texture)

Two individuals with identical genome-wide H may have very different
window-level H profiles — *flat* (uniform across the genome) vs. *spiky*
(peaks adjacent to valleys, often co-occurring with heterozygous large
inversions, recent introgression, or shared ROH on a partly outbred
background). This module emits the metrics that separate these cases.

## Metrics

### Per-chromosome H

```
H_chr(sample) = sum_over_callable_sites_in_chr( theta_pi_site )
                / N_callable_chr
```

Recovers per-site θπ by exponentiating the natural-log values that
`thetaStat print` emits. `N_callable_chr` is the per-chromosome count of
the callable mask (`CALLABLE_SITES`), cached on first run.

### Diversity dispersion index (DDI)

```
DDI(sample) = MAD( H_w ) / median( H_w )
```

Per-sample within-genome dispersion of windowed H, with MAD (not SD)
for robustness against ROH-induced hard zeros and right-skew. DDI ≈ 0
flags a *flat* genome; high DDI flags a *spiky* (mosaic) genome.

### χ_min

```
chi_min(sample) = min over windows w of
                  [ H_w(sample) / median_cohort_smoothed( H_w )_at_w ]
```

Cohort-relative diversity floor, scale-free (1.0 = cohort median for
that window). The cohort denominator is per-window with a rolling
median over ±⌊DDI_SMOOTH_WIN/2⌋ neighbours to dampen single-window
noise.

## Scripts

| Step | Script | Purpose |
|------|--------|---------|
| A05 | `STEP_A05_per_chromosome_heterozygosity.sh` | per-chromosome H aggregation |
| A05 helper | `aggregate_per_chrom_H.py` | long / wide / cohort summary |
| A06 | `STEP_A06_window_H_and_DDI.sh` | windowed H, cohort median, DDI, χ_min |
| A06 helper | `compute_window_metrics.py` | two-pass cohort metric computation |

A05 streams `thetaStat print` through `awk` into per-sample per-chrom
sums (sum of `exp(log_tP)` per chromosome + observed-site count), then
the helper Python script joins the cohort `n_callable_chr` cache to
produce H_chr. A06 reuses an existing matching pestPG if one is
present in `${DIR_HET}/03_theta/` (main or `multiscale/`) and otherwise
generates a non-overlapping pestPG at `DDI_WIN`/`DDI_STEP` into
`pestPG_extra/` from the existing `.thetas.idx`.

## Outputs

Under `${OUT_WINDIV}` = `${OUTBASE}/11_window_diversity_texture/`:

```
01_per_chrom_H/
  callable_sites_per_chrom.tsv             cohort callable count cache
  by_sample/<SAMPLE>.per_chrom_tP.tsv      per-sample tP sums
  per_sample_per_chromosome_H.tsv          long: sample, chrom, H_chr
  cohort_sample_x_chrom_H.tsv              wide matrix
  cohort_per_chromosome_summary.tsv        per-chrom median / IQR / outliers

02_window_metrics/win<W>_step<S>/
  pestPG_manifest.tsv                      sample → pestPG used
  per_sample_window_metrics.tsv            sample, DDI, chi_min, ...
  cohort_window_H_matrix.tsv.gz            long (sample, chrom, win_start, H_w)
  cohort_window_median.tsv                 per-window cohort median + smoothed
02_window_metrics/pestPG_extra/win<W>_step<S>/
  <SAMPLE>.win<W>.step<S>.pestPG           generated only if not present upstream
```

## Configuration

All flags read from `00_config.sh` (override via environment):

| Variable | Default | Meaning |
|----------|---------|---------|
| `DDI_WIN` | `50000` | Primary window size (bp), non-overlapping |
| `DDI_STEP` | `50000` | Step size (must equal `DDI_WIN` — independence) |
| `DDI_MIN_CALLABLE_FRAC` | `0.3` | Drop windows below this × WIN callable sites |
| `DDI_SMOOTH_WIN` | `11` | Cohort-median rolling smoother (set 1 to disable) |
| `DDI_MAD_CONSTANT` | `1.4826` | Reported MAD scale (DDI is invariant in this) |

ROH segments are real biological zeros and should NOT be excluded by
the `DDI_MIN_CALLABLE_FRAC` filter — the filter targets missing-data
windows, not zero-H windows. Validate by spot-checking dropped windows
against the ROH calls in [02_roh/](../02_roh/) if the fraction is
tightened.

## Sensitivity sweep

Default scale is 50 kb non-overlapping. Run additional sweeps for
supplementary analyses:

```bash
DDI_WIN=10000  DDI_STEP=10000  bash STEP_A06_window_H_and_DDI.sh
DDI_WIN=100000 DDI_STEP=100000 bash STEP_A06_window_H_and_DDI.sh
DDI_WIN=500000 DDI_STEP=500000 bash STEP_A06_window_H_and_DDI.sh
```

DDI should be approximately stable across scales for an individual; if
the sample rank order shifts substantially, that itself is diagnostic
(scale-dependent texture).

## Dependencies

- Upstream: `STEP_A02_run_heterozygosity.sh` (produces `.thetas.idx`
  and one main pestPG per sample).
- Tools: `thetaStat` (ANGSD), `awk`, `python3` with `numpy` + `pandas`.

## Notes

- The MAD constant 1.4826 makes MAD a Gaussian-consistent estimator of
  σ. For non-Gaussian H distributions the raw MAD (constant = 1) is
  more honest, but DDI = MAD / median is a ratio and the constant
  cancels — the choice only affects the reported `mad_H_w` column.
- The χ_min smoothing default (11 windows = ±5 neighbours at 50 kb =
  ±250 kb) was chosen to suppress single-window noise without smearing
  inversion-scale signals. Set `DDI_SMOOTH_WIN=1` for unsmoothed
  cohort-window medians.
- `thetaStat print` emits log-space per-site theta estimates; we
  `exp()` them in `awk` before summing. This is consistent with the
  per-site convention documented in
  [docs/methods/theta_pi_scaling.md](../docs/methods/theta_pi_scaling.md).
