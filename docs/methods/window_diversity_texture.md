# Window-level diversity texture — H_chr, DDI, χ_min

Methods text (paper-ready) for the per-sample within-genome H texture
metrics produced by [04_window_diversity_texture/](../../04_window_diversity_texture/).

## Methods

### Per-chromosome heterozygosity (H_chr)

For each sample, per-site nucleotide diversity θπ was recovered from
the per-sample `.thetas.idx` produced upstream (`STEP_A02`) by
`thetaStat print`, which emits natural-log per-site θπ at every
callable site under the cohort callable mask. Per-chromosome
heterozygosity was computed by exponentiating the log values, summing
per-site θπ across callable sites within each chromosome, and
normalising by the per-chromosome callable-site count:

```
H_chr(sample) = Σ_{i ∈ chr ∩ callable} θπ_i  /  N_callable_chr
```

The denominator was obtained by counting rows of the ANGSD callable
mask per chromosome (cached cohort-wide on first run). H_chr is
directly comparable to the genome-wide SFS-based H of `STEP_A02`
because both share the same callable-site denominator construction.

### Diversity dispersion index (DDI)

For each sample, windowed per-site heterozygosity was computed as
`H_w = tP / nSites` from ANGSD `thetaStat do_stat` over non-overlapping
50 kb physical windows (`-type 2`, anchored to chromosome start;
windows with fewer than 30 % of the nominal window size in callable
sites were excluded). Within-genome heterozygosity dispersion was then
summarised by

```
DDI(sample) = MAD( H_w )  /  median( H_w )
```

where MAD is the median absolute deviation. The ratio is robust to
ROH-induced hard zeros and to right-skew in the windowed H
distribution, is scale-free, and is comparable across individuals with
different genome-wide H. DDI ≈ 0 indicates a uniform window-level H
profile ("flat"); high DDI indicates a mosaic profile mixing high-H
peaks and low-H valleys ("spiky"). The MAD scale constant 1.4826
(Gaussian-consistent) was used for the reported MAD; the DDI ratio is
invariant in this constant.

### Cohort-relative diversity floor (χ_min)

To identify the most depleted region per individual relative to the
cohort, we computed

```
chi_min(sample) = min_w [ H_w(sample) / median_cohort_smoothed( H_w )_at_w ]
```

where `median_cohort_smoothed( H_w )_at_w` is the cohort median of
`H_w` at window *w*, smoothed by a centred rolling median over ±5
neighbouring windows (total span 11 windows = 550 kb at the 50 kb
scale) to suppress single-window noise. χ_min is unitless: 1.0 = at
cohort median for the most depleted window; values < 1 quantify how
far below cohort median the individual's diversity floor sits, with
the corresponding window position retained for follow-up.

### Implementation notes

- `tP` is a per-window sum across sites; the per-site comparable
  quantity is `tP / nSites`, in line with the repo-wide convention
  documented in [theta_pi_scaling.md](theta_pi_scaling.md).
- ROH segments are biological zeros, not missing data, and are
  retained in DDI / χ_min computation. The min-callable filter
  excludes windows where the data, not the diversity, is the zero.
- DDI / χ_min are computed at 50 kb non-overlapping windows for
  primary analyses. Sensitivity sweeps at 10 kb / 100 kb / 500 kb are
  available via the `DDI_WIN` / `DDI_STEP` environment overrides; rank
  stability of DDI across scales is itself a diagnostic of the
  texture's spatial structure.

## Key parameters

| Parameter | Value | Used in |
|-----------|-------|---------|
| Window size for DDI / χ_min | 50 kb non-overlapping | A06 |
| Min callable fraction per window | 0.3 × WIN (i.e. 15 kb @ 50 kb) | A06 |
| Cohort-median smoothing | 11 windows centred (±5) | A06 |
| MAD scale constant (reported) | 1.4826 | A06 |
| Per-site θπ source | `exp()` of `thetaStat print` log values | A05 |
| Per-chr H denominator | per-chrom row count of CALLABLE_SITES | A05 |

## Outputs referenced in Results

- `per_sample_per_chromosome_H.tsv` — long table for H_chr distributions.
- `cohort_sample_x_chrom_H.tsv` — wide sample × chromosome matrix.
- `per_sample_window_metrics.tsv` — one row per sample with median H_w,
  MAD H_w, DDI, χ_min, χ_min chromosome and start position, and the
  number of windows kept after the min-callable filter.
- `cohort_window_H_matrix.tsv.gz` — sample × window long table for
  downstream per-window inspection.
- `cohort_window_median.tsv` — per-window cohort median (raw and
  smoothed) used as the χ_min denominator.
