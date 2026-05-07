# catfish-diversity-analysis

Per-sample genetic diversity analysis for the *Clarias gariepinus* hatchery
cohort (n = 226): heterozygosity, nucleotide diversity θπ, and runs of
homozygosity.

This is one of three sibling **catfish-{population,diversity,variant}-analysis**
repos that together produce the population-genetic primitives consumed by
[`inversion-atlas`](https://github.com/Isoris/inversion-atlas) and any future
papers on this cohort.

## What this repo produces

Per-sample summaries of genetic diversity, callable-site-aware:

| Output | Description |
|---|---|
| Per-sample SAF (folded) | from ANGSD `-doSaf 1 -fold 1` |
| Per-sample θπ (pestPG, multi-scale) | window/step grid: 10 kb / 2 kb default |
| Per-sample heterozygosity (genome-wide + per-LG) | from per-sample SFS |
| Per-(sample, window) θπ matrix (TSV) | tidy long format, one row per (sample, window) |
| ROH calls per sample | bcftools/PLINK ROH on biSNP set |
| F<sub>ROH</sub> + diversity-contextualized F<sub>ROH</sub>|H | named framework from the manuscript |

The per-(sample, window) θπ matrix is the input that `inversion-atlas`
consumes for its page-12 θπ enrichment.

## Inputs (consumed from sibling repos)

| Input | From |
|---|---|
| BAMs + BAI | `catfish-variant-analysis` (or upstream MODULE_1 read prep) |
| Callable mask | `catfish-variant-analysis` |
| Reference FASTA | shared `00-samples/` on LANTA |
| biSNP list (for ROH only) | `catfish-variant-analysis` |

## Engines used

- ANGSD (system) — SAF, thetaStat, pestPG
- [`angsd_fixed_HWE`](https://github.com/Isoris/angsd_fixed_HWE) — patched
  ANGSD with fixed-F EM, used where standard HWE-based MAF estimation is
  inappropriate
- bcftools roh / PLINK (system) — ROH calling

## Layout

```text
catfish-diversity-analysis/
├── 00_config.sh                  root config
├── Modules/
│   ├── 01_saf_per_sample/        per-sample folded SAF
│   ├── 02_heterozygosity/        per-sample SFS → genome-wide H
│   ├── 03_theta_pi/              pestPG multi-scale + per-(sample,window) TSV
│   ├── 04_roh/                   ROH on biSNP Beagle GLs
│   └── 05_aggregated/            tidy matrices for downstream consumers
├── envs/
├── docs/
│   ├── module_contracts/
│   └── methods/
├── tests/
└── README.md
```

## Output contracts (what downstream consumers can rely on)

The atlas + inversion-analysis repos read these well-known paths:

```
${OUT_THETA_PI}/aggregated/theta_native.<CHROM>.<SCALE>.tsv.gz
    sample  chrom  window_idx  start_bp  end_bp  theta_pi  tP_sum  n_sites

${OUT_HETEROZYGOSITY}/per_sample/<SAMPLE>.het.tsv
    chrom  het_genome_wide  het_per_lg

${OUT_ROH}/per_sample/<SAMPLE>.roh.tsv
    chrom  start_bp  end_bp  length_bp  state
```

The θπ TSV uses **per-site** values (`tP / nSites`) by default — this is the
diversity-comparable estimate, not the raw pestPG `tP` window sum. See
`docs/methods/theta_pi_scaling.md` for the rationale (ANGSD GitHub issue
#329; Korunes & Samuk 2021 / pixy). The raw window sum is preserved as
`tP_sum` for diagnostic display.

## Status

Scaffold. Pipelines exist and have been run on LANTA but live outside any
git repo today (under `${BASE}/het_roh/`). They will be migrated into
`Modules/` over time, one module at a time, as part of the manuscript-prep
cleanup.

## Citation

Project umbrella DOI: TBD (Zenodo, will be issued at v1.0 tag).

## Cohort note

**This repo is for the 226-sample pure *Clarias gariepinus* hatchery cohort
only.** Do not use it for the F₁ hybrid (*C. gariepinus* × *C. macrocephalus*)
genome assembly cohort or any future *C. macrocephalus* wild cohort — those
are separate manuscripts and may need different parameter choices.
