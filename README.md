# catfish-diversity-analysis

Per-sample genetic diversity analysis for the *Clarias gariepinus* hatchery
cohort (n = 226): heterozygosity, nucleotide diversity θπ, and runs of
homozygosity.

This is one of three sibling **catfish-{population,diversity,variant}-analysis**
repos that together produce the population-genetic primitives consumed by
[`inversion-atlas`](https://github.com/Isoris/inversion-atlas) and any future
papers on this cohort.

## What this repo produces

Three independently-consumable result buckets. With these three, ~50% of
standard genetic-diversity and conservation analyses downstream can be
run without any further pipeline code:

| Bucket | Granularity | Description |
|---|---|---|
| **Heterozygosity** | per-sample, genome-wide | One H per individual from a per-sample folded SFS (ANGSD `-doSaf 1 -fold 1` → realSFS). The per-sample primitive for inbreeding metrics, conservation summaries, and Table 1. |
| **θπ (pestPG)** | cohort-level, windowed, multi-scale | Per-window θπ at multiple physical scales (default main 500 kb; multiscale set 5 kb/1 kb, 10 kb/2 kb, 50 kb/10 kb). The local-diversity primitive for selection scans, inversion-region contrast, and the inversion-atlas page-12 enrichment. |
| **ROH / F_ROH** | per-sample, regional + summary | ngsF-HMM ROH calls (10 reps × seeds 42–51, best-by-likelihood), parsed to BED, with F_ROH genome-wide and per-LG, length bins (short/medium/long), and H inside vs. outside ROH. The autozygosity primitive for inbreeding analysis and the F_ROH&#124;H framework. |

The buckets share inputs (BAMs, callable mask) but are independent
downstream — a consumer can take any one without needing the others.

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
│   ├── 01_het_theta_pi/          per-sample SAF → SFS → genome-wide H
│   │                             + windowed θπ (main + multiscale)
│   └── 02_roh/                   ngsF-HMM ROH + F_ROH + H in/out ROH
├── launchers/
│   ├── LAUNCH_module3.sh         orchestrates 01 → 02 sequentially or via SLURM
│   ├── STEP_B01_run_all_plots.sh aggregates plots across both modules
│   └── write_report.py           auto-generated Methods/Results markdown
├── envs/
├── docs/
│   └── methods/
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

Pipelines from the legacy `MODULE_3_heterozygosity_roh/` tree have been
landed under `Modules/01_het_theta_pi/` (heterozygosity + θπ; they share
the SAF→SFS pipeline) and `Modules/02_roh/` (ngsF-HMM, F_ROH). Top-level
launcher and plot orchestrator live in `launchers/`. Scripts run as-is on
LANTA — no logic was rewritten in the move, only re-homed.

The pestPG `tP / nSites` per-site adapter promised by `03_theta_pi`'s
output contract is still a TODO inside `01_het_theta_pi` — current
emissions are raw pestPG. See `docs/methods/theta_pi_scaling.md`.

## Citation

Project umbrella DOI: TBD (Zenodo, will be issued at v1.0 tag).

## Cohort note

**This repo is for the 226-sample pure *Clarias gariepinus* hatchery cohort
only.** Do not use it for the F₁ hybrid (*C. gariepinus* × *C. macrocephalus*)
genome assembly cohort or any future *C. macrocephalus* wild cohort — those
are separate manuscripts and may need different parameter choices.
