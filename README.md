# catfish-diversity-analysis

Genome-wide diversity, runs of homozygosity, and per-inversion band
diversity for the 226-sample pure *C. gariepinus* hatchery cohort.

Three primary result buckets are produced by the het/ROH chain (modules
[01](01_het_theta_pi/) + [02](02_roh/) + [STEP_B01](STEP_B01_run_all_plots.sh)):

- **Heterozygosity** — per-sample, genome-wide, from SAF → folded SFS.
- **θπ** — cohort-level windowed nucleotide diversity, multi-scale.
- **ROH / F_ROH** — ngsF-HMM ROH calls, FROH, H inside vs. outside ROH.

Module [03](03_pi_inversion_karyotype/) adds a separately-run downstream
analysis: per-inversion **band-level π** (π11 / π22 / π12) for inversion
candidates with karyotype calls from `catfish-inversion-analysis`.

Module [04](04_window_diversity_texture/) extends the per-sample H
into **within-genome texture**: per-chromosome H (`H_chr`), the
diversity dispersion index (`DDI = MAD(H_w)/median(H_w)`), and the
cohort-relative diversity floor (`χ_min`). Operates on the
`.thetas.idx` already produced by module 01 — no SAF/SFS
recomputation.

Module [05](05_roh_metrics_and_genes/) adds **ROH-derived metrics**
(`N_ROH`, `S_ROH`, `F_HOM`) and a **genes-in-ROH** annotation against
a GFF/GFF3, with per-sample views, a cohort recurrence table, and a
"private-ROH gene" set (analogue of the wild-private-ROH gene table).
Operates on the per-sample ROH BED + summary produced by module 02 —
no ngsF-HMM re-run.

## Cohort scope

This repo operates exclusively on the 226-sample pure *C. gariepinus*
hatchery cohort. It does **not** apply to the F₁ hybrid assembly cohort
or any *C. macrocephalus* wild cohort.

## Layout

```
catfish-diversity-analysis/
├── 00_config.sh                       single source of truth for paths / SLURM
├── LAUNCH.sh                          orchestrates het/ROH chain (steps 1–5)
├── STEP_B01_run_all_plots.sh          plots + stats + report
├── write_report.py                    methods/results report generator
│
├── 01_het_theta_pi/                   SAF → SFS → H + windowed θπ
│   ├── STEP_A01_prep_inputs.sh
│   ├── STEP_A02_run_heterozygosity.sh
│   ├── SLURM_A02_heterozygosity_worker.sh
│   └── README.md
│
├── 02_roh/                            ngsF-HMM ROH + F_ROH + H in/out ROH
│   ├── STEP_A03_run_ngsF_HMM.sh
│   ├── STEP_A04_parse_roh_and_het.sh
│   └── README.md
│
├── 03_pi_inversion_karyotype/         per-inversion band-π (π11 / π22 / π12)
│   ├── STEP_PI_B_make_karyotype_group_theta_pi.sh   driver
│   ├── STEP_PI_C_make_karyotype_bamlists.R
│   ├── STEP_PI_D_compute_pi11_pi22_pi12.R
│   ├── STEP_PI_E_average_sample_pestPG_by_karyotype.R
│   ├── STEP_PI_F_plot_band_pi_panel.R
│   └── README.md
│
├── 04_window_diversity_texture/       per-sample H texture: H_chr, DDI, χ_min
│   ├── STEP_A05_per_chromosome_heterozygosity.sh
│   ├── STEP_A06_window_H_and_DDI.sh
│   ├── aggregate_per_chrom_H.py
│   ├── compute_window_metrics.py
│   └── README.md
│
├── 05_roh_metrics_and_genes/          N_ROH / S_ROH / F_HOM + genes-in-ROH
│   ├── STEP_A07_roh_derived_metrics.sh
│   ├── STEP_A08_genes_in_roh.sh
│   ├── compute_roh_derived_metrics.py
│   ├── genes_in_roh.py
│   └── README.md
│
├── docs/methods/                      methods write-ups
│   ├── MODULE_3_methods.md            het / ROH / FROH methods
│   ├── theta_pi_scaling.md            tP vs. tP/nSites caveat
│   ├── window_diversity_texture.md    H_chr / DDI / χ_min methods
│   └── roh_metrics_and_genes.md       N_ROH / S_ROH / F_HOM / genes-in-ROH
│
└── utils/                             shared helpers (.ibd → BED, parsers)
```

## Running the het / ROH chain

[LAUNCH.sh](LAUNCH.sh) orchestrates the five-step pipeline that produces
the three primary result buckets:

```bash
bash LAUNCH.sh                       # everything sequentially
bash LAUNCH.sh --step 2 --slurm      # step 2 as a 226-sample SLURM array
bash LAUNCH.sh --from 3              # resume after the SLURM array completes
```

Step map:

| Step | Module | Script |
|------|--------|--------|
| 1 | [01_het_theta_pi](01_het_theta_pi/) | `STEP_A01_prep_inputs.sh` |
| 2 | [01_het_theta_pi](01_het_theta_pi/) | `STEP_A02_run_heterozygosity.sh` (or SLURM array) |
| 3 | [02_roh](02_roh/) | `STEP_A03_run_ngsF_HMM.sh` |
| 4 | [02_roh](02_roh/) | `STEP_A04_parse_roh_and_het.sh` |
| 5 | (root) | `STEP_B01_run_all_plots.sh` |
| 6 | [04_window_diversity_texture](04_window_diversity_texture/) | `STEP_A05_per_chromosome_heterozygosity.sh` |
| 7 | [04_window_diversity_texture](04_window_diversity_texture/) | `STEP_A06_window_H_and_DDI.sh` |
| 8 | [05_roh_metrics_and_genes](05_roh_metrics_and_genes/) | `STEP_A07_roh_derived_metrics.sh` |
| 9 | [05_roh_metrics_and_genes](05_roh_metrics_and_genes/) | `STEP_A08_genes_in_roh.sh` |

## Running the band-π module (independent)

Module 03 is **not** included in `LAUNCH.sh` — it depends on inversion
candidates and karyotype calls produced by the
`catfish-inversion-analysis` sibling repo, and is invoked separately:

```bash
# All inversions in INV_CANDIDATES
bash 03_pi_inversion_karyotype/STEP_PI_B_make_karyotype_group_theta_pi.sh

# A subset
bash 03_pi_inversion_karyotype/STEP_PI_B_make_karyotype_group_theta_pi.sh INV01 INV03
```

See [03_pi_inversion_karyotype/README.md](03_pi_inversion_karyotype/README.md)
for inputs, the shared 4-column sites file mechanism (which prevents the
silent π12 polarity flip), and the `variant_only` vs `all_callable`
site-basis caveat.

## Configuration

All paths, parameters, and SLURM defaults live in
[00_config.sh](00_config.sh). Source it from any script at the repo root:

```bash
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
```

Key environment variables every module reads: `BASE`, `REF`,
`CALLABLE_SITES`, `SAMPLE_MANIFEST`, `BAMLIST_QCPASS`, `SAMPLE_LIST`,
`OUTBASE`, `THREADS`. Module-specific overrides (e.g. `PI_WIN`,
`PI_STEP`, `MIN_BAND_N` for module 03) are documented in each module's
README.

## Sibling repos

- **catfish-variant-analysis** — produces the QC-passed BAMs consumed
  here (set via `${VARIANT_REPO}` in [00_config.sh](00_config.sh)).
- **catfish-population-analysis** — produces the BEAGLE genotype
  likelihoods and ancestry labels consumed by module 02 and the report
  (set via `${POPULATION_REPO}`).
- **catfish-inversion-analysis** — produces the
  `inversion_candidates.tsv` and `karyotype_calls.tsv` consumed by
  module 03.

## Methods documentation

- [docs/methods/MODULE_3_methods.md](docs/methods/MODULE_3_methods.md)
  — full Methods text for het / ROH / FROH (paper-ready).
- [docs/methods/theta_pi_scaling.md](docs/methods/theta_pi_scaling.md)
  — explains why ANGSD's `tP` must be divided by `nSites` to be
  comparable across windows or samples.
- [docs/methods/window_diversity_texture.md](docs/methods/window_diversity_texture.md)
  — paper-ready Methods for H_chr, DDI, and χ_min.
- [docs/methods/roh_metrics_and_genes.md](docs/methods/roh_metrics_and_genes.md)
  — paper-ready Methods for N_ROH, S_ROH, F_HOM, and genes-in-ROH.
