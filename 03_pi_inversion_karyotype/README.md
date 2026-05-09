# 03_pi_inversion_karyotype — Band-level π11 / π22 / π12

Per-inversion nucleotide diversity within and between the two
homokaryotype bands of each candidate inversion, plus an optional
per-sample θπ rollup by band.

## What this is (and is not)

For each inversion candidate with a stable karyotype call:

- **π11** — within-group diversity of `Homo_1` arrangement
- **π22** — within-group diversity of `Homo_2` arrangement
- **π12** — pairwise divergence between `Homo_1` and `Homo_2`,
  computed from allele frequencies in the two homokaryotype groups.
  **π12 is not the heterozygote-group diversity.** It does not use
  `Het` samples at all.

Sitewise formulas, where `p1` and `p2` are the frequency of the
**same** allele in `Homo_1` and `Homo_2`:

```
π11 = 2 * p1 * (1 - p1)
π22 = 2 * p2 * (1 - p2)
π12 = p1 * (1 - p2) + p2 * (1 - p1)
```

Naming policy: pre-polarization labels are `π11 / π22 / π12`. Only after
independent outgroup polarization should they be renamed to the
paper-style `πAA / πDD / πAD` (Ancestral / Derived).

## Why a shared 4-column sites file matters

`p1` and `p2` must be frequencies of the **same** allele. If ANGSD is
allowed to infer major/minor independently per group, on a balanced
site `p1` and `p2` may end up referring to opposite alleles, and the π12
formula silently breaks.

This module enforces shared coding by:

1. Running ANGSD on the **combined Homo_1 ∪ Homo_2 BAM list** in the
   inversion interval to infer major/minor.
2. Writing a 4-column sites file: `chr  pos  major  minor`.
3. `angsd sites index` on it.
4. Running ANGSD per band with `-doMajorMinor 3 -sites <file>`, which
   forces ANGSD to use the provided major/minor and skip per-group
   inference.

## Site-basis caveat (read this)

By default the shared sites file contains **only the variant sites**
that pass ANGSD's SNP-calling filter on the combined Homo_1 + Homo_2
sample (`-SNP_pval 1e-6 -doMaf 1`). That means:

- The reported π values are **per-variant-site**, not per-callable-site.
- Absolute values are **inflated** relative to the genome-wide θπ from
  module 01 (which is per-callable-site).
- Comparisons across bands and across windows **within this module**
  are still meaningful — they share the same denominator.

To get unbiased per-callable-site π, the caller must instead build a
shared sites file that includes monomorphic callable sites (which
contribute 0 to all three π). The output files carry a
`site_basis` column (`variant_only` or `all_callable`) so downstream
code can rescale or flag.

## Pipeline

```
inversion_candidates.tsv  +  karyotype_calls.tsv  +  manifest
  │
  └─ STEP_PI01_make_karyotype_group_theta_pi.sh    per-inversion driver
       │
       ├─ scripts/make_karyotype_bamlists.R         build Homo_1/Homo_2/Het BAM lists
       ├─ ANGSD on Homo_1 ∪ Homo_2                  shared 4-col sites file
       ├─ ANGSD per band with -doMajorMinor 3       per-band MAFs on shared sites
       ├─ scripts/compute_pi11_pi22_pi12.R          sitewise + windowed π
       ├─ scripts/average_sample_pestPG_by_karyotype.R   optional mean tP by band
       └─ scripts/plot_band_pi_panel.R              two-panel PDF
```

## Inputs

- `${INV_CANDIDATES}` (default `${DIR_INV}/inversion_candidates.tsv`)
  Header: `inversion_id  chr  start  end`
- `${KARYOTYPE_CALLS}` (default `${DIR_INV}/karyotype_calls.tsv`)
  Header: `sample  inversion_id  group`  with `group ∈ {Homo_1, Het, Homo_2}`
- `${SAMPLE_MANIFEST}` — same manifest used by STEP_A02; column 3 (filtered BAM).
- `${REF}`, `${CALLABLE_SITES}` — shared with the rest of the repo.
- `${PESTPG_DIR}` (optional, default `${DIR_HET}/03_theta`) — per-sample
  pestPG files from STEP_A02.

## Module-specific config (overridable)

| Variable          | Default                                | Meaning                             |
|-------------------|----------------------------------------|-------------------------------------|
| `INV_CANDIDATES`  | `${DIR_INV}/inversion_candidates.tsv`  | Inversion table                     |
| `KARYOTYPE_CALLS` | `${DIR_INV}/karyotype_calls.tsv`       | Per-sample karyotype calls          |
| `PESTPG_DIR`      | `${DIR_HET}/03_theta`                  | Per-sample pestPG location          |
| `MIN_BAND_N`      | `5`                                    | Skip inversion if Homo_1 or Homo_2 below this |
| `PI_WIN`          | `50000`                                | Window size for band-π (bp)         |
| `PI_STEP`         | `10000`                                | Step size (bp), supports overlap    |

`PI_WIN` / `PI_STEP` are independent of the genome-wide `WIN` / `STEP`
in `00_config.sh`; inversion intervals warrant a finer scale.

## Outputs

Under `${OUTBASE}/11_band_pi/`:

```
bamlists/
  {INV}.Homo_1.bamlist
  {INV}.Homo_2.bamlist
  {INV}.Het.bamlist
  {INV}.Homo_1_2.shared.bamlist
allele_freq/
  {INV}.shared.mafs.gz                 ← combined run output
  {INV}.shared.majorMinor.sites        ← 4-col indexed sites
  {INV}.shared.majorMinor.sites.idx
  {INV}.shared.majorMinor.sites.bin
  {INV}.Homo_1.mafs.gz                 ← per-band, shared coding
  {INV}.Homo_2.mafs.gz
pi_curves/
  {INV}.pi11_pi22_pi12.sitewise.tsv.gz
  {INV}.pi11_pi22_pi12.window.tsv
mean_tP/
  {INV}.mean_tP_by_band.tsv
plots/
  {INV}.band_pi_panel.pdf
logs/
  {INV}.*.log
```

### Sitewise output columns

`inversion_id  chromo  position  major  minor  p1  p2  pi11  pi22  pi12  site_basis`

### Windowed output columns

`inversion_id  chromo  win_start  win_end  win_center  n_sites  mean_pi11  mean_pi22  mean_pi12  sum_pi11  sum_pi22  sum_pi12  site_basis`

### mean_tP output columns

`inversion_id  chromo  WinStart  WinStop  WinCenter  group  n_samples  mean_tP  median_tP  sd_tP  se_tP`

The mean_tP file reports the **mean of per-sample within-individual
diversity** (per-site `tP / nSites` from each sample's pestPG). This is
**not** formal group π; the column is named `mean_tP` and the file is
under `mean_tP/`, not `pi_curves/`, deliberately.

## Usage

```bash
# All inversions in INV_CANDIDATES
bash 03_pi_inversion_karyotype/STEP_PI01_make_karyotype_group_theta_pi.sh

# A subset
bash 03_pi_inversion_karyotype/STEP_PI01_make_karyotype_group_theta_pi.sh INV01 INV03
```

## Cohort scope

This module operates on the 226-sample pure *C. gariepinus* hatchery
cohort (the same cohort as the rest of `catfish-diversity-analysis`).
It does **not** apply to the F₁ hybrid assembly cohort or any planned
*C. macrocephalus* wild cohort.
