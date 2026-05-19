# ROH derived metrics + genes-in-ROH — Methods

Methods text (paper-ready) for the per-sample ROH burden statistics
and the gene-level overlap analysis produced by
[05_roh_metrics_and_genes/](../../05_roh_metrics_and_genes/).

## Methods

### ROH-derived per-sample metrics (N_ROH, S_ROH, F_ROH, F_HOM)

From the per-sample ROH BED produced by `STEP_A04_parse_roh_and_het.sh`
(ngsF-HMM calls converted to per-sample BED with
`convert_ibd.pl` and length-filtered at the cohort minimum), we
computed the following per-sample summaries:

- **N_ROH**: count of ROH tracts, total and partitioned into short
  (≥ `ROH_MIN_LEN_BP` and < 1 Mb), medium (1–5 Mb), and long (≥ 5 Mb)
  length bins.
- **S_ROH**: sum of ROH tract lengths, total and per length bin, in bp
  and Mb. S_ROH is the unnormalised counterpart of F_ROH.
- **F_ROH**: S_ROH / callable_nonrepeat_bp, with the callable-non-
  repeat denominator taken from `STEP_A04` (genome span retained after
  intersecting the callable mask with the non-repeat mask).
- **F_HOM**: cohort-relative inbreeding deficit defined from observed
  per-sample heterozygosity `H_obs` (from the folded SFS produced in
  `STEP_A02`):

  ```
  F_HOM_mean   = 1 - H_obs(i) / mean   ( H_obs over cohort )
  F_HOM_median = 1 - H_obs(i) / median ( H_obs over cohort )
  ```

  Both centring conventions are reported. F_HOM = 0 corresponds to
  cohort-average heterozygosity; positive values indicate
  heterozygote deficit relative to the cohort, negative values
  indicate elevated heterozygosity relative to the cohort. This
  estimator is honest about what is computed: it is cohort-relative
  and does not assume Hardy–Weinberg expected genotypes (those would
  require cohort allele frequencies from the BEAGLE GLs, a follow-up
  if the classical F_HW is required).

### Genes in ROH (per sample, cohort recurrence, private)

Gene features were parsed from the *C. gariepinus* haplotype
reference GFF/GFF3 (`GFF_GENE_ANNOT`), retaining records with type
matching `GENE_FEATURE_TYPES` (default `gene`) and emitting BED with
gene identifier and human-readable name (`GENE_NAME_KEYS`
preferences, fallback to ID). The all-sample ROH tract BED was
intersected with this gene BED using `bedtools intersect -wa -wb`,
keeping overlaps of ≥ `ROH_GENE_MIN_OVERLAP_BP` bp. For each
(ROH × gene) overlap we recorded the overlap length, the fraction
of the gene covered, and whether the gene span was fully contained
within a single ROH tract.

From this raw intersection we constructed:

- a per-sample table for each individual listing every gene
  overlapping any of its ROH tracts, with the corresponding ROH
  coordinates, overlap length, and overlap fraction;
- a sample-level summary reporting the number of unique genes
  overlapped, the number fully contained, and the total gene bp in
  ROH;
- a cohort gene-by-sample long table;
- a per-gene cohort recurrence table reporting the number of
  samples in which the gene overlaps any ROH, the comma-separated
  sample list, the total cohort-wide overlap bp, and the fraction
  of the cohort affected;
- a "private-ROH gene" table restricted to genes overlapping ROH
  in exactly one cohort sample (analogous to the 80 wild-private-
  ROH genes table in the parent framework — here applied at the
  hatchery-cohort scale to surface broodstock-specific
  ROH-burdened genes).

## Key parameters

| Parameter | Value | Used in |
|-----------|-------|---------|
| Minimum tract length (`ROH_MIN_LEN_BP`) | 100 kb | A07 (binning), A08 (via BED) |
| Short / medium / long bin boundaries | 1 Mb, 5 Mb | A07 |
| H_obs cohort centres reported | mean and median | A07 |
| F_ROH denominator | callable ∩ non-repeat bp | A07 (joined from A04) |
| GFF feature filter | `gene` (configurable) | A08 |
| Gene-name attribute keys | `Name,gene_name,gene` → ID | A08 |
| Minimum gene-overlap bp recorded | 1 bp (configurable) | A08 |
| Gene-fully-contained criterion | gene span ⊆ single ROH tract | A08 |
| Private-ROH gene criterion | gene recurrent in 1 cohort sample | A08 |

## Reported per-sample columns

`per_sample_roh_derived.tsv`:

```
sample  ancestry_label
H_obs  F_HOM_mean  F_HOM_median
N_ROH_total  N_ROH_short  N_ROH_medium  N_ROH_long
S_ROH_total_bp  S_ROH_short_bp  S_ROH_medium_bp  S_ROH_long_bp
S_ROH_total_Mb  S_ROH_long_Mb
F_ROH  callable_nonrepeat_bp
```

`per_sample/<SAMPLE>.genes_in_roh.tsv`:

```
gene_id  gene_name  chrom  gene_start  gene_end  strand
roh_start  roh_end  roh_length
overlap_bp  overlap_frac_of_gene  fully_contained
```
