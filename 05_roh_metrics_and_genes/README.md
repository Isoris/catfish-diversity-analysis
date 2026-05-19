# 05_roh_metrics_and_genes — ROH-derived metrics + genes-in-ROH

Two analyses built on top of the ROH calls from [02_roh/](../02_roh/):

1. **STEP_A07 — derived per-sample ROH metrics**: counts (N_ROH), sums
   (S_ROH), the existing F_ROH from the parser, and a cohort-relative
   F_HOM (deficit of observed heterozygosity).
2. **STEP_A08 — genes-in-ROH annotation**: overlap each sample's ROH
   tracts with gene features parsed from a GFF/GFF3, then emit
   per-sample, cohort-recurrence, and "private ROH gene" tables
   (genes that overlap an ROH in exactly one sample).

Operates on the BED + summary that STEP_A04 already writes; nothing
upstream needs to be re-run.

## Metric definitions

### S_ROH

Per-sample sum of ROH tract lengths, in bp and Mb, total and per length
bin. Bins use the same thresholds as STEP_A04 / `MODULE_3_methods.md`:

| Bin | Range |
|------|-------|
| short  | [`ROH_MIN_LEN_BP`, 1 Mb) |
| medium | [1 Mb, 5 Mb) |
| long   | [5 Mb, ∞) |

### N_ROH

Tract count per sample, total and per length bin (same bins as S_ROH).

### F_HOM

Cohort-relative inbreeding deficit, computed from per-sample observed
heterozygosity `H_obs` (from the SFS in `genomewide_heterozygosity.tsv`):

```
F_HOM_mean   = 1 - H_obs / mean   (H_obs across cohort)
F_HOM_median = 1 - H_obs / median (H_obs across cohort)
```

Both centres are emitted so downstream choice between robust (median)
and conventional (mean) is explicit. Note this is a cohort-relative
estimator, not the classical PLINK F_HW (which would require HW
expected genotypes from cohort allele frequencies — a follow-up via
`-doMaf` on the BEAGLE GLs if needed).

## Scripts

| Step | Script | Purpose |
|------|--------|---------|
| A07 | `STEP_A07_roh_derived_metrics.sh` | driver for N_ROH/S_ROH/F_HOM |
| A07 helper | `compute_roh_derived_metrics.py` | joins ROH BED + ROH/H summaries |
| A08 | `STEP_A08_genes_in_roh.sh` | driver for GFF gene overlap |
| A08 helper | `genes_in_roh.py` | `gff-to-bed` and `aggregate` subcommands |

`STEP_A08` shells out to `bedtools intersect` once over
`roh_tracts_all.bed × genes.bed`, then aggregates in Python.

## Outputs

Under `${OUT_ROHEXT}` = `${OUTBASE}/12_roh_metrics_and_genes/`:

```
per_sample_roh_derived.tsv         A07: per-sample N_ROH/S_ROH/F_HOM/F_ROH
cohort_roh_metric_summary.tsv      A07: distribution summary of each metric

genes_in_roh/
  genes.bed                                       cached gene BED
  intersections/roh_x_genes.tsv                   raw bedtools output
  per_sample/<SAMPLE>.genes_in_roh.tsv            per-sample view
  per_sample_genes_in_roh_summary.tsv             one row per sample
  cohort_genes_in_roh.tsv                         long (gene × sample)
  cohort_gene_recurrence.tsv                      per-gene cohort recurrence
  cohort_private_roh_genes.tsv                    genes recurrent in 1 sample
                                                  (analogue of "80 wild-
                                                  private-ROH genes")
```

`per_sample/<SAMPLE>.genes_in_roh.tsv` columns:

```
gene_id  gene_name  chrom  gene_start  gene_end  strand
roh_start  roh_end  roh_length
overlap_bp  overlap_frac_of_gene  fully_contained
```

`cohort_gene_recurrence.tsv` columns:

```
gene_id  gene_name  gene_chr  gene_start  gene_end  strand
n_samples_with_overlap  samples  n_fully_contained
total_overlap_bp  frac_cohort
```

## Configuration

All flags read from `00_config.sh` (override via environment):

| Variable | Default | Meaning |
|----------|---------|---------|
| `GFF_GENE_ANNOT` | `${BASE}/00-samples/fClaHyb_Gar_LG.gff3` | path to GFF/GFF3 (gzip OK) |
| `GENE_FEATURE_TYPES` | `gene` | comma-separated GFF type filter (e.g. `gene,pseudogene`) |
| `GENE_NAME_KEYS` | `Name,gene_name,gene` | attribute keys tried in order for human-readable name |
| `ROH_GENE_MIN_OVERLAP_BP` | `1` | min bp of overlap to record |
| `ROH_MIN_LEN_BP` | `100000` (inherited) | min tract length for counts/bins |

## Dependencies

- Upstream: `02_roh/STEP_A04_parse_roh_and_het.sh` (BED + summary).
- Tools: `bedtools`, `sort`, `python3` with `numpy` + `pandas`.

## Notes

- `genes.bed` is cached after first GFF parse; delete it to rebuild
  (e.g. after changing `GENE_FEATURE_TYPES` or `GENE_NAME_KEYS`).
- "Private-ROH genes" are defined per-cohort: a gene whose overlap
  occurs in exactly one cohort sample. This is the operational form
  of the table referenced in the manuscript ("80 genes in wild
  private ROHs"); within the 226-sample hatchery cohort, this
  identifies broodstock-specific ROH-burdened genes.
- "Fully contained" means the gene's entire span lies inside a single
  ROH tract; this is stricter than `overlap_frac_of_gene = 1.0` only
  when a gene crosses two abutting tracts (rare).
- `F_HOM` here is cohort-relative because the upstream pipeline uses
  ANGSD SAF/SFS and does not maintain a global HW-expected H. If
  the classical F_HW is needed downstream, run `-doMaf 1` on the
  cohort BEAGLE GLs to obtain per-site `2pq`, sum across callable
  sites, divide H_obs by that. The current implementation is
  intentionally honest about what it computes.
