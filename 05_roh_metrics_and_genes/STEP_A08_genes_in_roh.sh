#!/usr/bin/env bash
# =============================================================================
# STEP_A08_genes_in_roh.sh
# =============================================================================
# Annotate per-sample ROH tracts with overlapping genes from a GFF/GFF3,
# producing both per-sample views and cohort-wide recurrence / private-ROH
# gene tables (analogous to the "80 genes in wild private ROHs" use case).
#
# Inputs
#   ${GFF_GENE_ANNOT}            GFF/GFF3 (gzip OK) — source of gene features
#   ${DIR_ROH}/roh_tracts_all.bed   per-sample ROH tracts from STEP_A04
#
# Outputs (under ${OUT_ROHEXT}/genes_in_roh/)
#   genes.bed                                       BED of gene features (cached)
#   intersections/roh_x_genes.tsv                   raw bedtools intersect
#   per_sample/<SAMPLE>.genes_in_roh.tsv            one TSV per sample
#   per_sample_genes_in_roh_summary.tsv             sample-level summary
#   cohort_genes_in_roh.tsv                         long: (gene, sample) pairs
#   cohort_gene_recurrence.tsv                      per-gene cohort recurrence
#   cohort_private_roh_genes.tsv                    genes recurrent in only 1
#                                                   sample (the "private ROH"
#                                                   set)
#
# Config (read from 00_config.sh; override via environment):
#   GFF_GENE_ANNOT         path to GFF/GFF3 (required)
#   GENE_FEATURE_TYPES     comma-separated GFF type column to keep
#                          (default "gene"; e.g. "gene,pseudogene")
#   GENE_NAME_KEYS         comma-separated GFF attribute keys for "name"
#                          (default "Name,gene_name,gene"; falls back to ID)
#   ROH_GENE_MIN_OVERLAP_BP minimum bp overlap to record (default 1)
#
# Usage:
#   bash STEP_A08_genes_in_roh.sh
# =============================================================================
set -euo pipefail
MODDIR="$(cd "$(dirname "$0")" && pwd)"
source "$(cd "${MODDIR}/.." && pwd)/00_config.sh"
hr_init_dirs

export OUT_ROHEXT="${OUTBASE}/12_roh_metrics_and_genes"
export GENES_DIR="${OUT_ROHEXT}/genes_in_roh"
mkdir -p "${GENES_DIR}" \
         "${GENES_DIR}/intersections" \
         "${GENES_DIR}/per_sample" \
         "${OUT_ROHEXT}/logs"

# Config knobs with safe defaults
export GFF_GENE_ANNOT="${GFF_GENE_ANNOT:-}"
export GENE_FEATURE_TYPES="${GENE_FEATURE_TYPES:-gene}"
export GENE_NAME_KEYS="${GENE_NAME_KEYS:-Name,gene_name,gene}"
export ROH_GENE_MIN_OVERLAP_BP="${ROH_GENE_MIN_OVERLAP_BP:-1}"

ROH_BED="${DIR_ROH}/roh_tracts_all.bed"
GENES_BED="${GENES_DIR}/genes.bed"
INTERSECT_TSV="${GENES_DIR}/intersections/roh_x_genes.tsv"

hr_log "=== Step A08: genes-in-ROH annotation ==="
hr_check_file "${ROH_BED}" "ROH tracts BED"
if [[ -z "${GFF_GENE_ANNOT}" || ! -f "${GFF_GENE_ANNOT}" ]]; then
  hr_die "GFF_GENE_ANNOT not set or missing: '${GFF_GENE_ANNOT}'. Export it
or add to 00_config.sh, e.g.:
  export GFF_GENE_ANNOT=\"\${BASE}/00-samples/fClaHyb_Gar_LG.gff3\""
fi
hr_check_cmd python3
hr_check_cmd bedtools
hr_check_cmd sort

# ── Build/reuse cached gene BED from GFF ─────────────────────────────────────
if [[ -s "${GENES_BED}" ]]; then
  hr_log "Using cached gene BED: ${GENES_BED}"
else
  hr_log "Parsing GFF → BED  (types=${GENE_FEATURE_TYPES})"
  python3 "${MODDIR}/genes_in_roh.py" gff-to-bed \
    --gff "${GFF_GENE_ANNOT}" \
    --feature-types "${GENE_FEATURE_TYPES}" \
    --name-keys "${GENE_NAME_KEYS}" \
    --out-bed "${GENES_BED}"
fi
hr_log "  $(wc -l < "${GENES_BED}") gene features"

# ── bedtools intersect: ROH tracts × genes ───────────────────────────────────
hr_log "bedtools intersect: ROH × genes"
# -wa -wb keeps both sides of each overlap; -sorted requires both inputs sorted
# by chrom. Both are already sorted (BED produced by genes_in_roh.py + the ROH
# BED from convert_ibd.pl), but we re-sort defensively into temp files.
ROH_SORTED="${GENES_DIR}/intersections/roh_sorted.bed"
GENES_SORTED="${GENES_DIR}/intersections/genes_sorted.bed"
sort -k1,1 -k2,2n "${ROH_BED}" > "${ROH_SORTED}"
sort -k1,1 -k2,2n "${GENES_BED}" > "${GENES_SORTED}"

# Output columns:
#   roh_chr roh_start roh_end roh_sample roh_length
#   gene_chr gene_start gene_end gene_id gene_name strand
bedtools intersect \
  -a "${ROH_SORTED}" \
  -b "${GENES_SORTED}" \
  -wa -wb -sorted \
  > "${INTERSECT_TSV}.tmp"
mv "${INTERSECT_TSV}.tmp" "${INTERSECT_TSV}"
hr_log "  $(wc -l < "${INTERSECT_TSV}") (ROH × gene) overlap rows"

# ── Aggregate into per-sample + cohort tables ────────────────────────────────
hr_log "Aggregating per-sample + cohort gene tables"
python3 "${MODDIR}/genes_in_roh.py" aggregate \
  --intersection "${INTERSECT_TSV}" \
  --min-overlap-bp "${ROH_GENE_MIN_OVERLAP_BP}" \
  --per-sample-dir "${GENES_DIR}/per_sample" \
  --out-per-sample-summary "${GENES_DIR}/per_sample_genes_in_roh_summary.tsv" \
  --out-cohort-long "${GENES_DIR}/cohort_genes_in_roh.tsv" \
  --out-cohort-recurrence "${GENES_DIR}/cohort_gene_recurrence.tsv" \
  --out-private "${GENES_DIR}/cohort_private_roh_genes.tsv"

hr_log "Outputs:"
hr_log "  ${GENES_DIR}/genes.bed                               (cached gene BED)"
hr_log "  ${INTERSECT_TSV}                                     (raw intersect)"
hr_log "  ${GENES_DIR}/per_sample/<SAMPLE>.genes_in_roh.tsv    (per-sample)"
hr_log "  ${GENES_DIR}/per_sample_genes_in_roh_summary.tsv     (sample-level)"
hr_log "  ${GENES_DIR}/cohort_genes_in_roh.tsv                 (gene × sample)"
hr_log "  ${GENES_DIR}/cohort_gene_recurrence.tsv              (per-gene)"
hr_log "  ${GENES_DIR}/cohort_private_roh_genes.tsv            (private to 1)"
hr_log "=== Step A08 complete ==="
