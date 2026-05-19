#!/usr/bin/env bash
# =============================================================================
# STEP_A07_roh_derived_metrics.sh
# =============================================================================
# Per-sample ROH-derived metrics. Builds on the BED + summary produced by
# STEP_A04_parse_roh_and_het.sh.
#
#   S_ROH        = sum of ROH tract lengths (bp; also reported in Mb)
#                  per-sample, all tracts and per length bin
#   N_ROH        = count of ROH tracts (all and per length bin)
#   F_ROH        = S_ROH / callable_nonrepeat_bp   (from STEP_A04 summary)
#   F_HOM        = 1 - H_obs / cohort_centre(H_obs)
#                  reported for both mean(H_obs) and median(H_obs) centres
#
# Length bins follow STEP_A04 / docs/methods/MODULE_3_methods.md:
#   short  = [ROH_MIN_LEN_BP, 1 Mb)
#   medium = [1 Mb, 5 Mb)
#   long   = [5 Mb, ∞)
#
# Output (under ${OUT_ROHEXT}/):
#   per_sample_roh_derived.tsv    one row per sample (NROH, SROH, FROH, FHOM)
#   cohort_roh_metric_summary.tsv distribution summary of each metric
#
# Usage:
#   bash STEP_A07_roh_derived_metrics.sh
# =============================================================================
set -euo pipefail
MODDIR="$(cd "$(dirname "$0")" && pwd)"
source "$(cd "${MODDIR}/.." && pwd)/00_config.sh"
hr_init_dirs

export OUT_ROHEXT="${OUTBASE}/12_roh_metrics_and_genes"
mkdir -p "${OUT_ROHEXT}" "${OUT_ROHEXT}/logs"

ROH_BED="${DIR_ROH}/roh_tracts_all.bed"
ROH_SUMMARY="${DIR_ROH}/catfish_roh.per_sample_roh.tsv"
HET_SUMMARY="${DIR_HET}/04_summary/genomewide_heterozygosity.tsv"

hr_log "=== Step A07: ROH-derived metrics (N_ROH, S_ROH, F_HOM) ==="
hr_check_file "${ROH_BED}" "ROH tracts BED"
hr_check_file "${ROH_SUMMARY}" "per-sample ROH summary from STEP_A04"
hr_check_file "${HET_SUMMARY}" "genome-wide H summary from STEP_A02"
hr_check_cmd python3

python3 "${MODDIR}/compute_roh_derived_metrics.py" \
  --roh-bed "${ROH_BED}" \
  --roh-summary "${ROH_SUMMARY}" \
  --het-summary "${HET_SUMMARY}" \
  --min-roh-bp "${ROH_MIN_LEN_BP}" \
  --bin-medium-bp 1000000 \
  --bin-long-bp 5000000 \
  --out-per-sample "${OUT_ROHEXT}/per_sample_roh_derived.tsv" \
  --out-summary "${OUT_ROHEXT}/cohort_roh_metric_summary.tsv"

hr_log "Outputs:"
hr_log "  ${OUT_ROHEXT}/per_sample_roh_derived.tsv"
hr_log "  ${OUT_ROHEXT}/cohort_roh_metric_summary.tsv"
hr_log "=== Step A07 complete ==="
