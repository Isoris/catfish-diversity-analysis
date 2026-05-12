#!/usr/bin/env bash
# =============================================================================
# STEP_A05_per_chromosome_heterozygosity.sh
# =============================================================================
# Aggregate per-site theta_pi from each sample's .thetas.idx (produced by
# STEP_A02_run_heterozygosity.sh) into a per-chromosome heterozygosity table.
#
#   H_chr(sample) = sum_over_sites_in_chr( theta_pi_site ) / N_callable_chr
#
# theta_pi_site is recovered by exp() of the natural-log values emitted by
# `thetaStat print`. N_callable_chr is read from a cohort-level cache built
# from CALLABLE_SITES (column 1 line counts per chromosome).
#
# Outputs (under ${OUT_WINDIV}/01_per_chrom_H/):
#   per_sample_per_chromosome_H.tsv   long: sample, chrom, n_sites_obs,
#                                     n_callable_chr, sum_theta_pi, H_chr
#   cohort_sample_x_chrom_H.tsv       wide matrix: rows=sample, cols=chrom
#   cohort_per_chromosome_summary.tsv per-chrom median, IQR, range, n_outliers
#   callable_sites_per_chrom.tsv      cache of N_callable per chromosome
#
# Usage:
#   bash STEP_A05_per_chromosome_heterozygosity.sh            # all samples
#   bash STEP_A05_per_chromosome_heterozygosity.sh CGA009     # single sample
# =============================================================================
set -euo pipefail
MODDIR="$(cd "$(dirname "$0")" && pwd)"
source "$(cd "${MODDIR}/.." && pwd)/00_config.sh"
hr_init_dirs

SINGLE_SAMPLE="${1:-}"

export OUT_WINDIV="${OUTBASE}/11_window_diversity_texture"
export OUT_PERCHROM="${OUT_WINDIV}/01_per_chrom_H"
mkdir -p "${OUT_PERCHROM}" "${OUT_WINDIV}/logs"

hr_log "=== Step A05: per-chromosome heterozygosity (H_chr) ==="
hr_check_file "${SAMPLE_LIST}" "sample list"
hr_check_file "${CALLABLE_SITES}" "ANGSD callable sites"
hr_check_cmd thetaStat
hr_check_cmd awk
hr_check_cmd python3

# ── Build / reuse cohort-level cache: callable sites per chromosome ────────
CALLABLE_CACHE="${OUT_PERCHROM}/callable_sites_per_chrom.tsv"
if [[ ! -s "${CALLABLE_CACHE}" ]]; then
  hr_log "Building callable-sites cache from ${CALLABLE_SITES}"
  # Use default awk FS so the cache is robust to tab- or space-separated
  # ANGSD sites files. Chromosome is column 1 in both layouts.
  {
    echo -e "chrom\tn_callable"
    awk 'BEGIN{OFS="\t"} $1 !~ /^#/ {n[$1]++} END{
      for (c in n) print c, n[c];
    }' "${CALLABLE_SITES}" | sort -k1,1
  } > "${CALLABLE_CACHE}"
  hr_log "  cache: ${CALLABLE_CACHE}"
else
  hr_log "Using existing callable-sites cache: ${CALLABLE_CACHE}"
fi

# ── Per-sample aggregation via thetaStat print | awk → per-chrom rows ─────
PER_SAMPLE_DIR="${OUT_PERCHROM}/by_sample"
mkdir -p "${PER_SAMPLE_DIR}"

run_one_sample() {
  local SAMPLE="$1"
  local THETAS="${DIR_HET}/03_theta/${SAMPLE}.thetas.idx"
  local OUT="${PER_SAMPLE_DIR}/${SAMPLE}.per_chrom_tP.tsv"

  if [[ ! -s "${THETAS}" ]]; then
    hr_err "Missing thetas.idx for ${SAMPLE}: ${THETAS}"
    return 1
  fi
  if [[ -s "${OUT}" ]]; then
    hr_log "  ${SAMPLE} already aggregated, skipping (delete ${OUT} to rerun)"
    return 0
  fi

  hr_log "  ${SAMPLE}: aggregating per-site tP per chromosome"
  # thetaStat print emits log-space per-site theta estimates:
  #   chr  pos  log_tW  log_tP  log_tF  log_tH  log_tL
  # We sum exp(log_tP) and count sites, grouped by chrom. The trailing END
  # block emits the final chromosome.
  thetaStat print "${THETAS}" 2>> "${OUT_WINDIV}/logs/${SAMPLE}.thetaprint.log" \
    | awk -v OFS='\t' -v sample="${SAMPLE}" '
        BEGIN { prev=""; n=0; s=0; print "sample","chrom","n_sites_obs","sum_theta_pi" }
        /^#/ || /^chr\t/ || /^Chr\t/ { next }
        {
          if ($1 != prev) {
            if (prev != "") print sample, prev, n, s;
            prev=$1; n=0; s=0;
          }
          n++;
          s += exp($4);
        }
        END {
          if (prev != "") print sample, prev, n, s;
        }
      ' > "${OUT}.tmp" && mv "${OUT}.tmp" "${OUT}"
}

if [[ -n "${SINGLE_SAMPLE}" ]]; then
  run_one_sample "${SINGLE_SAMPLE}"
else
  while read -r SAMPLE; do
    [[ -z "${SAMPLE}" || "${SAMPLE}" =~ ^# ]] && continue
    run_one_sample "${SAMPLE}"
  done < "${SAMPLE_LIST}"
fi

# ── Consolidate to long, wide, and per-chrom summary ──────────────────────
hr_log "Consolidating per-sample tables → long / wide / per-chrom summary"
python3 "${MODDIR}/aggregate_per_chrom_H.py" \
  --per-sample-dir "${PER_SAMPLE_DIR}" \
  --callable-cache "${CALLABLE_CACHE}" \
  --out-long "${OUT_PERCHROM}/per_sample_per_chromosome_H.tsv" \
  --out-wide "${OUT_PERCHROM}/cohort_sample_x_chrom_H.tsv" \
  --out-summary "${OUT_PERCHROM}/cohort_per_chromosome_summary.tsv"

hr_log "Outputs:"
hr_log "  ${OUT_PERCHROM}/per_sample_per_chromosome_H.tsv"
hr_log "  ${OUT_PERCHROM}/cohort_sample_x_chrom_H.tsv"
hr_log "  ${OUT_PERCHROM}/cohort_per_chromosome_summary.tsv"
hr_log "=== Step A05 complete ==="
