#!/usr/bin/env bash
# =============================================================================
# STEP_A06_window_H_and_DDI.sh
# =============================================================================
# Per-sample windowed heterozygosity texture metrics:
#
#   H_w(sample)   = tP / nSites           per non-overlapping window
#   DDI(sample)   = MAD(H_w) / median(H_w)        within-genome dispersion
#   chi_min(samp) = min_w [ H_w / median_cohort(H_w) ]
#                   (cohort denominator per-window, optionally smoothed
#                    over ±DDI_SMOOTH_WIN/2 neighbours)
#
# Primary scale is non-overlapping (default 50 kb / 50 kb, configurable via
# DDI_WIN / DDI_STEP in 00_config.sh). If a matching pestPG for the sample
# does not already exist, this script generates it from the existing
# .thetas.idx (no SAF/SFS recomputation).
#
# Outputs (under ${OUT_WINDIV}/02_window_metrics/win<W>_step<S>/):
#   per_sample_window_metrics.tsv    sample, median_H_w, mad_H_w, DDI,
#                                    chi_min, chi_min_chrom, chi_min_pos,
#                                    n_windows_kept, n_windows_dropped
#   cohort_window_H_matrix.tsv.gz    long (sample, chrom, win_start,
#                                    n_sites, H_w); compressed
#   cohort_window_median.tsv         chrom, win_start, n_samples,
#                                    median_H_w_cohort, median_H_w_smoothed
#
# Usage:
#   bash STEP_A06_window_H_and_DDI.sh                  # default DDI_WIN
#   DDI_WIN=10000 DDI_STEP=10000 bash STEP_A06_..sh    # 10 kb sweep
# =============================================================================
set -euo pipefail
MODDIR="$(cd "$(dirname "$0")" && pwd)"
source "$(cd "${MODDIR}/.." && pwd)/00_config.sh"
hr_init_dirs

# ── Module-specific config (with sane defaults) ──────────────────────────────
# Primary scale for DDI / chi_min. Non-overlapping windows are required for
# the metric to make statistical sense (independent observations for MAD).
export DDI_WIN="${DDI_WIN:-50000}"
export DDI_STEP="${DDI_STEP:-50000}"

# Minimum callable sites per window (as fraction of nominal window size).
# Default 0.3 = drop windows with <15 kb callable at 50 kb.
export DDI_MIN_CALLABLE_FRAC="${DDI_MIN_CALLABLE_FRAC:-0.3}"

# Cohort-median smoothing window for chi_min denominator (number of
# windows, must be odd). Default 11 = ±5 neighbours. Set to 1 to disable.
export DDI_SMOOTH_WIN="${DDI_SMOOTH_WIN:-11}"

# MAD constant (R default 1.4826 is Gaussian-consistent; for highly
# non-Gaussian H distributions raw MAD with constant=1 is more honest).
# DDI is a ratio so the constant cancels; we expose this for clarity.
export DDI_MAD_CONSTANT="${DDI_MAD_CONSTANT:-1.4826}"

export OUT_WINDIV="${OUTBASE}/11_window_diversity_texture"
export OUT_WINMET="${OUT_WINDIV}/02_window_metrics/win${DDI_WIN}_step${DDI_STEP}"
mkdir -p "${OUT_WINMET}" "${OUT_WINDIV}/logs"

hr_log "=== Step A06: window H + DDI + chi_min ==="
hr_log "  WIN/STEP        = ${DDI_WIN} / ${DDI_STEP}"
hr_log "  min callable    = ${DDI_MIN_CALLABLE_FRAC} of WIN"
hr_log "  smooth window   = ${DDI_SMOOTH_WIN} windows"
hr_log "  MAD constant    = ${DDI_MAD_CONSTANT}"

hr_check_file "${SAMPLE_LIST}" "sample list"
hr_check_cmd thetaStat
hr_check_cmd python3

# ── Ensure per-sample pestPG at DDI_WIN/DDI_STEP exists ──────────────────────
PESTPG_DIR="${DIR_HET}/03_theta"
EXTRA_PESTPG_DIR="${OUT_WINDIV}/02_window_metrics/pestPG_extra/win${DDI_WIN}_step${DDI_STEP}"
mkdir -p "${EXTRA_PESTPG_DIR}"

ensure_pestpg() {
  local SAMPLE="$1"
  local THETAS="${PESTPG_DIR}/${SAMPLE}.thetas.idx"

  # Look for an existing matching pestPG first (main or multiscale)
  for candidate in \
      "${PESTPG_DIR}/${SAMPLE}.win${DDI_WIN}.step${DDI_STEP}.pestPG" \
      "${PESTPG_DIR}/multiscale/${SAMPLE}.win${DDI_WIN}.step${DDI_STEP}.pestPG" \
      "${EXTRA_PESTPG_DIR}/${SAMPLE}.win${DDI_WIN}.step${DDI_STEP}.pestPG"; do
    if [[ -s "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done

  # None found — compute into EXTRA_PESTPG_DIR
  if [[ ! -s "${THETAS}" ]]; then
    hr_err "Missing thetas.idx for ${SAMPLE}: ${THETAS}"
    return 1
  fi
  hr_log "  ${SAMPLE}: generating pestPG at ${DDI_WIN}/${DDI_STEP}"
  thetaStat do_stat "${THETAS}" \
      -win "${DDI_WIN}" -step "${DDI_STEP}" -type 2 \
      -outnames "${EXTRA_PESTPG_DIR}/${SAMPLE}.win${DDI_WIN}.step${DDI_STEP}" \
      >> "${OUT_WINDIV}/logs/${SAMPLE}.thetaStat.log" 2>&1
  echo "${EXTRA_PESTPG_DIR}/${SAMPLE}.win${DDI_WIN}.step${DDI_STEP}.pestPG"
}

PESTPG_MANIFEST="${OUT_WINMET}/pestPG_manifest.tsv"
: > "${PESTPG_MANIFEST}"
echo -e "sample\tpestPG_path" >> "${PESTPG_MANIFEST}"

while read -r SAMPLE; do
  [[ -z "${SAMPLE}" || "${SAMPLE}" =~ ^# ]] && continue
  PATH_PESTPG="$(ensure_pestpg "${SAMPLE}")" || continue
  echo -e "${SAMPLE}\t${PATH_PESTPG}" >> "${PESTPG_MANIFEST}"
done < "${SAMPLE_LIST}"

hr_log "pestPG manifest: ${PESTPG_MANIFEST} ($(($(wc -l < ${PESTPG_MANIFEST}) - 1)) samples)"

# ── Compute window metrics (two-pass: cohort medians → per-sample DDI) ────
hr_log "Computing per-window H, cohort medians, DDI, chi_min..."
python3 "${MODDIR}/compute_window_metrics.py" \
  --pestpg-manifest "${PESTPG_MANIFEST}" \
  --win "${DDI_WIN}" \
  --min-callable-frac "${DDI_MIN_CALLABLE_FRAC}" \
  --smooth-window "${DDI_SMOOTH_WIN}" \
  --mad-constant "${DDI_MAD_CONSTANT}" \
  --out-per-sample "${OUT_WINMET}/per_sample_window_metrics.tsv" \
  --out-matrix "${OUT_WINMET}/cohort_window_H_matrix.tsv.gz" \
  --out-cohort-median "${OUT_WINMET}/cohort_window_median.tsv"

hr_log "Outputs:"
hr_log "  ${OUT_WINMET}/per_sample_window_metrics.tsv"
hr_log "  ${OUT_WINMET}/cohort_window_H_matrix.tsv.gz"
hr_log "  ${OUT_WINMET}/cohort_window_median.tsv"
hr_log "=== Step A06 complete ==="
