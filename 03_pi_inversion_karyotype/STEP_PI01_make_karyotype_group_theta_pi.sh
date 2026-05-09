#!/usr/bin/env bash
# =============================================================================
# STEP_PI01_make_karyotype_group_theta_pi.sh
# =============================================================================
# Per-inversion band-level nucleotide diversity / divergence:
#
#   pi11  = within-group diversity of Homo_1   (= 2 * p1 * (1 - p1))
#   pi22  = within-group diversity of Homo_2   (= 2 * p2 * (1 - p2))
#   pi12  = pairwise divergence Homo_1 vs Homo_2
#                                              (= p1*(1-p2) + p2*(1-p1))
#
# IMPORTANT: pi12 is NOT the heterozygote group. It is computed from
# allele frequencies in Homo_1 and Homo_2 separately, with allele coding
# fixed across both groups via a shared 4-column sites file.
#
# Naming policy:
#   - Pre-polarization (this module): pi11 / pi22 / pi12
#   - Post-polarization (downstream, with outgroup):  piAA / piDD / piAD
#
# Optionally also averages existing per-sample .pestPG tracks by karyotype
# group as "mean_per_sample_tP" (this is *not* group pi; label accordingly).
#
# Usage:
#   bash STEP_PI01_make_karyotype_group_theta_pi.sh                  # all inversions
#   bash STEP_PI01_make_karyotype_group_theta_pi.sh INV01 INV03      # subset
# =============================================================================
set -euo pipefail

MODDIR="$(cd "$(dirname "$0")" && pwd)"
source "$(cd "${MODDIR}/.." && pwd)/00_config.sh"
hr_init_dirs

# ── Module-specific config (with sane defaults) ──────────────────────────────
# Inversion candidates and karyotype calls live alongside the inversion repo;
# default to the het_roh inversion-support area but allow override.
export INV_CANDIDATES="${INV_CANDIDATES:-${DIR_INV}/inversion_candidates.tsv}"
export KARYOTYPE_CALLS="${KARYOTYPE_CALLS:-${DIR_INV}/karyotype_calls.tsv}"
# Per-sample pestPG directory (from STEP_A02). Optional.
export PESTPG_DIR="${PESTPG_DIR:-${DIR_HET}/03_theta}"

# Min samples per band (skip inversion if either Homo_1 or Homo_2 below this)
export MIN_BAND_N="${MIN_BAND_N:-5}"

# Window/step for band-pi (independent of genome-wide WIN/STEP from config)
# Default to a finer scale appropriate for inversion intervals.
export PI_WIN="${PI_WIN:-50000}"
export PI_STEP="${PI_STEP:-10000}"

# ── Module output tree ───────────────────────────────────────────────────────
export OUT_BANDPI="${OUTBASE}/11_band_pi"
mkdir -p "${OUT_BANDPI}/bamlists" \
         "${OUT_BANDPI}/allele_freq" \
         "${OUT_BANDPI}/pi_curves" \
         "${OUT_BANDPI}/mean_tP" \
         "${OUT_BANDPI}/plots" \
         "${OUT_BANDPI}/logs"

LOG_PI="${OUT_BANDPI}/logs"

hr_log "=== STEP_PI01: Karyotype-band nucleotide diversity (pi11/pi22/pi12) ==="
hr_log "Inversion candidates : ${INV_CANDIDATES}"
hr_log "Karyotype calls      : ${KARYOTYPE_CALLS}"
hr_log "pestPG dir (optional): ${PESTPG_DIR}"
hr_log "Output               : ${OUT_BANDPI}"
hr_log "Window / step        : ${PI_WIN} / ${PI_STEP}"
hr_log "Min N per band       : ${MIN_BAND_N}"

hr_check_file "${INV_CANDIDATES}"   "inversion candidates table"
hr_check_file "${KARYOTYPE_CALLS}"  "karyotype calls table"
hr_check_file "${REF}"              "reference"
hr_check_file "${CALLABLE_SITES}"   "ANGSD callable sites"
hr_check_file "${SAMPLE_MANIFEST}"  "sample manifest"

hr_check_cmd angsd
hr_check_cmd "${RSCRIPT_BIN}"

# Optional subset of inversion IDs from CLI
INV_FILTER=("$@")
inv_in_filter() {
  [[ ${#INV_FILTER[@]} -eq 0 ]] && return 0
  local needle="$1"
  for x in "${INV_FILTER[@]}"; do [[ "$x" == "$needle" ]] && return 0; done
  return 1
}

# ── Per-inversion loop ───────────────────────────────────────────────────────
# inversion_candidates.tsv: header inversion_id chr start end
{
  read -r _hdr  # discard header
  while IFS=$'\t' read -r INV CHR START END _rest; do
    [[ -z "${INV}" || "${INV}" =~ ^# ]] && continue
    inv_in_filter "${INV}" || continue

    hr_log "----- Inversion ${INV}: ${CHR}:${START}-${END} -----"

    INV_OUT="${OUT_BANDPI}"
    BL_HOMO1="${INV_OUT}/bamlists/${INV}.Homo_1.bamlist"
    BL_HOMO2="${INV_OUT}/bamlists/${INV}.Homo_2.bamlist"
    BL_HET="${INV_OUT}/bamlists/${INV}.Het.bamlist"
    BL_SHARED="${INV_OUT}/bamlists/${INV}.Homo_1_2.shared.bamlist"

    # ── 1. Build per-band bamlists from karyotype + manifest ────────────────
    hr_log "  [1/8] Build bamlists"
    "${RSCRIPT_BIN}" "${MODDIR}/scripts/make_karyotype_bamlists.R" \
      --karyotype  "${KARYOTYPE_CALLS}" \
      --manifest   "${SAMPLE_MANIFEST}" \
      --inversion  "${INV}" \
      --out-homo1  "${BL_HOMO1}" \
      --out-homo2  "${BL_HOMO2}" \
      --out-het    "${BL_HET}" \
      2> "${LOG_PI}/${INV}.bamlists.log"

    N1=$(wc -l < "${BL_HOMO1}" 2>/dev/null || echo 0)
    N2=$(wc -l < "${BL_HOMO2}" 2>/dev/null || echo 0)
    NH=$(wc -l < "${BL_HET}"   2>/dev/null || echo 0)
    hr_log "    N(Homo_1)=${N1}  N(Het)=${NH}  N(Homo_2)=${N2}"

    if (( N1 < MIN_BAND_N || N2 < MIN_BAND_N )); then
      hr_log "    SKIP ${INV}: too few samples in Homo_1 (${N1}) or Homo_2 (${N2}); min=${MIN_BAND_N}"
      continue
    fi

    # ── 2. Combined Homo_1 + Homo_2 bamlist (for shared allele coding only) ─
    cat "${BL_HOMO1}" "${BL_HOMO2}" > "${BL_SHARED}"

    # ── 3. ANGSD on shared list to define major/minor over the interval ─────
    hr_log "  [2/8] ANGSD on shared (Homo_1 + Homo_2) to infer major/minor"
    SHARED_PFX="${INV_OUT}/allele_freq/${INV}.shared"
    angsd \
      -b "${BL_SHARED}" \
      -r "${CHR}:${START}-${END}" \
      -ref "${REF}" \
      -anc "${REF}" \
      -GL 1 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -doCounts 1 \
      -SNP_pval 1e-6 \
      -minQ "${MINQ}" \
      -minMapQ "${MINMAPQ}" \
      -C "${CLIP}" \
      -sites "${CALLABLE_SITES}" \
      -nThreads "${THREADS}" \
      -out "${SHARED_PFX}" \
      2> "${LOG_PI}/${INV}.shared_angsd.log"

    if [[ ! -s "${SHARED_PFX}.mafs.gz" ]]; then
      hr_err "    No shared mafs produced for ${INV}; skipping"
      continue
    fi

    # ── 4. Build 4-column sites file: chr pos major minor ──────────────────
    hr_log "  [3/8] Build 4-column shared sites file"
    SITES4="${INV_OUT}/allele_freq/${INV}.shared.majorMinor.sites"
    # mafs.gz columns: chromo position major minor knownEM/unknownEM ...
    zcat "${SHARED_PFX}.mafs.gz" \
      | awk 'NR>1 && $3!="N" && $4!="N" {print $1"\t"$2"\t"$3"\t"$4}' \
      > "${SITES4}"
    N_SITES4=$(wc -l < "${SITES4}")
    hr_log "    ${N_SITES4} shared sites with fixed major/minor"

    if (( N_SITES4 < 1 )); then
      hr_err "    Empty shared sites file for ${INV}; skipping"
      continue
    fi

    # ── 5. Index the 4-column sites file ───────────────────────────────────
    angsd sites index "${SITES4}" 2>> "${LOG_PI}/${INV}.shared_angsd.log"

    # ── 6. Per-band ANGSD with -doMajorMinor 3 + shared sites ─────────────
    run_band_angsd() {
      local BAND="$1" BL="$2"
      local PFX="${INV_OUT}/allele_freq/${INV}.${BAND}"
      hr_log "  [4/8] ANGSD on ${BAND} (-doMajorMinor 3, shared sites)"
      angsd \
        -b "${BL}" \
        -r "${CHR}:${START}-${END}" \
        -ref "${REF}" \
        -anc "${REF}" \
        -GL 1 \
        -doMajorMinor 3 \
        -sites "${SITES4}" \
        -doMaf 1 \
        -doCounts 1 \
        -minQ "${MINQ}" \
        -minMapQ "${MINMAPQ}" \
        -C "${CLIP}" \
        -nThreads "${THREADS}" \
        -out "${PFX}" \
        2> "${LOG_PI}/${INV}.${BAND}_angsd.log"
      [[ -s "${PFX}.mafs.gz" ]] || hr_die "    ${BAND} mafs missing for ${INV}"
    }
    run_band_angsd "Homo_1" "${BL_HOMO1}"
    run_band_angsd "Homo_2" "${BL_HOMO2}"

    # ── 7. Compute pi11/pi22/pi12 in R ─────────────────────────────────────
    hr_log "  [5/8] Compute pi11 / pi22 / pi12"
    "${RSCRIPT_BIN}" "${MODDIR}/scripts/compute_pi11_pi22_pi12.R" \
      --homo1        "${INV_OUT}/allele_freq/${INV}.Homo_1.mafs.gz" \
      --homo2        "${INV_OUT}/allele_freq/${INV}.Homo_2.mafs.gz" \
      --inv          "${INV}" \
      --chr          "${CHR}" \
      --start        "${START}" \
      --end          "${END}" \
      --win          "${PI_WIN}" \
      --step         "${PI_STEP}" \
      --site-basis   "variant_only" \
      --out-sitewise "${INV_OUT}/pi_curves/${INV}.pi11_pi22_pi12.sitewise.tsv.gz" \
      --out-window   "${INV_OUT}/pi_curves/${INV}.pi11_pi22_pi12.window.tsv" \
      2> "${LOG_PI}/${INV}.pi.log"

    # ── 8. Optional: average per-sample pestPG by band ──────────────────────
    if [[ -d "${PESTPG_DIR}" ]]; then
      hr_log "  [6/8] Average per-sample pestPG by karyotype band"
      "${RSCRIPT_BIN}" "${MODDIR}/scripts/average_sample_pestPG_by_karyotype.R" \
        --karyotype  "${KARYOTYPE_CALLS}" \
        --pestpg-dir "${PESTPG_DIR}" \
        --inversion  "${INV}" \
        --chr        "${CHR}" \
        --start      "${START}" \
        --end        "${END}" \
        --out        "${INV_OUT}/mean_tP/${INV}.mean_tP_by_band.tsv" \
        2> "${LOG_PI}/${INV}.mean_tP.log" \
        || hr_log "    pestPG averaging failed for ${INV} (non-fatal)"
    else
      hr_log "  [6/8] pestPG dir not found, skipping mean_tP step"
    fi

    # ── 9. Two-panel plot ───────────────────────────────────────────────────
    hr_log "  [7/8] Plot band pi panel"
    "${RSCRIPT_BIN}" "${MODDIR}/scripts/plot_band_pi_panel.R" \
      --inv         "${INV}" \
      --chr         "${CHR}" \
      --start       "${START}" \
      --end         "${END}" \
      --window-tsv  "${INV_OUT}/pi_curves/${INV}.pi11_pi22_pi12.window.tsv" \
      --mean-tp-tsv "${INV_OUT}/mean_tP/${INV}.mean_tP_by_band.tsv" \
      --out-pdf     "${INV_OUT}/plots/${INV}.band_pi_panel.pdf" \
      2> "${LOG_PI}/${INV}.plot.log" \
      || hr_log "    plot failed for ${INV} (non-fatal)"

    hr_log "  [8/8] ${INV} done."
  done
} < "${INV_CANDIDATES}"

hr_log "=== STEP_PI01 complete ==="
hr_log "Outputs: ${OUT_BANDPI}"
