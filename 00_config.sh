#!/usr/bin/env bash
# =============================================================================
# 00_config.sh — catfish-diversity-analysis
# =============================================================================
# Single source of truth for paths, parameters, and SLURM defaults.
# Source from any script at the root:
#
#   source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
#
# This repo is for the 226-sample pure *C. gariepinus* hatchery cohort.
# Three result buckets: heterozygosity, theta_pi, ROH.
# =============================================================================

set -euo pipefail

# ── Project root ─────────────────────────────────────────────────────────────
export BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
export PROJECT="${BASE}"   # alias for legacy MODULE_3 references

# ── This repo's location ────────────────────────────────────────────────────
export REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Sibling repos ───────────────────────────────────────────────────────────
export VARIANT_REPO="${VARIANT_REPO:-${BASE}/catfish-variant-analysis}"
export POPULATION_REPO="${POPULATION_REPO:-${BASE}/catfish-population-analysis}"

# ── Reference / cohort-stable inputs ────────────────────────────────────────
export REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
export REF_FAI="${REF}.fai"
export ANC="${REF}"

# Callable mask: ANGSD format (1-based) and BED (0-based)
export CALLABLE_SITES="${BASE}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.1-based.angsd"
export CALLABLE_BED="${BASE}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.0based.bed"
export REPEAT_BED="${BASE}/fClaHyb_Gar_LG.mask_regions.softacgt.bed"

# Sample manifest + lists
export SAMPLE_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
export BAMLIST_QCPASS="${BASE}/het_roh/01_inputs_check/bamlist_qcpass.txt"
export SAMPLE_LIST="${BASE}/het_roh/01_inputs_check/samples_qcpass.txt"

# ── ngsF-HMM inputs (BEAGLE GLs from population-analysis MODULE_2A) ────────
export BEAGLE_GL="${BASE}/het_roh/01_inputs_check/main_qcpass.beagle.gz"
export POS_FILE="${BASE}/het_roh/01_inputs_check/main_qcpass.pos"
export SAMPLES_IND="${BASE}/het_roh/01_inputs_check/samples.ind"

# ── ngsF-HMM executable ────────────────────────────────────────────────────
export NGSFHMM_BIN="/project/lt200308-agbsci/01-catfish_assembly/ngsF-HMM/ngsF-HMM"
export NGSFHMM_WRAPPER="/project/lt200308-agbsci/01-catfish_assembly/ngsF-HMM/ngsF-HMM.sh"
export NGSFHMM_MODULES=(
  "HTSlib/1.17-cpeGNU-23.03"
  "Boost/1.81.0-cpeGNU-23.03"
  "GSL/2.7-cpeGNU-23.03"
)

# ── Ancestry / structure metadata (from population-analysis MODULE_2B) ─────
export NGSRELATE_PAIRS="${BASE}/popstruct_thin/05_ngsrelate/catfish_226_for_natora.txt"
export PRUNED81_SAMPLES="${BASE}/popstruct_thin/05_ngsrelate/catfish_first_degree_pairwise_toKeep.txt"
export BESTK_TABLE="${BASE}/popstruct_thin/05_ngsadmix_global/best_seed_by_K.tsv"
export NGSADMIX_DIR="${BASE}/popstruct_thin/05_ngsadmix_global/runs_thin500"
export EVALADMIX_DIR="${BASE}/popstruct_thin/05_ngsadmix_global/evaladmix_thin500"
export Q_MATRIX="${BASE}/popstruct_thin/05_ngsadmix_global/runs_thin500/thin500_K08_best.qopt"
export COVTREE_ORDER="${BASE}/popstruct_thin/05_pcaangsd_global/thin_500/K4/catfish.wg.byRF.thin_500.pcangsd.tree"
export ANCESTRY_LABELS="${BASE}/popstruct_thin/05_ngsadmix_global/sample_main_ancestry_by_K.tsv"

# ── Linkage groups ─────────────────────────────────────────────────────────
export CHROM_LIST=(C_gar_LG01 C_gar_LG02 C_gar_LG03 C_gar_LG04 C_gar_LG05
                   C_gar_LG06 C_gar_LG07 C_gar_LG08 C_gar_LG09 C_gar_LG10
                   C_gar_LG11 C_gar_LG12 C_gar_LG13 C_gar_LG14 C_gar_LG15
                   C_gar_LG16 C_gar_LG17 C_gar_LG18 C_gar_LG19 C_gar_LG20
                   C_gar_LG21 C_gar_LG22 C_gar_LG23 C_gar_LG24 C_gar_LG25
                   C_gar_LG26 C_gar_LG27 C_gar_LG28)
export N_CHROM=${#CHROM_LIST[@]}

# ── Outputs: three result buckets ─────────────────────────────────────────
export OUTBASE="${OUTBASE:-${BASE}/het_roh}"
export OUTROOT="${OUTROOT:-${OUTBASE}}"

export DIR_INPUTS="${OUTBASE}/01_inputs_check"
export OUT_HETEROZYGOSITY="${OUTBASE}/02_heterozygosity"
export DIR_NGSFHMM="${OUTBASE}/03_ngsF_HMM"
export OUT_ROH="${OUTBASE}/04_roh_summary"
export DIR_INV="${OUTBASE}/05_inversion_support"
export DIR_PLOTS_CORE="${OUTBASE}/06_plots_core"
export DIR_PLOTS_META="${OUTBASE}/07_plots_metadata"
export DIR_STATS="${OUTBASE}/08_stats"
export OUT_TABLES="${OUTBASE}/09_final_tables"
export DIR_REPORT="${OUTBASE}/10_report"
export DIR_LOGS="${OUTBASE}/logs"

# Aliases used by legacy MODULE_3 scripts (DIR_HET, DIR_ROH, etc.)
export DIR_HET="${OUT_HETEROZYGOSITY}"
export DIR_ROH="${OUT_ROH}"
export DIR_TABLES="${OUT_TABLES}"
export LOG_DIR="${DIR_LOGS}"
export OUT_THETA_PI="${OUT_HETEROZYGOSITY}/03_theta"

# ── ANGSD QC filters ───────────────────────────────────────────────────────
export MINQ=20
export MINMAPQ=30
export CLIP=50

# ── realSFS convergence ───────────────────────────────────────────────────
export REALSFS_MAXITER=2000
export REALSFS_TOLE="1e-16"

# ── thetaStat windows ─────────────────────────────────────────────────────
# Main scale: 500 kb non-overlapping
export WIN=500000
export STEP=500000
# Multiscale (local diversity landscape)
export RUN_EXTRA_THETA_SCALES=1
export THETA_SCALES=("5000_1000" "10000_2000" "50000_10000")
export THETA_SCALE_LABELS=("5kb_1kb" "10kb_2kb" "50kb_10kb")

# ── ROH parameters ────────────────────────────────────────────────────────
export NGSFHMM_REPS=10
export NGSFHMM_SEED_BASE=42
export ROH_MIN_LEN_BP=100000
export ROH_MIN_SNPS=20
export ROH_MAX_GAP_BP=50000

# ── Compute / SLURM ───────────────────────────────────────────────────────
export THREADS="${SLURM_CPUS_PER_TASK:-8}"
export RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"
export MAMBA_ENV="${MAMBA_ENV:-assembly}"
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-lt200308}"
export SLURM_PARTITION="${SLURM_PARTITION:-compute}"

# ── Helpers ────────────────────────────────────────────────────────────────
hr_timestamp() { date '+%F %T'; }
hr_log()  { echo "[$(hr_timestamp)] [div] $*"; }
hr_err()  { echo "[$(hr_timestamp)] [div] [ERROR] $*" >&2; }
hr_die()  { hr_err "$@"; exit 1; }
hr_check_file() {
  local f="$1" label="${2:-file}"
  [[ -e "$f" ]] || hr_die "Missing ${label}: ${f}"
}
hr_check_cmd() { command -v "$1" &>/dev/null || hr_die "Command not found: $1"; }
hr_init_dirs() {
  mkdir -p \
    "$DIR_INPUTS" "$OUT_HETEROZYGOSITY" "$DIR_NGSFHMM" "$OUT_ROH" \
    "$DIR_INV" "$DIR_PLOTS_CORE" "$DIR_PLOTS_META" \
    "$DIR_STATS" "$OUT_TABLES" "$DIR_REPORT" "$DIR_LOGS" \
    "${OUT_HETEROZYGOSITY}/01_saf" "${OUT_HETEROZYGOSITY}/02_sfs" \
    "${OUT_HETEROZYGOSITY}/03_theta" "${OUT_HETEROZYGOSITY}/03_theta/multiscale" \
    "${OUT_HETEROZYGOSITY}/04_summary"
}

# Short-name aliases for any non-MODULE_3 scripts
log()       { hr_log "$@"; }
die()       { hr_die "$@"; }
init_dirs() { hr_init_dirs; }

config_print() {
  cat <<EOF
=== catfish-diversity-analysis ===
BASE             : ${BASE}
REPO_DIR         : ${REPO_DIR}
OUTBASE          : ${OUTBASE}
N_CHROM          : ${N_CHROM}
WIN/STEP         : ${WIN}/${STEP}
THETA_SCALES     : ${THETA_SCALES[*]}
NGSFHMM_REPS     : ${NGSFHMM_REPS}
==================================
EOF
}

if [[ "${1:-}" == "print" ]]; then config_print; fi
