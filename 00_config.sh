#!/usr/bin/env bash
# =============================================================================
# 00_config.sh — catfish-diversity-analysis
# =============================================================================
# Single source of truth for paths, parameters, and SLURM defaults used by
# every module in this repo. Source from any script:
#
#   source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
#       — or, from a Module subdir —
#   source "$(dirname "${BASH_SOURCE[0]}")/../../00_config.sh"
#
# This repo is for the 226-sample pure *C. gariepinus* hatchery cohort.
# Do not use for other cohorts without auditing every parameter.
# =============================================================================

set -euo pipefail

# ── Project root ─────────────────────────────────────────────────────────────
export BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── This repo's location ────────────────────────────────────────────────────
export REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MODULES_DIR="${REPO_DIR}/Modules"

# ── Sibling repos ───────────────────────────────────────────────────────────
export VARIANT_REPO="${VARIANT_REPO:-${BASE}/catfish-variant-analysis}"
export POPULATION_REPO="${POPULATION_REPO:-${BASE}/catfish-population-analysis}"

# ── Reference / cohort-stable inputs ────────────────────────────────────────
export REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
export REF_FAI="${REF}.fai"
export ANC="${REF}"                     # ref-as-ancestor → folded SFS
export SAMPLES_IND="${BASE}/het_roh/01_inputs_check/samples.ind"
export BAMLIST="${BASE}/het_roh/01_inputs_check/bamlist_qcpass.txt"
export SAMPLE_LIST="${BASE}/01_inputs_check/samples_226_pure_gariepinus.txt"
export CALLABLE_SITES="${BASE}/01_inputs_check/callable_sites.bed.gz"

# Linkage groups
export CHROM_LIST=(C_gar_LG01 C_gar_LG02 C_gar_LG03 C_gar_LG04 C_gar_LG05
                   C_gar_LG06 C_gar_LG07 C_gar_LG08 C_gar_LG09 C_gar_LG10
                   C_gar_LG11 C_gar_LG12 C_gar_LG13 C_gar_LG14 C_gar_LG15
                   C_gar_LG16 C_gar_LG17 C_gar_LG18 C_gar_LG19 C_gar_LG20
                   C_gar_LG21 C_gar_LG22 C_gar_LG23 C_gar_LG24 C_gar_LG25
                   C_gar_LG26 C_gar_LG27 C_gar_LG28)
export N_CHROM=${#CHROM_LIST[@]}

# ── Outputs ─────────────────────────────────────────────────────────────────
export OUTROOT="${OUTROOT:-${BASE}/results/catfish-diversity-analysis}"
export OUT_SAF="${OUTROOT}/01_saf_per_sample"
export OUT_HETEROZYGOSITY="${OUTROOT}/02_heterozygosity"
export OUT_THETA_PI="${OUTROOT}/03_theta_pi"
export OUT_ROH="${OUTROOT}/04_roh"
export OUT_AGGREGATED="${OUTROOT}/05_aggregated"
export LOG_DIR="${OUTROOT}/logs"

# ── θπ window-grid parameters ───────────────────────────────────────────────
# Default scale matches the inversion-atlas page-12 enrichment: 10 kb windows
# stepped 2 kb. The pestPG `tP` column is a per-window SUM, not a per-site
# density — modules in 03_theta_pi/ divide by nSites to get the diversity-
# comparable per-site estimate. See docs/methods/theta_pi_scaling.md.
export PESTPG_SCALE="${PESTPG_SCALE:-win10000.step2000}"
export THETA_WIN_BP=10000
export THETA_STEP_BP=2000
# Multi-scale set for heterozygosity (when sweeping window scales):
export THETA_SCALE_SET=("win10000.step2000" "win50000.step10000")

# ── ROH parameters ──────────────────────────────────────────────────────────
export ROH_MIN_LEN_BP=100000
export ROH_MIN_SNPS=20
export ROH_MAX_GAP_BP=50000

# ── Compute / SLURM ─────────────────────────────────────────────────────────
export RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"
export MAMBA_ENV="${MAMBA_ENV:-assembly}"
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-lt200308}"
export SLURM_PARTITION="${SLURM_PARTITION:-compute}"

# ── Helpers ─────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%F %T')] [div] $*"; }
die() { echo "[$(date '+%F %T')] [div] [ERROR] $*" >&2; exit 1; }

init_dirs() {
  mkdir -p "${OUTROOT}" \
           "${OUT_SAF}" "${OUT_HETEROZYGOSITY}" "${OUT_THETA_PI}" \
           "${OUT_ROH}" "${OUT_AGGREGATED}" \
           "${LOG_DIR}"
}

config_print() {
  cat <<EOF
=== catfish-diversity-analysis ===
BASE             : ${BASE}
REPO_DIR         : ${REPO_DIR}
OUTROOT          : ${OUTROOT}
N_CHROM          : ${N_CHROM}
PESTPG_SCALE     : ${PESTPG_SCALE}
RSCRIPT_BIN      : ${RSCRIPT_BIN}
==================================
EOF
}

if [[ "${1:-}" == "print" ]]; then config_print; fi
