#!/usr/bin/env bash
# =============================================================================
# STEP_C01_emit_atlas_payloads.sh
# =============================================================================
# Runs the four diversity-atlas adapters and writes their JSONs into
# ${OUTBASE}/13_atlas_payloads/. Each adapter degrades to a schema-valid
# v0-stub when its upstream inputs aren't yet present, so this driver is
# safe to run before all sibling pipelines exist.
#
# Optional: --copy-to-atlas <DIVERSITY_ATLAS_REPO> also copies the four
# emitted JSONs into <repo>/data/ so the atlas picks them up immediately
# (avoids a manual cp step during dev).
#
# Usage:
#   bash STEP_C01_emit_atlas_payloads.sh
#   bash STEP_C01_emit_atlas_payloads.sh --copy-to-atlas /path/to/diversity-atlas
# =============================================================================

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO_DIR}/00_config.sh"

ADAPTERS_DIR="${REPO_DIR}/06_atlas_adapters"
OUT_DIR="${OUTBASE}/13_atlas_payloads"
mkdir -p "${OUT_DIR}"

COPY_TARGET=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --copy-to-atlas) COPY_TARGET="$2"; shift 2 ;;
    -h|--help) sed -n '1,30p' "$0"; exit 0 ;;
    *) hr_die "unknown arg: $1" ;;
  esac
done

PY="${PY:-python3}"

hr_log "OUT_DIR = ${OUT_DIR}"

# ---------------------------------------------------------------------------
# 1) texture_metrics — module 04 outputs
# ---------------------------------------------------------------------------
hr_log "[1/4] texture_metrics"
TEXTURE_DIR="${OUTBASE}/11_window_diversity_texture/02_window_metrics"
"${PY}" "${ADAPTERS_DIR}/adapt_texture_metrics.py" \
  --per-sample      "${TEXTURE_DIR}/per_sample_summary.tsv" \
  --cohort-median   "${TEXTURE_DIR}/cohort_median_per_window.tsv" \
  --long-matrix     "${TEXTURE_DIR}/long_window_matrix.tsv.gz" \
  --het-summary     "${OUT_HETEROZYGOSITY}/04_summary/genomewide_heterozygosity.tsv" \
  --win-bp          50000 \
  --step-bp         50000 \
  --min-callable    15000 \
  --out             "${OUT_DIR}/texture_metrics.json"

# ---------------------------------------------------------------------------
# 2) roh_gene_overlap — module 05 outputs
# ---------------------------------------------------------------------------
hr_log "[2/4] roh_gene_overlap"
ROH_GENES_DIR="${OUTBASE}/12_roh_metrics_and_genes/genes_in_roh"
ROH_BED="${OUT_ROH}/catfish_roh.bed"
"${PY}" "${ADAPTERS_DIR}/adapt_roh_gene_overlap.py" \
  --roh-bed             "${ROH_BED}" \
  --genes-bed           "${ROH_GENES_DIR}/genes.bed" \
  --cohort-recurrence   "${ROH_GENES_DIR}/cohort_recurrence.tsv" \
  --k8-assign           "${OUTBASE}/14_atlas_aux/k8_assign.tsv" \
  --family-assign       "${OUTBASE}/14_atlas_aux/family_assign.tsv" \
  --froh-summary        "${OUT_ROH}/per_sample_froh.tsv" \
  --constraint-proxy    "${OUTBASE}/14_atlas_aux/gene_constraint.tsv" \
  --n-top-peaks         24 \
  --out                 "${OUT_DIR}/roh_gene_overlap.json"

# ---------------------------------------------------------------------------
# 3) divergence_network — upstream pairwise FST pipeline (not yet shipped)
# ---------------------------------------------------------------------------
hr_log "[3/4] divergence_network"
DIVERGE_DIR="${OUTBASE}/15_divergence_network"
"${PY}" "${ADAPTERS_DIR}/adapt_divergence_network.py" \
  --groups-tsv      "${DIVERGE_DIR}/groups_K8.tsv" \
  --edges-tsv       "${DIVERGE_DIR}/edges_K8.tsv" \
  --grouping        "K=8" \
  --fst-estimator   "Weir-Cockerham" \
  --n-bootstrap     100 \
  --min-callable    100000 \
  --out             "${OUT_DIR}/divergence_network.json"

# ---------------------------------------------------------------------------
# 4) functional_burden — upstream snpEff/VESM pipeline (not yet shipped)
# ---------------------------------------------------------------------------
hr_log "[4/4] functional_burden"
BURDEN_DIR="${OUTBASE}/16_functional_burden"
"${PY}" "${ADAPTERS_DIR}/adapt_functional_burden.py" \
  --per-sample          "${BURDEN_DIR}/per_sample_burden.tsv" \
  --per-group           "${BURDEN_DIR}/per_group_K8.tsv" \
  --per-group-key       "family=${BURDEN_DIR}/per_group_family.tsv" \
  --per-group-key       "F_ROH_quartile=${BURDEN_DIR}/per_group_froh_quartile.tsv" \
  --variant-inventory   "${BURDEN_DIR}/variant_inventory.tsv" \
  --snpeff-totals       "${BURDEN_DIR}/snpeff_totals.tsv" \
  --gerp-inventory      "${BURDEN_DIR}/gerp_inventory.tsv" \
  --transcripts-json    "${BURDEN_DIR}/transcripts.json" \
  --msa-links-json      "${BURDEN_DIR}/msa_links.json" \
  --splice-events       "${BURDEN_DIR}/splice_events.tsv" \
  --pairwise-ks-json    "${BURDEN_DIR}/pairwise_ks.json" \
  --pin-pis-method      "NG86" \
  --pi0-pi4-method      "codon-degeneracy" \
  --block-size-bp       5000000 \
  --out                 "${OUT_DIR}/functional_burden.json"

# ---------------------------------------------------------------------------
# Optional publish
# ---------------------------------------------------------------------------
if [[ -n "${COPY_TARGET}" ]]; then
  hr_log "publishing to ${COPY_TARGET}/data/"
  if [[ ! -d "${COPY_TARGET}/data" ]]; then
    hr_die "target ${COPY_TARGET}/data/ does not exist — pass the path to the diversity-atlas repo root"
  fi
  for fname in texture_metrics.json roh_gene_overlap.json \
               divergence_network.json functional_burden.json; do
    cp -v "${OUT_DIR}/${fname}" "${COPY_TARGET}/data/${fname}"
  done
fi

hr_log "STEP_C01 done."
