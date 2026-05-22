#!/usr/bin/env bash
# run_clustering.sh
#
# Orchestrates the GCF clustering pipeline:
#   01_extract_domains.py  →  bgc_domain_arrays.tsv
#   02_jaccard_matrix.py   →  jaccard_pairs.tsv
#   04_dbscan.py           →  gcf_assignments.tsv
#
# 03_backbone_identity.py is halted pending access to antiSMASH GenBank files.
#
# Usage (from scripts/gcf_clustering/):
#   bash run_clustering.sh
#
# Optional overrides (environment variables):
#   COMBGC_INPUT   path to combgc_complete_summary.tsv
#   RESULTS_DIR    path to output directory
#   PYTHON         python executable (default: python3)

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

COMBGC_INPUT="${COMBGC_INPUT:-${ANALYSIS_DIR}/R/copied_from_funcscan_results/combgc_complete_summary.tsv}"
RESULTS_DIR="${RESULTS_DIR:-${SCRIPT_DIR}/../results}"
PYTHON="${PYTHON:-python3}"

DOMAIN_ARRAYS="${RESULTS_DIR}/bgc_domain_arrays.tsv"
JACCARD_PAIRS="${RESULTS_DIR}/jaccard_pairs.tsv"
GCF_ASSIGNMENTS="${RESULTS_DIR}/gcf_assignments.tsv"

LOG_FILE="${RESULTS_DIR}/run_clustering.log"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"; }
hr()  { echo "─────────────────────────────────────────" | tee -a "${LOG_FILE}"; }

mkdir -p "${RESULTS_DIR}"
: > "${LOG_FILE}"   # truncate log

hr
log "GCF clustering pipeline starting"
log "COMBGC_INPUT : ${COMBGC_INPUT}"
log "RESULTS_DIR  : ${RESULTS_DIR}"
log "PYTHON       : $(${PYTHON} --version 2>&1)"
hr

# ---------------------------------------------------------------------------
# Step 01 — extract domains
# ---------------------------------------------------------------------------
log "Step 01: extracting Pfam domain arrays ..."
START=$(date +%s)

${PYTHON} "${SCRIPT_DIR}/01_extract_domains.py" \
    --input  "${COMBGC_INPUT}" \
    --output "${DOMAIN_ARRAYS}" \
    2>&1 | tee -a "${LOG_FILE}"

log "Step 01 done in $(( $(date +%s) - START ))s"
hr

# ---------------------------------------------------------------------------
# Step 02 — Jaccard matrix
# ---------------------------------------------------------------------------
log "Step 02: computing pairwise Jaccard similarities ..."
START=$(date +%s)

${PYTHON} "${SCRIPT_DIR}/02_jaccard_matrix.py" \
    --input  "${DOMAIN_ARRAYS}" \
    --output "${JACCARD_PAIRS}" \
    2>&1 | tee -a "${LOG_FILE}"

log "Step 02 done in $(( $(date +%s) - START ))s"
hr

# ---------------------------------------------------------------------------
# Step 03 — halted (backbone sequence identity via hmmalign)
# ---------------------------------------------------------------------------
log "Step 03: HALTED — backbone_identity skipped (antiSMASH GBK files not local)."
log "         Similarity score will use Jaccard only: sqrt(0.2 * Jaccard)."
hr

# ---------------------------------------------------------------------------
# Step 04 — DBSCAN clustering
# ---------------------------------------------------------------------------
log "Step 04: running DBSCAN clustering ..."
START=$(date +%s)

${PYTHON} "${SCRIPT_DIR}/04_dbscan.py" \
    --domains "${DOMAIN_ARRAYS}" \
    --pairs   "${JACCARD_PAIRS}" \
    --output  "${GCF_ASSIGNMENTS}" \
    2>&1 | tee -a "${LOG_FILE}"

log "Step 04 done in $(( $(date +%s) - START ))s"
hr

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
log "Pipeline complete."
log "Results in: ${RESULTS_DIR}"
log "  bgc_domain_arrays.tsv  — filtered BGCs with domain arrays"
log "  jaccard_pairs.tsv      — sparse pairwise Jaccard similarities"
log "  gcf_assignments.tsv    — GCF assignments (main output)"
log "  run_clustering.log     — this log"
hr