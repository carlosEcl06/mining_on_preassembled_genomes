#!/usr/bin/env bash
# run_clustering.sh
#
# Orchestrates the GCF clustering pipeline:
#   00_extract_bgc_fastas.py  →  per-BGC protein FASTAs (from antiSMASH GBKs)
#   01_extract_domains.py     →  bgc_domain_arrays.tsv
#   02_jaccard_matrix.py      →  jaccard_pairs.tsv
#   03_backbone_identity.py   →  backbone_identity.tsv
#   04_dbscan.py              →  gcf_assignments.tsv
#
# Usage (from scripts/gcf_clustering/):
#   bash run_clustering.sh
#
# Optional overrides (environment variables):
#   COMBGC_INPUT       path to combgc_complete_summary.tsv
#   RESULTS_DIR        path to output directory
#   ANTISMASH_DIR      path to funcscan antiSMASH results (per-sample GBKs)
#   FASTA_DIR          path to write per-BGC protein FASTAs
#   PFAM_HMM           path to Pfam-A.hmm (pressed with hmmpress)
#   PYTHON             python executable (default: python3)

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

COMBGC_INPUT="${COMBGC_INPUT:-${ANALYSIS_DIR}/R/copied_from_funcscan_results/combgc_complete_summary.tsv}"
RESULTS_DIR="${RESULTS_DIR:-${SCRIPT_DIR}/../results}"
PYTHON="${PYTHON:-python3}"

ANTISMASH_DIR="${ANTISMASH_DIR:-${ANALYSIS_DIR}/funcscan/results/bgc/antismash}"
FASTA_DIR="${FASTA_DIR:-${ANALYSIS_DIR}/rawdata/fasta/proteins}"
PFAM_HMM="${PFAM_HMM:-${ANALYSIS_DIR}/rawdata/pfam/Pfam-A.hmm}"
BACKBONE_IDENTITY="${RESULTS_DIR}/backbone_identity.tsv"

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
log "ANTISMASH_DIR: ${ANTISMASH_DIR}"
log "FASTA_DIR    : ${FASTA_DIR}"
log "PFAM_HMM     : ${PFAM_HMM}"
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
# Step 00 — extract per-BGC protein FASTAs from antiSMASH GBKs
# ---------------------------------------------------------------------------
log "Step 00: extracting per-BGC protein FASTAs from antiSMASH GBKs ..."
START=$(date +%s)

${PYTHON} "${SCRIPT_DIR}/00_extract_bgc_fastas.py" \
    --domains       "${DOMAIN_ARRAYS}" \
    --antismash-dir "${ANTISMASH_DIR}" \
    --output-dir    "${FASTA_DIR}" \
    2>&1 | tee -a "${LOG_FILE}"

log "Step 00 done in $(( $(date +%s) - START ))s"
hr

# ---------------------------------------------------------------------------
# Step 03 — backbone sequence identity (hmmalign)
# ---------------------------------------------------------------------------
log "Step 03: estimating backbone sequence identity (hmmalign) ..."
START=$(date +%s)

${PYTHON} "${SCRIPT_DIR}/03_backbone_identity.py" \
    --domains   "${DOMAIN_ARRAYS}" \
    --pairs     "${JACCARD_PAIRS}" \
    --fasta-dir "${FASTA_DIR}" \
    --pfam-hmm  "${PFAM_HMM}" \
    --output    "${BACKBONE_IDENTITY}" \
    2>&1 | tee -a "${LOG_FILE}"

log "Step 03 done in $(( $(date +%s) - START ))s"
hr

# ---------------------------------------------------------------------------
# Step 04 — DBSCAN clustering
# ---------------------------------------------------------------------------
log "Step 04: running DBSCAN clustering ..."
START=$(date +%s)

${PYTHON} "${SCRIPT_DIR}/04_dbscan.py" \
    --domains  "${DOMAIN_ARRAYS}" \
    --pairs    "${JACCARD_PAIRS}" \
    --identity "${BACKBONE_IDENTITY}" \
    --output   "${GCF_ASSIGNMENTS}" \
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
log "  backbone_identity.tsv  — pairwise backbone sequence identity"
log "  gcf_assignments.tsv    — GCF assignments (main output)"
log "  run_clustering.log     — this log"
hr