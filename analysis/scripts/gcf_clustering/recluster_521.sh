#!/usr/bin/env bash
# recluster_521.sh
#
# Re-runs the GCF clustering pipeline on the deduplicated 521-genome cohort. This
# is a full re-cluster, not a post-hoc filter of the 554-genome assignments:
# duplicate BGCs can bias DBSCAN's density-based clustering itself, not just the
# downstream per-genome statistics.
#
# Requires the funcscan outputs (analysis/funcscan/results), HMMER on PATH, and
# the pipeline's Python dependencies (requirements.txt). hmmfetch and hmmalign
# need a pressed Pfam-A; point PFAM_HMM at it.
#
# Usage (from scripts/gcf_clustering/):
#   PFAM_HMM=/path/to/Pfam-A.hmm bash recluster_521.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RESULTS_DIR="${ANALYSIS_DIR}/scripts/results_521"

mkdir -p "${RESULTS_DIR}"

# --- filter combgc_complete_summary.tsv to the 521-genome cohort ---
COMBGC_FULL="${ANALYSIS_DIR}/R/copied_from_funcscan_results/combgc_complete_summary.tsv"
COMBGC_521="${RESULTS_DIR}/combgc_complete_summary_521.tsv"
awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1 || ($1 in keep)' \
    "${ANALYSIS_DIR}/cohort/ids_521.txt" "${COMBGC_FULL}" > "${COMBGC_521}"
test -s "${COMBGC_521}"

export COMBGC_INPUT="${COMBGC_521}"
export RESULTS_DIR
export FASTA_DIR="${RESULTS_DIR}/fastas_proteins"

cd "${SCRIPT_DIR}"
bash run_clustering.sh
test -s "${RESULTS_DIR}/gcf_assignments.tsv"
