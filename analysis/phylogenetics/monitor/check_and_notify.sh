#!/usr/bin/env bash
# check_and_notify.sh — hands-free PBS job monitor for the lab server.
#
# Invoked periodically by Windows Task Scheduler (via wsl.exe), independent of
# whether an interactive Claude Code session is open — that's what makes this
# "hands-free": it fires on a timer as long as the laptop is on, no terminal
# needed. See ../../../.claude adjacent notes / the project plan for the full
# design rationale (cloud-routine SSH was ruled out — see plan history).
#
# Design: plain bash does the actual SSH/qstat work (deterministic, no LLM
# judgment needed, so no permission risk in an unattended run). A tiny,
# narrowly-scoped `claude -p` call — allowed ONLY the PushNotification tool,
# nothing else — fires exactly when a tracked job's state changed, and stays
# silent otherwise. Routine checks with nothing to report cost zero LLM usage.
#
# Setup:
#   - List jobs to watch in tracked_jobs.tsv (job_id<TAB>label), one per line,
#     git-tracked — edit it as new PBS jobs are submitted during the pipeline.
#   - .state.json (gitignored) holds last-seen status per job_id.
#   - .monitor.log (gitignored) is a running log for debugging.
#
# Scheduled via Windows Task Scheduler task "LabServerPBSMonitor" (created
# 2026-09-16), firing daily at 08:00/10:00/12:00/14:00/16:00/18:00.
set -uo pipefail   # NOT -e: a single job's ssh/qstat hiccup must not abort the whole check

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRACKED_JOBS="${SCRIPT_DIR}/tracked_jobs.tsv"
STATE_FILE="${SCRIPT_DIR}/.state.json"
LOG_FILE="${SCRIPT_DIR}/.monitor.log"

SSH_KEY="${HOME}/.ssh/id_ed25519_labserver"
SSH_HOST="clima@143.107.143.27"
SSH_PORT="2222"
SSH_OPTS=(-i "${SSH_KEY}" -p "${SSH_PORT}" -o BatchMode=yes -o ConnectTimeout=15 -o StrictHostKeyChecking=accept-new)

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "${LOG_FILE}"; }

log "--- check started ---"

if [[ ! -s "${TRACKED_JOBS}" ]]; then
    log "No tracked jobs yet (${TRACKED_JOBS} empty/missing); nothing to check."
    exit 0
fi

[[ -f "${STATE_FILE}" ]] || echo '{}' > "${STATE_FILE}"

CHANGES=""

while IFS=$'\t' read -r job_id label; do
    [[ -z "${job_id}" || "${job_id}" == \#* ]] && continue

    current=$(ssh "${SSH_OPTS[@]}" "${SSH_HOST}" "qstat -f ${job_id} 2>&1" 2>/dev/null \
        | grep -E "job_state|Exit_status" | tr -d ' \t' | tr '\n' ';')
    if [[ -z "${current}" ]]; then
        current="GONE"   # left the queue entirely: finished+purged, or never existed
    fi

    prev=$(python3 -c "
import json
try:
    d = json.load(open('${STATE_FILE}'))
except Exception:
    d = {}
print(d.get('${job_id}', ''))
" 2>/dev/null)

    if [[ "${current}" != "${prev}" ]]; then
        CHANGES="${CHANGES}${label} (${job_id}): ${prev:-<none>} -> ${current}\n"
        log "CHANGE ${job_id} (${label}): '${prev:-<none>}' -> '${current}'"
    fi

    python3 -c "
import json
try:
    d = json.load(open('${STATE_FILE}'))
except Exception:
    d = {}
d['${job_id}'] = '''${current}'''
json.dump(d, open('${STATE_FILE}', 'w'))
"
done < "${TRACKED_JOBS}"

if [[ -n "${CHANGES}" ]]; then
    log "Changes detected, invoking claude -p for notification:"
    log "${CHANGES}"
    claude -p "One or more tracked PBS jobs on the lab server changed state since the last check. Send exactly one PushNotification (status: proactive) summarizing this — under 200 characters, one line, no markdown, lead with what changed and what (if anything) needs a decision. Changes:
${CHANGES}" \
        --model claude-haiku-4-5-20251001 \
        --allowedTools PushNotification \
        >> "${LOG_FILE}" 2>&1
else
    log "No changes."
fi

log "--- check finished ---"
