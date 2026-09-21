#!/usr/bin/env bash
# check_and_notify.sh — hands-free PBS job monitor for the lab server.
#
# Invoked periodically by Windows Task Scheduler (via wsl.exe), independent of
# whether an interactive Claude Code session is open — that's what makes this
# "hands-free": it fires on a timer as long as the laptop is on, no terminal
# needed. See ../../../.claude adjacent notes / the project plan for the full
# design rationale (cloud-routine SSH was ruled out — see plan history).
#
# Design: plain bash does the SSH/qstat work and, when a tracked job or file
# changed state, posts a plain-text message to a private ntfy.sh topic (phone
# push). Silent otherwise; no LLM involved.
#
# Setup:
#   - List jobs to watch in tracked_jobs.tsv (job_id<TAB>label), one per line,
#     git-tracked — edit it as new PBS jobs are submitted during the pipeline.
#   - List files to watch for existence in watched_files.tsv (path<TAB>label) —
#     for Nextflow-orchestrated runs, where individual PBS job IDs are
#     dispatched dynamically and can't be enumerated in advance; watching the
#     expected final output file appear is a cleaner completion signal.
#   - .state.json (gitignored) holds last-seen status per job_id/file.
#   - .monitor.log (gitignored) is a running log for debugging.
#
# Scheduled via Windows Task Scheduler task "LabServerPBSMonitor" (created
# 2026-09-16), firing daily at 08:00/10:00/12:00/14:00/16:00/18:00.
set -uo pipefail   # NOT -e: a single job's ssh/qstat hiccup must not abort the whole check

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRACKED_JOBS="${SCRIPT_DIR}/tracked_jobs.tsv"
WATCHED_FILES="${SCRIPT_DIR}/watched_files.tsv"
STATE_FILE="${SCRIPT_DIR}/.state.json"
LOG_FILE="${SCRIPT_DIR}/.monitor.log"

SSH_KEY="${HOME}/.ssh/id_ed25519_labserver"
SSH_HOST="clima@143.107.143.27"
SSH_PORT="2222"
SSH_OPTS=(-i "${SSH_KEY}" -p "${SSH_PORT}" -o BatchMode=yes -o ConnectTimeout=15 -o StrictHostKeyChecking=accept-new)

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "${LOG_FILE}"; }

log "--- check started ---"

if [[ ! -s "${TRACKED_JOBS}" && ! -s "${WATCHED_FILES}" ]]; then
    log "Nothing tracked yet (${TRACKED_JOBS}, ${WATCHED_FILES} both empty/missing); nothing to check."
    exit 0
fi

[[ -f "${STATE_FILE}" ]] || echo '{}' > "${STATE_FILE}"

CHANGES=""

update_state() {
    # $1=key $2=value — merges into STATE_FILE (shared by both job-id and file watchers)
    python3 -c "
import json
try:
    d = json.load(open('${STATE_FILE}'))
except Exception:
    d = {}
d['$1'] = '''$2'''
json.dump(d, open('${STATE_FILE}', 'w'))
"
}

get_prev_state() {
    # $1=key
    python3 -c "
import json
try:
    d = json.load(open('${STATE_FILE}'))
except Exception:
    d = {}
print(d.get('$1', ''))
" 2>/dev/null
}

if [[ -s "${TRACKED_JOBS}" ]]; then
while IFS=$'\t' read -r job_id label; do
    [[ -z "${job_id}" || "${job_id}" == \#* ]] && continue

    current=$(ssh "${SSH_OPTS[@]}" "${SSH_HOST}" "qstat -xf ${job_id} 2>&1" 2>/dev/null \
        | grep -E "job_state|Exit_status" | tr -d ' \t' | tr '\n' ';')
    if [[ -z "${current}" ]]; then
        current="GONE"   # left the queue entirely: finished+purged, or never existed
    fi

    prev=$(get_prev_state "job:${job_id}")

    if [[ "${current}" != "${prev}" ]]; then
        CHANGES="${CHANGES}${label} (${job_id}): ${prev:-<none>} -> ${current}\n"
        log "CHANGE ${job_id} (${label}): '${prev:-<none>}' -> '${current}'"
    fi

    update_state "job:${job_id}" "${current}"
done < "${TRACKED_JOBS}"
fi

if [[ -s "${WATCHED_FILES}" ]]; then
while IFS=$'\t' read -r remote_path label; do
    [[ -z "${remote_path}" || "${remote_path}" == \#* ]] && continue

    if ssh "${SSH_OPTS[@]}" "${SSH_HOST}" "test -f '${remote_path}'" 2>/dev/null; then
        current="EXISTS"
    else
        current="MISSING"
    fi

    prev=$(get_prev_state "file:${remote_path}")

    if [[ "${current}" != "${prev}" ]]; then
        CHANGES="${CHANGES}${label} (${remote_path}): ${prev:-<none>} -> ${current}\n"
        log "CHANGE ${remote_path} (${label}): '${prev:-<none>}' -> '${current}'"
    fi

    update_state "file:${remote_path}" "${current}"
done < "${WATCHED_FILES}"
fi

if [[ -n "${CHANGES}" ]]; then
    log "Changes detected, notifying via ntfy:"
    log "${CHANGES}"
    # Phone push via ntfy.sh. The private topic name lives in .ntfy_topic (gitignored, chmod 600); anyone who
    # knows it can read/post, so never commit it. (The old `claude -p` + PushNotification path could not reach the
    # phone: that tool only pushes when a Remote Control session is connected, and it skips when a terminal is active.)
    NTFY_TOPIC_FILE="${SCRIPT_DIR}/.ntfy_topic"
    if [[ -s "${NTFY_TOPIC_FILE}" ]]; then
        NTFY_TOPIC="$(tr -d '[:space:]' < "${NTFY_TOPIC_FILE}")"
        if curl -sS -m 20 --fail -o /dev/null \
            -H "Title: Lab server PBS" -H "Tags: bar_chart" \
            -d "$(printf '%b' "${CHANGES}")" "https://ntfy.sh/${NTFY_TOPIC}" >> "${LOG_FILE}" 2>&1; then
            log "ntfy push sent."
        else
            log "ntfy push FAILED (see curl error above)."
        fi
    else
        log "No ${NTFY_TOPIC_FILE}; cannot push."
    fi
else
    log "No changes."
fi

log "--- check finished ---"
