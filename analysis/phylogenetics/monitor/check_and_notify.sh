#!/usr/bin/env bash
# check_and_notify.sh — hands-free PBS job monitor for the lab server.
#
# Invoked periodically by the user's cron on the Mint laptop, independent of
# whether an interactive Claude Code session is open — that's what makes this
# "hands-free": it fires on a timer as long as the laptop is awake, no terminal
# needed (the old Windows Task Scheduler setup was retired in the 2026-09-20 migration). See ../../../.claude adjacent notes / the project plan for the full
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
#   - List a glob to watch for its highest "iteration_N" match in watched_globs.tsv
#     (remote_glob<TAB>label), e.g. a Gubbins run dir's `prefix.iteration_*.tre` —
#     for a long-running job whose PBS job_state stays "R" for its whole multi-day
#     run (a tracked_jobs.tsv entry alone is blind to progress within that run, so a
#     Gubbins iteration finishing looked like "nothing changed" — user report,
#     2026-09-25). Reports the highest completed iteration number as its tracked
#     value, so an iteration finishing is a real, notifiable state change and shows
#     up in heartbeats too, independent of which underlying PBS job is running it
#     (survives a kill+resume onto a new job ID, since it just reads the directory).
#   - .state.json (gitignored) holds last-seen status per job_id/file/glob.
#   - .monitor.log (gitignored) is a running log for debugging.
#
# Scheduled via user cron (`crontab -l`): hourly 08:00-22:00 (`0 8-22 * * *`) plus an @reboot run.
# Cron does not fire while the laptop is suspended; on wake it runs the missed job immediately, hence the
# network wait below.
#
# Heartbeat: at 08:00, 11:00, 14:00 and 17:00, send a fixed status push even when nothing changed
# (user request, 2026-09-23, widened to 4 slots same day). Every other hour stays silent-on-no-
# change as before, so this doesn't turn into hourly noise.
#
# Catch-up: Debian/Ubuntu cron only auto-replays a missed slot if the clock jump on wake is <= 3h
# (see cron(8), "Special considerations..."); a longer sleep (e.g. overnight) is treated as a clock
# correction and the missed slot is dropped entirely, no replay. And even a <=3h catch-up runs the
# job at the ACTUAL current time, which may already be a later hour than the target by the time it
# fires. So the heartbeat can't just check "is the current hour one of the slot hours" — it tracks,
# per calendar day, in .state.json (`heartbeat:<date>:<hour>`), whether each slot has been sent yet;
# any check that runs after a slot's hour has passed and finds it unsent (whatever the reason:
# cron's own >3h catch-up gap, or this script's own network-wait delay) sends it then, once,
# however late.
set -uo pipefail   # NOT -e: a single job's ssh/qstat hiccup must not abort the whole check

HOUR="$(date +%H)"
TODAY="$(date +%F)"
HEARTBEAT_HOURS=("08" "11" "14" "17")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRACKED_JOBS="${SCRIPT_DIR}/tracked_jobs.tsv"
WATCHED_FILES="${SCRIPT_DIR}/watched_files.tsv"
WATCHED_GLOBS="${SCRIPT_DIR}/watched_globs.tsv"
STATE_FILE="${SCRIPT_DIR}/.state.json"
LOG_FILE="${SCRIPT_DIR}/.monitor.log"

SSH_KEY="${HOME}/.ssh/id_ed25519_labserver"
SSH_HOST="clima@143.107.143.27"
SSH_PORT="2222"
SSH_OPTS=(-i "${SSH_KEY}" -p "${SSH_PORT}" -o BatchMode=yes -o ConnectTimeout=15 -o StrictHostKeyChecking=accept-new)

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "${LOG_FILE}"; }

log "--- check started ---"

if [[ ! -s "${TRACKED_JOBS}" && ! -s "${WATCHED_FILES}" && ! -s "${WATCHED_GLOBS}" ]]; then
    log "Nothing tracked yet (${TRACKED_JOBS}, ${WATCHED_FILES}, ${WATCHED_GLOBS} all empty/missing); nothing to check."
    exit 0
fi

[[ -f "${STATE_FILE}" ]] || echo '{}' > "${STATE_FILE}"

# After a suspend/resume, cron fires the missed job the instant the lid opens, before Wi-Fi/DNS is back.
# Wait (up to ~2 min) for the network so ssh/ntfy don't fail spuriously.
for _ in $(seq 1 24); do
    getent hosts ntfy.sh >/dev/null 2>&1 && break
    sleep 5
done

CHANGES=""
STATUS_LINES=""   # every tracked item's current state, regardless of change — used for the heartbeat body
PENDING_KEYS=()   # state updates are applied only after the notification succeeded, so a failed push is retried next tick
PENDING_VALS=()

update_state() {
    # $1=key $2=value — merges into STATE_FILE (shared by both job-id and file watchers)
    python3 -c "
import json, sys
try:
    d = json.load(open('${STATE_FILE}'))
except Exception:
    d = {}
d[sys.argv[1]] = sys.argv[2]
json.dump(d, open('${STATE_FILE}', 'w'))
" "$1" "$2"
}

get_prev_state() {
    # $1=key
    python3 -c "
import json, sys
try:
    d = json.load(open('${STATE_FILE}'))
except Exception:
    d = {}
print(d.get(sys.argv[1], ''))
" "$1" 2>/dev/null
}

# Which heartbeat slots are due today (their hour has arrived or passed) and not yet sent today.
CUR_HOUR_INT=$((10#${HOUR}))
DUE_HEARTBEATS=()
for h in "${HEARTBEAT_HOURS[@]}"; do
    h_int=$((10#${h}))
    if (( CUR_HOUR_INT >= h_int )); then
        sent="$(get_prev_state "heartbeat:${TODAY}:${h}")"
        [[ "${sent}" != "sent" ]] && DUE_HEARTBEATS+=("${h}")
    fi
done

if [[ -s "${TRACKED_JOBS}" ]]; then
while IFS=$'\t' read -r job_id label; do
    [[ -z "${job_id}" || "${job_id}" == \#* ]] && continue

    out=$(ssh "${SSH_OPTS[@]}" "${SSH_HOST}" "qstat -xf ${job_id} 2>&1" 2>/dev/null); rc=$?
    if [[ ${rc} -eq 255 ]]; then
        log "ssh FAILED for ${job_id} (network/server unreachable); state unknown, skipping."
        continue
    fi
    current=$(echo "${out}" | grep -E "job_state|Exit_status" | tr -d ' \t' | tr '\n' ';')
    if [[ -z "${current}" ]]; then
        current="GONE"   # left the queue entirely: finished+purged, or never existed
    fi

    prev=$(get_prev_state "job:${job_id}")

    if [[ "${current}" != "${prev}" ]]; then
        CHANGES="${CHANGES}${label} (${job_id}): ${prev:-<none>} -> ${current}\n"
        log "CHANGE ${job_id} (${label}): '${prev:-<none>}' -> '${current}'"
    fi
    STATUS_LINES="${STATUS_LINES}${label}: ${current}\n"

    PENDING_KEYS+=("job:${job_id}"); PENDING_VALS+=("${current}")
done < "${TRACKED_JOBS}"
fi

if [[ -s "${WATCHED_FILES}" ]]; then
while IFS=$'\t' read -r remote_path label; do
    [[ -z "${remote_path}" || "${remote_path}" == \#* ]] && continue

    ssh "${SSH_OPTS[@]}" "${SSH_HOST}" "test -f '${remote_path}'" 2>/dev/null; rc=$?
    if [[ ${rc} -eq 0 ]]; then
        current="EXISTS"
    elif [[ ${rc} -eq 1 ]]; then
        current="MISSING"
    else
        log "ssh FAILED for ${remote_path} (rc=${rc}); state unknown, skipping."
        continue
    fi

    prev=$(get_prev_state "file:${remote_path}")

    if [[ "${current}" != "${prev}" ]]; then
        CHANGES="${CHANGES}${label} (${remote_path}): ${prev:-<none>} -> ${current}\n"
        log "CHANGE ${remote_path} (${label}): '${prev:-<none>}' -> '${current}'"
    fi
    STATUS_LINES="${STATUS_LINES}${label}: ${current}\n"

    PENDING_KEYS+=("file:${remote_path}"); PENDING_VALS+=("${current}")
done < "${WATCHED_FILES}"
fi

if [[ -s "${WATCHED_GLOBS}" ]]; then
while IFS=$'\t' read -r remote_glob label; do
    [[ -z "${remote_glob}" || "${remote_glob}" == \#* ]] && continue

    # highest "iteration_N" among the glob's matches, or empty if none exist yet. The pipeline's
    # last stage (tail) exits 0 even on no matches, so a non-255 ssh rc here just means "checked
    # successfully, found nothing" — not a failure.
    out=$(ssh "${SSH_OPTS[@]}" "${SSH_HOST}" "ls -1 ${remote_glob} 2>/dev/null | grep -oE 'iteration_[0-9]+' | grep -oE '[0-9]+' | sort -n | tail -1" 2>/dev/null); rc=$?
    if [[ ${rc} -eq 255 ]]; then
        log "ssh FAILED for glob ${remote_glob}; state unknown, skipping."
        continue
    fi
    current="${out:-none}"

    prev=$(get_prev_state "iter:${remote_glob}")

    if [[ "${current}" != "${prev}" ]]; then
        CHANGES="${CHANGES}${label}: iteration ${prev:-<none>} -> ${current}\n"
        log "CHANGE ${remote_glob} (${label}): '${prev:-<none>}' -> '${current}'"
    fi
    STATUS_LINES="${STATUS_LINES}${label}: iteration ${current}\n"

    PENDING_KEYS+=("iter:${remote_glob}"); PENDING_VALS+=("${current}")
done < "${WATCHED_GLOBS}"
fi

commit_state() {
    local i
    for i in "${!PENDING_KEYS[@]}"; do update_state "${PENDING_KEYS[$i]}" "${PENDING_VALS[$i]}"; done
}

# Phone push via ntfy.sh. The private topic name lives in .ntfy_topic (gitignored, chmod 600); anyone who
# knows it can read/post, so never commit it. (The old `claude -p` + PushNotification path could not reach the
# phone: that tool only pushes when a Remote Control session is connected, and it skips when a terminal is active.)
push_ntfy() {
    # $1=title $2=tags $3=body -> 0 on success, 1 on failure (caller decides whether that blocks commit_state)
    local title="$1" tags="$2" body="$3"
    local NTFY_TOPIC_FILE="${SCRIPT_DIR}/.ntfy_topic"
    if [[ ! -s "${NTFY_TOPIC_FILE}" ]]; then
        log "No ${NTFY_TOPIC_FILE}; cannot push."
        return 1
    fi
    local NTFY_TOPIC; NTFY_TOPIC="$(tr -d '[:space:]' < "${NTFY_TOPIC_FILE}")"
    if curl -sS -m 20 --retry 3 --retry-delay 5 --retry-all-errors --fail -o /dev/null \
        -H "Title: ${title}" -H "Tags: ${tags}" \
        -d "$(printf '%b' "${body}")" "https://ntfy.sh/${NTFY_TOPIC}" >> "${LOG_FILE}" 2>&1; then
        log "ntfy push sent."
        return 0
    else
        log "ntfy push FAILED (see curl error above)."
        return 1
    fi
}

if [[ -n "${CHANGES}" ]]; then
    log "Changes detected, notifying via ntfy:"
    log "${CHANGES}"
    if push_ntfy "Lab server PBS" "bar_chart" "${CHANGES}"; then
        commit_state
    else
        log "State NOT saved, will retry next check."
    fi
elif [[ ${#DUE_HEARTBEATS[@]} -gt 0 ]]; then
    log "No changes, but heartbeat slot(s) due and unsent today: ${DUE_HEARTBEATS[*]} — sending a fixed status ping."
    HB_LINES=""
    for h in "${DUE_HEARTBEATS[@]}"; do
        if [[ "${h}" == "${HOUR}" ]]; then
            HB_LINES="${HB_LINES}${h}:00 check, on time.\n"
        else
            HB_LINES="${HB_LINES}${h}:00 check, sent late at $(date +%H:%M) — machine was likely asleep.\n"
        fi
    done
    HEARTBEAT="Monitor alive, no changes.\n${HB_LINES}${STATUS_LINES:-(nothing tracked)}"
    if push_ntfy "Lab server PBS — heartbeat" "white_check_mark" "${HEARTBEAT}"; then
        for h in "${DUE_HEARTBEATS[@]}"; do update_state "heartbeat:${TODAY}:${h}" "sent"; done
    else
        log "Heartbeat push failed; that slot stays unsent and will be retried next check."
    fi
    # tracked-item state is independent of the heartbeat outcome — nothing there to retry either way
    commit_state
else
    log "No changes."
    commit_state
fi

log "--- check finished ---"
