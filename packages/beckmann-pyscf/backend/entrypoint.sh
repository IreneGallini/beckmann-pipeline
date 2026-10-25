#!/usr/bin/env bash
# Container entrypoint. Two layouts, picked by BECKMANN_REDIS_URL:
#
# - unset (the image default; Cloud Run): a single gunicorn process with
#   threads. In sync mode (BECKMANN_API_MODE=sync, also the image default)
#   each scan runs inside its own /predict request, guarded by a
#   per-process lock -- so exactly one worker process (-w 1). gthread, not
#   the sync worker class: gthread's --timeout is a worker-liveness check,
#   not a per-request limit, so a 30-minute streaming request is fine; a
#   sync worker would be killed after --timeout seconds mid-scan.
#
# - set: the async multi-worker layout -- redis-server (unless
#   BECKMANN_LOCAL_REDIS=0, i.e. an external Redis), one job dispatcher,
#   and gunicorn -w WEB_CONCURRENCY. Only for hosts that keep CPU on
#   between requests; pair it with BECKMANN_API_MODE=async.
#
# In both, if any child exits the whole container exits (wait -n), so the
# platform restarts everything together. This script is PID 1: it forwards
# SIGTERM (Cloud Run sends one before stopping an instance) to its children.
set -euo pipefail

PORT="${PORT:-8080}"
WEB_CONCURRENCY="${WEB_CONCURRENCY:-2}"

pids=()

if [[ -n "${BECKMANN_REDIS_URL:-}" ]]; then
    if [[ "${BECKMANN_LOCAL_REDIS:-1}" == "1" ]]; then
        redis-server --bind 127.0.0.1 --port 6379 \
            --save "" --appendonly no \
            --maxmemory "${BECKMANN_REDIS_MAXMEMORY:-256mb}" --maxmemory-policy noeviction \
            --dir /tmp &
        pids+=($!)
        for _ in $(seq 1 50); do
            python -c "import redis, os; redis.Redis.from_url(os.environ['BECKMANN_REDIS_URL']).ping()" 2>/dev/null && break
            sleep 0.2
        done
    fi

    python dispatcher.py &
    pids+=($!)

    gunicorn -w "$WEB_CONCURRENCY" -b "0.0.0.0:$PORT" \
        --timeout 60 --graceful-timeout 30 \
        --access-logfile - wsgi:app &
    pids+=($!)
else
    gunicorn -w 1 --worker-class gthread --threads 4 -b "0.0.0.0:$PORT" \
        --timeout 120 --graceful-timeout 30 \
        --access-logfile - wsgi:app &
    pids+=($!)
fi

shutdown() {
    kill -TERM "${pids[@]}" 2>/dev/null || true
    wait || true
}
trap shutdown TERM INT

status=0
wait -n || status=$?
echo "entrypoint: a child process exited (status $status), stopping container" >&2
shutdown
exit "$status"
