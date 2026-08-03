"""
gunicorn entrypoint: `gunicorn -w 1 -b 0.0.0.0:5001 wsgi:app`

MUST run as exactly one worker (`-w 1`). jobs.py's job store (JOBS/
_PROCESSES/_PENDING) lives in plain Python memory in whichever process
imports it -- with more than one gunicorn worker, each worker would have
its own independent copy, so a `/status/<job_id>` poll would 404 whenever
it happened to land on a different worker than the one that handled the
matching `/submit`. This is safe to run with one worker despite being a
prediction service: the actual heavy compute happens in a separate
`multiprocessing.Process` per job (see jobs.py), not in the request-handling
worker itself, so a single gunicorn worker can still comfortably handle
many concurrent lightweight `/submit`/`/status` requests. The real
concurrency limit on heavy compute is `settings.MAX_CONCURRENT_JOBS`, not
the gunicorn worker count.

Going beyond one worker (or beyond one host) requires replacing jobs.py's
in-memory store with something shared across processes (Redis, a database)
-- not done here, no infra decision has been made yet to build against.
"""
from app import app

if __name__ == "__main__":
    app.run()
