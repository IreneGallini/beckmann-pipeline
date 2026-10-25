"""
gunicorn entrypoint. How to run it depends on the mode (see jobs.py's and
entrypoint.sh's docstrings for the full design):

- **Sync mode** (`BECKMANN_API_MODE=sync`, the Docker image / Cloud Run
  default): `gunicorn -w 1 --worker-class gthread --threads 4 -b 0.0.0.0:8080 wsgi:app`.
  Each scan runs inside its own streaming /predict request; one worker
  process so the one-scan-at-a-time lock covers the whole container, and
  gthread so long requests aren't killed by gunicorn's --timeout.

- **Async mode, no Redis** (the `python app.py` default):
  `gunicorn -w 1 -b 0.0.0.0:5001 wsgi:app`. Exactly one worker: the job
  store lives in that worker's memory, so with more, a `/status/<job_id>`
  poll would 404 whenever it landed on a different worker than the
  matching `/submit`. The dispatcher runs as a thread inside that worker.

- **Async mode with Redis**: `python dispatcher.py` (exactly one) plus
  `gunicorn -w N -b 0.0.0.0:5001 wsgi:app` (any N). Job records, the
  pending queue and rate-limit counts all live in Redis.

In every mode the heavy compute runs in a separate `multiprocessing.Process`
per job (its own process group), never in the request-handling process
itself.
"""
from app import app

if __name__ == "__main__":
    app.run()
