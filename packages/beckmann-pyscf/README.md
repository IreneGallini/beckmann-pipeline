# beckmann-pyscf

Open-source, HPC-free Beckmann rearrangement predictor: SMILES → AIMNet2 geometry → PySCF-native wCNmax scan → R/F prediction, with a Flask web UI. No Gaussian, no NBO7, no Citadel — verified by `backend/tests/test_no_hpc_dependency.py`, not just intended.

Depends on `beckmann-core` (pinned version) for the shared oxime/conformer/AIMNet2/geometry primitives and the wCNmax-minimum R/F rule. `backend/beckmann_pyscf/engine/` is the validated PySCF wCNmax computation engine, relocated unchanged from the original prototype.

## Local development

```bash
pip install -e .   # pulls in beckmann-core automatically
cd backend && python app.py   # http://localhost:5001, debug off by default

python -m pytest tests/ backend/ -v -m "not slow" -c backend/pytest.ini
```

`python app.py` is for local development only. Set `FLASK_DEBUG=1` to opt into Werkzeug's interactive debugger locally — never set it in production (it allows arbitrary code execution from the browser on an unhandled exception).

## Deployment

Each job (Auto3D + AIMNet2 + a 7-point PySCF scan) runs in its own OS process, not a thread — this is what lets a stuck/slow job actually be terminated on timeout instead of just being reported as timed out while it keeps burning CPU. See `backend/jobs.py`'s module docstring for the full design.

**Run via gunicorn, not `python app.py`:**

```bash
gunicorn -w 1 -b 0.0.0.0:5001 wsgi:app
```

**`-w 1` is required, not a tunable default.** The job store (`JOBS`/
`_PROCESSES`/`_PENDING` in `jobs.py`) is plain in-process Python memory —
more than one gunicorn worker would each have its own independent copy, so
a `/status/<job_id>` poll would 404 whenever it landed on a different
worker than the one that handled the matching `/submit`. A single worker
can still comfortably handle many concurrent `/submit`/`/status` requests,
since those return immediately — the actual heavy compute concurrency is
bounded separately by `BECKMANN_MAX_CONCURRENT` (below), not by the
gunicorn worker count. Going beyond one worker (or one host) requires
replacing the in-memory store with something shared across processes
(Redis, a database) — not built yet, no hosting decision has been made to
build it against.

### Environment variables

All optional, all have CPU-only-hosting defaults (see `beckmann_pyscf/settings.py`):

| Variable | Default | What it controls |
|---|---|---|
| `BECKMANN_JOB_TIMEOUT` | `1800` (30 min) | Wall-clock ceiling per job before its process is terminated and the job reports `timeout` |
| `BECKMANN_MAX_CONCURRENT` | `2` | How many jobs run as real OS processes at once |
| `BECKMANN_MAX_QUEUE` | `10` | How many additional jobs wait in a FIFO queue before `/submit` starts returning 503 |
| `BECKMANN_JOB_RETENTION` | `3600` (1 hr) | How long a finished job's record stays pollable before being swept |
| `BECKMANN_RATE_LIMIT` | `"5 per hour"` | Per-IP limit on `/submit` (flask-limiter syntax) |
| `BECKMANN_LOG_LEVEL` | `INFO` | Python `logging` level |
| `BECKMANN_CORS_ORIGINS` | *(empty = same-origin only)* | Comma-separated allowed origins, only needed if the frontend is ever served from a different domain than this API |

### Health check

`GET /health` → `{"status": "ok"}`, no rate limit, no compute — the target for whatever hosting platform's liveness/readiness probe ends up in front of this.

### Docker

```bash
# from the monorepo root
docker build -f packages/beckmann-pyscf/backend/Dockerfile -t beckmann-pyscf .
docker run -p 5001:5001 --memory=4g --cpus=2 beckmann-pyscf
```

Size `--memory`/`--cpus` (or the equivalent resource requests/limits on
whatever container platform you use) to comfortably fit
`BECKMANN_MAX_CONCURRENT` simultaneous AIMNet2/PySCF jobs, not just one —
each is a real, memory-hungry scientific computation. **Not built/tested
in the environment this Dockerfile was written in** (no Docker available
there) — build and run it yourself before trusting it.

### Not yet decided (yours to make)

Hosting provider, domain/HTTPS, and the Redis migration needed to run more
than one worker process. See the monorepo root `CLAUDE.md`'s Deployment
section.
