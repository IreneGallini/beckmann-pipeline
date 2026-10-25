# beckmann-pyscf

Open-source, HPC-free Beckmann rearrangement predictor: SMILES → AIMNet2 geometry → PySCF-native wCNmax scan → R/F prediction, with a Flask web UI. No Gaussian, no NBO7, no Citadel: verified by `backend/tests/test_no_hpc_dependency.py`, not just intended.

Depends on `beckmann-core` (pinned version) for the shared oxime/conformer/AIMNet2/geometry primitives and the wCNmax-minimum R/F rule. `backend/beckmann_pyscf/engine/` is the validated PySCF wCNmax computation engine, relocated unchanged from the original prototype.

## Install the CLI

No clone needed. Python 3.11+, in a fresh virtualenv or conda env:

```bash
pip install \
  "beckmann-core @ git+https://github.com/IreneGallini/beckmann-pipeline.git@v0.1.0#subdirectory=packages/beckmann-core" \
  "beckmann-pyscf @ git+https://github.com/IreneGallini/beckmann-pipeline.git@v0.1.0#subdirectory=packages/beckmann-pyscf"

beckmann-pyscf --help
```

Both packages come from the same tagged commit (`v0.1.0`). `beckmann-core` has to be listed explicitly, because it isn't on PyPI, so pip can't resolve `beckmann-pyscf`'s `beckmann-core==0.1.0` requirement by itself. To track the latest development version instead, replace `@v0.1.0` with `@main` in both URLs. On CPU-only machines, run `pip install torch==2.12.0 --index-url https://download.pytorch.org/whl/cpu` first to avoid downloading the much larger CUDA build.

## `beckmann-pyscf` CLI

The recommended way to run this locally, once installed (see **Install the CLI** above, or `pip install -e .` from a clone; either registers the `beckmann-pyscf` command). The Flask web app documented further below stays available separately as a hosted prototype for external collaborators without a local Python environment: if you have this repo cloned, use the CLI.

```bash
beckmann-pyscf predict --smiles "O=C1CCC2=C1C=CC=C2" --name test1 --plot
# [1/4]..[4/4] progress -> ./beckmann_pyscf_runs/test1/{optimized.sdf,wcnmax_series.csv,summary.txt,wcnmax_vs_rno.png}

beckmann-pyscf predict --csv molecules.csv --plot   # 'id'/'SMILES' columns, same shape as data/input/benchmark.csv;
                                                     # one bad row is skipped (printed to stderr), doesn't abort the batch

beckmann-pyscf conformers --smiles "O=C1CCC2=C1C=CC=C2" --name test1   # SMILES -> Auto3D conformers only
beckmann-pyscf optimize --conformers-sdf <path from above>              # conformers SDF -> AIMNet2-optimized geometry only
beckmann-pyscf scan --sdf <path from above> --plot                      # optimized SDF -> PySCF wCNmax scan + prediction only
```

`predict` runs the full pipeline in one call (or, with `--csv`, once per row). `conformers`/`optimize`/`scan` run one stage at a time, each printing its own result (energy, resolved oxime atom map, etc.) so a stage's output can be inspected before trusting a full `predict` run: useful for catching a bad atom-map resolution or an unexpected conformer before it silently propagates downstream. Every subcommand takes an explicit input path from the previous stage's printed output (no implicit state tracking), and defaults its own `--out` to `./beckmann_pyscf_runs/<name>/...` when not given (per `--csv` row, `<name>` is that row's `id`). `scan` also accepts `--ci`/`--ni`/`--oi`/`--c-aryl`/`--c-alkyl` to override the auto-detected atom map, and `--r-min`/`--r-max`/`--r-step` to adjust the scan window/resolution.

Full walkthrough with worked examples: the monorepo root `INSTRUCTIONS.md`.

## Local development

```bash
pip install -e ../beckmann-core -e '.[test]'   # [test] adds pytest + fakeredis
cd backend && python app.py   # http://localhost:5001, debug off by default, in-process job store

python -m pytest tests/ backend/ -v -m "not slow" -c backend/pytest.ini
```

`python app.py` is for local development only. Set `FLASK_DEBUG=1` to opt into Werkzeug's interactive debugger locally, never set it in production (it allows arbitrary code execution from the browser on an unhandled exception).

## Deployment

**Hosted for $0:** the API runs on Google Cloud Run, the page on GitHub Pages, and the image on Docker Hub. A per-scan usage gate and an hourly guard keep Cloud Run inside its free tier. Setup, the $0 safeguards and verification are in the monorepo root's [`DEPLOY.md`](../../DEPLOY.md). The app has two API modes, picked by `BECKMANN_API_MODE`, because they need different things from the host.

### Sync mode (`BECKMANN_API_MODE=sync`: the Docker image and Cloud Run)

```
 browser ── POST /predict ──▶ gunicorn -w 1 (gthread) ──▶ one OS process (own process group)
         ◀── " " every 20 s … then the final JSON record ──     Auto3D → AIMNet2 → PySCF scan
```

- The whole scan runs **inside one HTTP request** (`jobs.run_blocking()`). That's what Cloud Run's request-based billing needs, because it throttles CPU outside requests, so background work would stall.
- `/predict` validates the SMILES synchronously, so bad input gets a plain 400. Then it streams a space every `BECKMANN_SYNC_HEARTBEAT` seconds, then the final record (`{"status": "done", "result": …}` or `error`/`invalid_input`/`timeout`). Leading whitespace is valid JSON, so the client parses the whole body. Because the HTTP status is committed (200) before the scan starts, scan failures are reported in `status`.
- **The heartbeats do two jobs.** They keep a 20–30-minute connection from looking idle. And once the client disconnects, the next write fails and the scan's process group is killed, so no orphaned scan keeps running.
- **One scan per container:** a second `/predict` gets 503 "busy" immediately. Run with `-w 1 --worker-class gthread`: gthread's `--timeout` is a worker-liveness check, not a per-request limit, while a plain sync worker would be killed mid-scan.
- `/submit` and `/status` return 404 in this mode.

### Async mode (`BECKMANN_API_MODE=async`: `python app.py`, other hosts)

```
            HTTP                  enqueue / read status
 browser ─────────▶ gunicorn ×N ───────────────────────▶  Redis
                    (wsgi.py)                              ▲  job records (TTL), pending queue,
                                                           │  running set, rate-limit counts
                                     pop pending / write   │
                    dispatcher.py ─────────────────────────┘
                    (exactly one)
                         │ one OS process (own process group) per job,
                         ▼ at most BECKMANN_MAX_CONCURRENT at once
                    Auto3D → AIMNet2 → PySCF scan
```

`/submit` returns a job id at once, the job runs in the background, and the frontend polls `/status/<id>`. **This mode needs a host that keeps CPU on between requests**, so on Cloud Run it would need instance-based billing.
- **gunicorn workers** only call `jobs.submit()`/`jobs.get_status()`, which read and write Redis. Any worker can answer any `/status` poll.
- **`dispatcher.py`** is the only process that starts, watches, times out (`os.killpg` on the job's process group) and collects job processes. Run exactly one per Redis instance. When it starts, any job still marked `running` from a previous dispatcher is reported as `error` ("server restarted"). On SIGTERM it kills its running jobs before exiting.
- **Capacity**: `/submit` returns 503 once `BECKMANN_MAX_CONCURRENT` jobs are running and `BECKMANN_MAX_QUEUE` more are waiting. The check-and-enqueue is a single Redis transaction (WATCH/MULTI), so racing workers can't overfill the queue.
- **Without `BECKMANN_REDIS_URL`** the same code uses an in-process store, with the dispatcher as a thread in the web process. That **requires `gunicorn -w 1`**, because with more workers a `/status` poll 404s whenever it lands on a different worker than the matching `/submit`.

Status values (both modes): `queued` (async only) → `running` → `done` | `error` | `invalid_input` | `timeout`.

### Running it

```bash
# Sync mode, exactly as in the container / on Cloud Run
PORT=8080 BECKMANN_API_MODE=sync gunicorn -w 1 --worker-class gthread --threads 4 -b 0.0.0.0:8080 wsgi:app

# Async, multi-worker (any Redis: local redis-server, managed, ...)
export BECKMANN_REDIS_URL=redis://127.0.0.1:6379/0
python dispatcher.py &                                  # exactly one
gunicorn -w 3 -b 0.0.0.0:5001 wsgi:app                  # any number of workers

# Async, single process, no Redis
gunicorn -w 1 -b 0.0.0.0:5001 wsgi:app
```

### Environment variables

All optional (see `beckmann_pyscf/settings.py`). Cloud Run's values are set in `.github/workflows/deploy-cloud-run.yml`.

| Variable | Default | What it controls |
|---|---|---|
| `BECKMANN_API_MODE` | `async` (`sync` in the Docker image) | `sync`: scan inside the `/predict` request. `async`: `/submit` plus `/status` polling (see above) |
| `BECKMANN_REDIS_URL` | *(empty = in-process store, `-w 1` only)* | Async mode's shared job store **and** rate-limit storage, e.g. `redis://127.0.0.1:6379/0`. Setting it also makes `entrypoint.sh` start redis-server and the dispatcher |
| `BECKMANN_TRUSTED_PROXIES` | `0` | Number of reverse proxies whose `X-Forwarded-For` entries to trust. **`1` on Cloud Run** (reasoning and a check: `DEPLOY.md`). Leave at `0` when clients connect directly, since the header is client-controlled |
| `BECKMANN_DEBUG_CLIENT_IP` | *(off)* | `1` serves `GET /debug/client-ip` (the raw XFF and the IP the rate limiter keys on), for verifying `BECKMANN_TRUSTED_PROXIES`. Turn it off again afterwards |
| `BECKMANN_JOB_TIMEOUT` | `1800` (30 min) | Wall-clock ceiling per scan before its process group is killed and the job reports `timeout`. Cloud Run uses `3300`, under its 3600 s request timeout |
| `BECKMANN_SYNC_HEARTBEAT` | `20` | Seconds between `/predict` heartbeats |
| `BECKMANN_MAX_CONCURRENT` | `2` | Async mode: how many jobs run as real OS processes at once. Sync mode always runs one per container |
| `BECKMANN_MAX_QUEUE` | `10` | Async mode: how many additional jobs can wait before `/submit` returns 503 |
| `BECKMANN_JOB_RETENTION` | `3600` (1 hr) | Async mode: how long a finished job's record stays pollable |
| `BECKMANN_RATE_LIMIT` | `"5 per hour"` | Per-IP limit on `/submit` and `/predict` (flask-limiter syntax) |
| `BECKMANN_LOG_LEVEL` | `INFO` | Python `logging` level |
| `BECKMANN_CORS_ORIGINS` | *(empty = same-origin only)* | Comma-separated allowed origins, only needed if the frontend is ever served from a different domain than this API |
| `BECKMANN_USAGE_GATE` | *(off)* | `1` (Cloud Run): before each scan, read this month's billable Cloud Run time from Cloud Monitoring, and refuse (503) if a worst-case scan could exceed the free-tier allowance. Fails closed |
| `BECKMANN_CLOUD_RUN_CPU` / `BECKMANN_CLOUD_RUN_MEMORY_GIB` | `1` / `4` | Instance size the allowance is computed from. Must match Cloud Run's `--cpu`/`--memory` |
| `BECKMANN_FREE_TIER_FRACTION` | `0.75` | Share of Cloud Run's monthly free tier the service may use |
| `BECKMANN_STARTUP_ALLOWANCE` | `120` | Seconds reserved per scan for billed cold start/shutdown |
| `BECKMANN_FRONTEND_URL` | *(empty)* | When set (Cloud Run: the GitHub Pages URL), `/` redirects there and no static files are served, so Cloud Run's billed outbound data stays tiny |
| `PYSCF_MAX_MEMORY` | `4000` (PySCF's own default; `2000` in the image) | PySCF's per-calculation memory budget in MB. Caps in-core integral storage so the SCF fits next to torch/AIMNet2 |
| `PORT` | `8080` (image) | Port gunicorn binds. Cloud Run sets it |
| `WEB_CONCURRENCY` | `2` (image) | Async + Redis layout only: gunicorn worker count |
| `BECKMANN_LOCAL_REDIS` | `1` (image) | Async + Redis layout only: `0` = don't start the in-container redis-server (use an external `BECKMANN_REDIS_URL`) |

### Health check

`GET /health` returns `{"status": "ok"}`, or 503 `{"status": "degraded"}` if Redis is configured but unreachable. It has no rate limit, does no compute, and stays responsive during a sync-mode scan.

### Docker

```bash
# from the monorepo root
docker build -f packages/beckmann-pyscf/backend/Dockerfile -t beckmann-pyscf .
docker run -p 8080:8080 --memory=4g --cpus=1 beckmann-pyscf
# -> http://localhost:8080   (sync mode, same as Cloud Run)
```

The image installs CPU-only torch, includes the AIMNet2 weights at build time (so no job downloads anything), and runs as an unprivileged user. `entrypoint.sh` starts the sync layout by default, or the async Redis layout when `BECKMANN_REDIS_URL` is set.

### Still open

Custom domain (the `*.run.app` URL is HTTPS out of the box) and GPU hosting. GPU would also need `generate_conformers()`'s hardcoded `use_gpu=False` revisited and the timeout defaults re-tuned.
