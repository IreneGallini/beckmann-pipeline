"""
Deployment-tunable settings, all env-var-driven with CPU-only-hosting
defaults (see CLAUDE.md's Deployment section for the reasoning behind each
default). Nothing here requires a code change to retune for a different
host -- just set the env var.
"""
import os


def _int_env(name: str, default: int) -> int:
    val = os.environ.get(name)
    return int(val) if val else default


# Wall-clock ceiling for one job. Measured directly on a real substrate
# (mol_006, packages/beckmann-pyscf's own scan diagnostics): the 7 PySCF
# SCF calls alone took ~870s once each was moved into its own
# thread-capped subprocess (pyscf_subprocess.py) to work around a
# PyTorch/PySCF native runtime conflict -- plus Auto3D (~1 min) and
# AIMNet2 selection/relaxation on top, a real job on CPU-only hosting
# comfortably needs more than the previous 600s default. 1800s (30 min)
# gives real headroom above the ~950s observed total before we assume
# something is genuinely stuck.
JOB_TIMEOUT_SECONDS = _int_env("BECKMANN_JOB_TIMEOUT", 1800)

# How many jobs run as actual OS processes at once. Each one is a real
# AIMNet2/PySCF workload -- keep this low on a modest CPU instance.
MAX_CONCURRENT_JOBS = _int_env("BECKMANN_MAX_CONCURRENT", 2)

# Once MAX_CONCURRENT_JOBS are running, additional submissions wait in a
# bounded FIFO queue instead of starting immediately. Beyond this depth,
# /submit rejects new work with 503 rather than growing the queue forever.
MAX_QUEUE_DEPTH = _int_env("BECKMANN_MAX_QUEUE", 10)

# How long a finished job's record (done/error/timeout/invalid_input)
# stays in the job store before the reaper sweeps it out. Long enough for
# a user to keep polling or come back later; short enough not to leak
# memory indefinitely on a long-running server.
JOB_RETENTION_SECONDS = _int_env("BECKMANN_JOB_RETENTION", 3600)

# Per-IP rate limit on /submit (flask-limiter syntax, e.g. "5 per hour").
# Every submission triggers real, expensive compute -- this is the
# single-user-abuse guard; MAX_CONCURRENT_JOBS/MAX_QUEUE_DEPTH are the
# separate whole-server guard regardless of how many distinct IPs hit it.
RATE_LIMIT = os.environ.get("BECKMANN_RATE_LIMIT", "5 per hour")

LOG_LEVEL = os.environ.get("BECKMANN_LOG_LEVEL", "INFO")

# Comma-separated list of allowed CORS origins. Empty (the default) means
# same-origin only -- correct today, since frontend/ and backend/ are
# served from the same Flask app. Set this only if the frontend ever moves
# to a different origin/domain than the API.
_origins = os.environ.get("BECKMANN_CORS_ORIGINS", "")
ALLOWED_ORIGINS = [o.strip() for o in _origins.split(",") if o.strip()]

# Shared job store + rate-limit storage, e.g. "redis://127.0.0.1:6379/0".
# Empty (the default) keeps everything in this one process's memory, which
# is fine for local dev but requires `gunicorn -w 1` (see jobs.py). Set it
# to run several gunicorn workers plus one `python dispatcher.py`.
REDIS_URL = os.environ.get("BECKMANN_REDIS_URL", "")

# How many reverse proxies sit in front of the app and append to
# X-Forwarded-For. Must be >= 1 behind a proxy (Cloud Run -- see DEPLOY.md --
# nginx, ...) -- otherwise every request appears to come from the proxy's
# own IP and the per-IP rate limit becomes one limit shared by everyone.
# 0 (the default) trusts no forwarding headers: correct when clients
# connect directly, and the only safe choice then, since the header is
# client-controlled.
TRUSTED_PROXIES = _int_env("BECKMANN_TRUSTED_PROXIES", 0)

# "async" (default): /submit returns a job id at once, the job runs in the
# background, the frontend polls /status. Needs a host that keeps CPU on
# between requests.
# "sync": /predict runs the whole scan inside one streaming HTTP request,
# and /submit + /status are disabled. For hosts that throttle CPU outside
# requests -- Cloud Run with request-based billing (see DEPLOY.md).
API_MODE = os.environ.get("BECKMANN_API_MODE", "async").strip().lower()
if API_MODE not in ("async", "sync"):
    raise ValueError(f"BECKMANN_API_MODE must be 'async' or 'sync', got {API_MODE!r}")

# Seconds between the whitespace heartbeats /predict streams while a scan
# runs: keeps the long connection from being dropped as idle, and makes a
# disconnected client show up (as a failed write) within this long.
SYNC_HEARTBEAT_SECONDS = _int_env("BECKMANN_SYNC_HEARTBEAT", 20)

# Serves GET /debug/client-ip when "1": shows the raw X-Forwarded-For and
# the client IP the rate limiter actually keys on, to verify
# BECKMANN_TRUSTED_PROXIES after a deploy. Off (404) otherwise.
DEBUG_CLIENT_IP = os.environ.get("BECKMANN_DEBUG_CLIENT_IP", "") == "1"

# --- Free-tier usage gate (Cloud Run; see DEPLOY.md, "Staying at $0") ---
# When "1", every /predict first reads this month's billable Cloud Run time
# (beckmann_pyscf/usage.py) and refuses with 503 if this scan, at its
# worst case (JOB_TIMEOUT_SECONDS + startup allowance), could push usage
# past FREE_TIER_FRACTION of the free tier. Fails closed: if usage can't be
# read, no scan runs.
USAGE_GATE = os.environ.get("BECKMANN_USAGE_GATE", "") == "1"
SERVICE_NAME = os.environ.get("K_SERVICE", "beckmann-pyscf")  # Cloud Run sets K_SERVICE
CLOUD_RUN_CPU = float(os.environ.get("BECKMANN_CLOUD_RUN_CPU", "1"))
CLOUD_RUN_MEMORY_GIB = float(os.environ.get("BECKMANN_CLOUD_RUN_MEMORY_GIB", "4"))
FREE_TIER_FRACTION = float(os.environ.get("BECKMANN_FREE_TIER_FRACTION", "0.75"))
# Cold start + shutdown are billed too; reserve this much on top of the
# scan's own worst case.
STARTUP_ALLOWANCE_SECONDS = _int_env("BECKMANN_STARTUP_ALLOWANCE", 120)

# Set on Cloud Run when the page is served elsewhere (GitHub Pages): "/"
# redirects there and /static/* isn't served, so nobody downloads the
# 24 MB Ketcher bundle from Cloud Run (whose outbound data is billed
# beyond 1 GiB/month).
FRONTEND_URL = os.environ.get("BECKMANN_FRONTEND_URL", "").strip()
