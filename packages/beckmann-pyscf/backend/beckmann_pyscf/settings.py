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
