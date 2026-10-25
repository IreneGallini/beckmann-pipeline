"""
Job dispatcher entrypoint for multi-worker deployments:
`BECKMANN_REDIS_URL=redis://... python dispatcher.py`

Run exactly ONE of these per Redis instance, alongside any number of
gunicorn workers (`gunicorn -w N wsgi:app`). The gunicorn workers only
enqueue jobs and read their status from Redis; this process is the only
one that actually starts, watches, times out and collects the per-job
compute processes. See jobs.py's module docstring.

A separate file rather than `python -m jobs`: running jobs.py itself as
__main__ would give the dispatcher a second, distinct copy of the module
(`__main__` vs `jobs`), with its own process bookkeeping.
"""
import logging
import os
import warnings

warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import jobs
from beckmann_pyscf import settings

if __name__ == "__main__":
    logging.basicConfig(
        level=getattr(logging, settings.LOG_LEVEL, logging.INFO),
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    )
    jobs.run_dispatcher()
