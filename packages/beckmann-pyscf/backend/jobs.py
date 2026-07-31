"""
In-memory job store + background-thread runner for the async submit/poll
pattern. Behavior is unchanged from the original app.py's inline JOBS dict
(same single-process, in-memory store -- not multi-worker safe; documented
v1 limitation, no heavier task queue needed unless load testing says
otherwise). Pulled into its own module so app.py stays a thin route file.
"""
import threading
import uuid

from beckmann_pyscf.pipeline import predict

# job_id -> {"status": "running"|"done"|"error", "result": ..., "error": ...}
JOBS: dict[str, dict] = {}


def _run_job(job_id: str, smiles: str) -> None:
    try:
        result = predict(smiles)
        JOBS[job_id] = {"status": "done", "result": result}
    except Exception as e:
        JOBS[job_id] = {"status": "error", "error": str(e)}


def submit(smiles: str) -> str:
    """Starts a background prediction job, returns its job_id."""
    job_id = str(uuid.uuid4())
    JOBS[job_id] = {"status": "running"}
    threading.Thread(target=_run_job, args=(job_id, smiles), daemon=True).start()
    return job_id


def get_status(job_id: str) -> dict | None:
    """The job's current {"status", ...} dict, or None if job_id is unknown."""
    return JOBS.get(job_id)
