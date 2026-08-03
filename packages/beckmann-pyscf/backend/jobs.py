"""
Job store + executor for the async submit/poll pattern. Each job runs in
its own OS process (multiprocessing, spawn context), not a thread -- only a
process can be forcibly reclaimed on timeout, which is the whole point:
AIMNet2/PySCF computations run on real CPU for real minutes, and a stuck
one has to actually stop, not just get reported as stopped while it keeps
burning CPU in the background.

Single-process-only job store (JOBS/_PROCESSES/_PENDING are plain Python
dicts/deques in this process's memory) -- documented v1 limitation, same as
the original thread-based design. Must run behind exactly one gunicorn
worker (see wsgi.py) until this is replaced with a shared store (Redis)
capable of coordinating across multiple worker processes.

Statuses: "queued" | "running" | "done" | "error" | "invalid_input" | "timeout".
"invalid_input" and "timeout" are distinct from a generic "error" on
purpose -- a bad SMILES and a stuck SCF convergence are different problems
and should read differently to the user.
"""
import collections
import logging
import multiprocessing
import os
import queue as queue_module
import signal
import threading
import time
import uuid

from beckmann_pyscf import settings

logger = logging.getLogger("beckmann_pyscf.jobs")

_ctx = multiprocessing.get_context("spawn")

# job_id -> public status dict, exactly what /status exposes via jsonify.
# Never put a Process/Queue/other non-JSON-serializable object in here.
JOBS: dict[str, dict] = {}

# Internal bookkeeping, kept OUT of JOBS so it never accidentally leaks
# into a jsonify() response.
_PROCESSES: dict[str, multiprocessing.Process] = {}
_STARTED_AT: dict[str, float] = {}
_FINISHED_AT: dict[str, float] = {}
_PENDING: collections.deque[tuple[str, str, object]] = collections.deque()
_RESULT_QUEUE = _ctx.Queue()
_LOCK = threading.Lock()

_TERMINAL_STATUSES = {"done", "error", "invalid_input", "timeout"}

_reaper_started = False


def _worker(job_id: str, smiles: str, result_queue) -> None:
    """Runs IN the child process -- imports beckmann_pyscf.pipeline here
    (not at this module's top level) so merely importing jobs.py in the
    parent doesn't eagerly pull in torch/rdkit/pyscf before a job is even
    submitted. Never touches JOBS/_PROCESSES/etc. directly; the only way
    results get back to the parent is this queue."""
    from beckmann_pyscf.pipeline import predict

    try:
        result = predict(smiles)
        result_queue.put({"job_id": job_id, "status": "done", "result": result})
    except ValueError as e:
        result_queue.put({"job_id": job_id, "status": "invalid_input", "error": str(e)})
    except Exception as e:
        result_queue.put({"job_id": job_id, "status": "error", "error": str(e)})


def _run_in_own_group(target, job_id: str, smiles: str, result_queue) -> None:
    """The actual multiprocessing entrypoint -- wraps `target` so every job
    process becomes the leader of its own new process group (os.setpgrp())
    before anything else runs. A real smoke test found that AIMNet2/torch
    spawn their own internal multiprocessing workers as grandchildren of our
    tracked process -- os.kill()/proc.terminate() aimed at just that one pid
    doesn't reach them, so they survive a "successful" timeout kill. Giving
    the job its own pgid lets the reaper target the whole group at once
    (os.killpg) instead."""
    os.setpgrp()
    target(job_id, smiles, result_queue)


def _start(job_id: str, smiles: str, target=_worker) -> None:
    """Caller must hold _LOCK. `target` defaults to the real _worker (which
    runs the actual predict() pipeline); tests override it with a fast
    dummy so the job-store mechanics (queueing, timeout, expiry) can be
    exercised in milliseconds instead of the real multi-minute pipeline.
    Must be a top-level, picklable function either way -- multiprocessing's
    spawn context re-imports it by (module, name) in the child. daemon=False:
    the real predict() pipeline (Auto3D/AIMNet2) spawns its own child
    processes internally, and Python forbids a daemonic process from having
    children -- confirmed by a real "daemonic processes are not allowed to
    have children" job failure during manual smoke testing. Non-daemon means
    an abnormal (SIGKILL) death of the parent won't auto-reap these, but
    normal shutdown (SIGTERM/gunicorn worker restart) still lets them be
    joined/terminated explicitly, which _reap_once() and timeout handling
    already do."""
    process = _ctx.Process(
        target=_run_in_own_group, args=(target, job_id, smiles, _RESULT_QUEUE), daemon=False
    )
    process.start()
    _PROCESSES[job_id] = process
    _STARTED_AT[job_id] = time.time()
    JOBS[job_id] = {"status": "running"}
    logger.info("job %s started (pid=%s)", job_id, process.pid)


def submit(smiles: str, target=_worker) -> str:
    """Starts (or queues) a background prediction job, returns its job_id.
    Raises RuntimeError if the pending queue is already full -- app.py
    turns that into a 503. `target` is a testing hook -- see _start()."""
    job_id = str(uuid.uuid4())
    with _LOCK:
        if len(_PENDING) >= settings.MAX_QUEUE_DEPTH:
            raise RuntimeError("Server is at capacity, try again shortly")
        if len(_PROCESSES) < settings.MAX_CONCURRENT_JOBS:
            _start(job_id, smiles, target=target)
        else:
            JOBS[job_id] = {"status": "queued"}
            _PENDING.append((job_id, smiles, target))
            logger.info("job %s queued (%d ahead of it)", job_id, len(_PENDING) - 1)
    return job_id


def get_status(job_id: str) -> dict | None:
    """The job's current {"status", ...} dict, or None if job_id is unknown
    (never submitted, or already expired past JOB_RETENTION_SECONDS)."""
    return JOBS.get(job_id)


def _reap_once() -> None:
    """One pass: drain finished results, kill timed-out processes, promote
    queued jobs into freed slots, sweep expired terminal records."""
    while True:
        try:
            msg = _RESULT_QUEUE.get_nowait()
        except queue_module.Empty:
            break
        job_id = msg.pop("job_id")
        with _LOCK:
            JOBS[job_id] = msg
            _FINISHED_AT[job_id] = time.time()
            proc = _PROCESSES.pop(job_id, None)
            _STARTED_AT.pop(job_id, None)
        if proc is not None:
            proc.join(timeout=5)
        logger.info("job %s finished: %s", job_id, msg["status"])

    now = time.time()
    with _LOCK:
        timed_out = [
            job_id for job_id, started in _STARTED_AT.items()
            if now - started > settings.JOB_TIMEOUT_SECONDS
        ]
    for job_id in timed_out:
        with _LOCK:
            proc = _PROCESSES.pop(job_id, None)
            _STARTED_AT.pop(job_id, None)
        if proc is not None:
            pid = proc.pid
            # Target the whole process group (pid == pgid: _run_in_own_group
            # made this job its own group leader via os.setpgrp()), not just
            # the single tracked pid. A real smoke test found AIMNet2/torch's
            # own internal multiprocessing workers -- grandchildren of our
            # tracked process -- survive a plain os.kill()/proc.terminate()
            # aimed at just the direct child. Also don't trust
            # proc.is_alive()'s word before escalating: the same smoke test
            # saw it report False (skipping straight past a kill attempt)
            # while the OS process group was still alive minutes later, a
            # multiprocessing/spawn-context bookkeeping quirk on macOS.
            if pid is not None:
                try:
                    os.killpg(pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            proc.join(timeout=5)
            if pid is not None:
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                except OSError:
                    logger.exception("job %s: os.killpg(%s, SIGKILL) failed", job_id, pid)
            proc.join(timeout=5)
            if proc.is_alive():
                logger.error("job %s timeout cleanup FAILED to terminate pid %s",
                              job_id, pid)
        with _LOCK:
            JOBS[job_id] = {"status": "timeout"}
            _FINISHED_AT[job_id] = time.time()
        logger.warning("job %s timed out after %ss", job_id, settings.JOB_TIMEOUT_SECONDS)

    with _LOCK:
        while _PENDING and len(_PROCESSES) < settings.MAX_CONCURRENT_JOBS:
            job_id, smiles, target = _PENDING.popleft()
            _start(job_id, smiles, target=target)

    now = time.time()
    with _LOCK:
        expired = [
            job_id for job_id, finished in _FINISHED_AT.items()
            if now - finished > settings.JOB_RETENTION_SECONDS
        ]
        for job_id in expired:
            JOBS.pop(job_id, None)
            _FINISHED_AT.pop(job_id, None)


def _reaper_loop(poll_interval: float) -> None:
    while True:
        try:
            _reap_once()
        except Exception:
            logger.exception("reaper pass failed")
        time.sleep(poll_interval)


def start_reaper(poll_interval: float = 2.0) -> None:
    """Starts the single background reaper thread. Idempotent -- app.py
    calls this once at startup; must NOT happen as a module-level import
    side effect, since multiprocessing (spawn) re-executes this module's
    top level in every child worker process too, which would otherwise
    spin up a redundant, useless reaper in each one."""
    global _reaper_started
    if _reaper_started:
        return
    _reaper_started = True
    threading.Thread(target=_reaper_loop, args=(poll_interval,), daemon=True).start()
