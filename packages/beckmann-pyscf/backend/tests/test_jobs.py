"""
Tests for jobs.py's job-store mechanics: concurrency cap, queueing, real
process termination on timeout, and expiry sweeping. Uses fast dummy
worker functions (module-level, so multiprocessing's spawn context can
pickle/re-import them in the child) instead of the real multi-minute
predict() pipeline -- these tests exercise the store's own logic, not the
chemistry.

Each test resets jobs.py's module-level state first (reset_jobs_state
fixture) since JOBS/_PROCESSES/_PENDING/etc. are plain module globals
shared across the whole test session.
"""
import os
import time

import pytest

import jobs
from beckmann_pyscf import settings


def _done_worker(job_id, smiles, result_queue):
    result_queue.put({"job_id": job_id, "status": "done", "result": {"smiles": smiles}})


def _invalid_worker(job_id, smiles, result_queue):
    result_queue.put({"job_id": job_id, "status": "invalid_input", "error": "bad smiles"})


def _slow_worker(job_id, smiles, result_queue):
    import time as _time
    _time.sleep(30)
    result_queue.put({"job_id": job_id, "status": "done", "result": {}})


def _sigterm_ignoring_worker(job_id, smiles, result_queue):
    # Regression target: a real smoke test found that proc.terminate()'s
    # SIGTERM alone can leave the process alive (and, separately, that
    # proc.is_alive() can report False without the process actually being
    # gone) -- this worker ignores SIGTERM outright, forcing the reaper's
    # SIGKILL escalation to be the only thing that can end it.
    import signal as _signal
    import time as _time
    _signal.signal(_signal.SIGTERM, _signal.SIG_IGN)
    _time.sleep(30)
    result_queue.put({"job_id": job_id, "status": "done", "result": {}})


@pytest.fixture(autouse=True)
def reset_jobs_state():
    jobs.JOBS.clear()
    jobs._PROCESSES.clear()
    jobs._STARTED_AT.clear()
    jobs._FINISHED_AT.clear()
    jobs._PENDING.clear()
    yield
    # Clean up any processes a test left running (e.g. a killed test before
    # its own cleanup ran) so they don't linger past the test session.
    for proc in jobs._PROCESSES.values():
        if proc.is_alive():
            proc.terminate()
            proc.join(timeout=5)


def _wait_for_status(job_id, want, timeout=10):
    """Poll jobs.get_status via manual _reap_once() calls (deterministic,
    no dependency on the background reaper thread's own sleep timing)."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        jobs._reap_once()
        job = jobs.get_status(job_id)
        if job and job["status"] == want:
            return job
        time.sleep(0.1)
    raise AssertionError(f"job {job_id} never reached status {want!r}, last: {jobs.get_status(job_id)}")


def test_submit_runs_immediately_under_the_concurrency_cap(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    job_id = jobs.submit("CC(C)=O", target=_done_worker)
    assert jobs.get_status(job_id)["status"] == "running"
    _wait_for_status(job_id, "done")


def test_submit_queues_once_at_the_concurrency_cap(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 1)
    monkeypatch.setattr(settings, "MAX_QUEUE_DEPTH", 5)
    job1 = jobs.submit("CC(C)=O", target=_slow_worker)
    job2 = jobs.submit("CC(C)=O", target=_slow_worker)
    assert jobs.get_status(job1)["status"] == "running"
    assert jobs.get_status(job2)["status"] == "queued"
    # Cleanup: both are slow workers that would otherwise sleep 30s.
    jobs._PROCESSES[job1].terminate()
    jobs._PROCESSES[job1].join(timeout=5)


def test_submit_rejects_once_the_queue_is_full(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 1)
    monkeypatch.setattr(settings, "MAX_QUEUE_DEPTH", 1)
    job1 = jobs.submit("CC(C)=O", target=_slow_worker)
    jobs.submit("CC(C)=O", target=_slow_worker)  # fills the 1-deep queue
    with pytest.raises(RuntimeError):
        jobs.submit("CC(C)=O", target=_slow_worker)
    jobs._PROCESSES[job1].terminate()
    jobs._PROCESSES[job1].join(timeout=5)


def test_timeout_actually_terminates_the_process(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    monkeypatch.setattr(settings, "JOB_TIMEOUT_SECONDS", 1)
    job_id = jobs.submit("CC(C)=O", target=_slow_worker)
    proc = jobs._PROCESSES[job_id]
    assert proc.is_alive()

    time.sleep(1.5)  # let real wall-clock time pass the 1s timeout
    jobs._reap_once()

    assert jobs.get_status(job_id)["status"] == "timeout"
    assert not proc.is_alive(), "process must actually be terminated, not just reported as timed out"


def test_timeout_escalates_to_sigkill_when_process_ignores_sigterm(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    monkeypatch.setattr(settings, "JOB_TIMEOUT_SECONDS", 1)
    job_id = jobs.submit("CC(C)=O", target=_sigterm_ignoring_worker)
    pid = jobs._PROCESSES[job_id].pid

    time.sleep(1.5)
    jobs._reap_once()

    assert jobs.get_status(job_id)["status"] == "timeout"
    # os.kill(pid, 0) sends no real signal, just probes liveness -- it
    # raises ProcessLookupError once the pid is actually gone from the OS's
    # perspective. This is the real invariant that matters (not proc.
    # is_alive(), which a real smoke test caught giving a false negative).
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)


def test_invalid_input_status_is_distinct_from_error(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    job_id = jobs.submit("garbage", target=_invalid_worker)
    job = _wait_for_status(job_id, "invalid_input")
    assert job["error"] == "bad smiles"


def test_queued_job_is_promoted_when_a_slot_frees_up(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 1)
    monkeypatch.setattr(settings, "MAX_QUEUE_DEPTH", 5)
    job1 = jobs.submit("CC(C)=O", target=_done_worker)
    job2 = jobs.submit("CC(C)=O", target=_done_worker)
    assert jobs.get_status(job2)["status"] == "queued"

    _wait_for_status(job1, "done")
    _wait_for_status(job2, "done")


def test_expired_terminal_jobs_are_swept(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    job_id = jobs.submit("CC(C)=O", target=_done_worker)
    _wait_for_status(job_id, "done")  # observe it complete before shortening retention

    monkeypatch.setattr(settings, "JOB_RETENTION_SECONDS", 0)
    time.sleep(0.1)
    jobs._reap_once()
    assert jobs.get_status(job_id) is None, "expired job record should have been swept"
