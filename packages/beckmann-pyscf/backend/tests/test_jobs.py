"""
Tests for jobs.py's job-store mechanics: concurrency cap, queueing, real
process termination on timeout, and expiry sweeping. Uses fast dummy
worker functions (module-level, so multiprocessing's spawn context can
pickle/re-import them in the child) instead of the real multi-minute
predict() pipeline -- these tests exercise the store's own logic, not the
chemistry.

Every test runs twice, once against the in-process MemoryStore and once
against RedisStore (backed by fakeredis), via the parametrized `store`
fixture: both must give identical job-store behavior. The tests call
jobs._reap_once() directly, i.e. they play the dispatcher themselves.
"""
import os
import time

import fakeredis
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


@pytest.fixture(autouse=True, params=["memory", "redis"])
def store(request):
    if request.param == "memory":
        s = jobs.MemoryStore()
    else:
        s = jobs.RedisStore(fakeredis.FakeRedis(decode_responses=True))
    jobs.reset_for_tests(s)
    yield s
    # Clean up any processes a test left running (e.g. a killed test before
    # its own cleanup ran) so they don't linger past the test session.
    jobs.reset_for_tests()


def _submit(smiles, target):
    """submit() plus one dispatcher pass -- with RedisStore, submit() only
    enqueues, and starting the job is the dispatcher's job."""
    job_id = jobs.submit(smiles, target=target)
    jobs._promote()
    return job_id


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
    job_id = _submit("CC(C)=O", target=_done_worker)
    assert jobs.get_status(job_id)["status"] == "running"
    _wait_for_status(job_id, "done")


def test_submit_queues_once_at_the_concurrency_cap(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 1)
    monkeypatch.setattr(settings, "MAX_QUEUE_DEPTH", 5)
    job1 = _submit("CC(C)=O", target=_slow_worker)
    job2 = _submit("CC(C)=O", target=_slow_worker)
    assert jobs.get_status(job1)["status"] == "running"
    assert jobs.get_status(job2)["status"] == "queued"
    # Cleanup: both are slow workers that would otherwise sleep 30s.
    jobs._PROCESSES[job1].terminate()
    jobs._PROCESSES[job1].join(timeout=5)


def test_submit_rejects_once_the_queue_is_full(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 1)
    monkeypatch.setattr(settings, "MAX_QUEUE_DEPTH", 1)
    job1 = _submit("CC(C)=O", target=_slow_worker)
    _submit("CC(C)=O", target=_slow_worker)  # fills the 1-deep queue
    with pytest.raises(RuntimeError):
        _submit("CC(C)=O", target=_slow_worker)
    jobs._PROCESSES[job1].terminate()
    jobs._PROCESSES[job1].join(timeout=5)


def test_timeout_actually_terminates_the_process(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    monkeypatch.setattr(settings, "JOB_TIMEOUT_SECONDS", 1)
    job_id = _submit("CC(C)=O", target=_slow_worker)
    proc = jobs._PROCESSES[job_id]
    assert proc.is_alive()

    time.sleep(1.5)  # let real wall-clock time pass the 1s timeout
    jobs._reap_once()

    assert jobs.get_status(job_id)["status"] == "timeout"
    assert not proc.is_alive(), "process must actually be terminated, not just reported as timed out"


def test_timeout_escalates_to_sigkill_when_process_ignores_sigterm(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    monkeypatch.setattr(settings, "JOB_TIMEOUT_SECONDS", 1)
    job_id = _submit("CC(C)=O", target=_sigterm_ignoring_worker)
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
    job_id = _submit("garbage", target=_invalid_worker)
    job = _wait_for_status(job_id, "invalid_input")
    assert job["error"] == "bad smiles"


def test_queued_job_is_promoted_when_a_slot_frees_up(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 1)
    monkeypatch.setattr(settings, "MAX_QUEUE_DEPTH", 5)
    job1 = _submit("CC(C)=O", target=_done_worker)
    job2 = _submit("CC(C)=O", target=_done_worker)
    assert jobs.get_status(job2)["status"] == "queued"

    _wait_for_status(job1, "done")
    _wait_for_status(job2, "done")


def test_expired_terminal_jobs_are_swept(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    # 1s, not 0: RedisStore expresses retention as a key TTL, which Redis
    # requires to be positive.
    monkeypatch.setattr(settings, "JOB_RETENTION_SECONDS", 1)
    job_id = _submit("CC(C)=O", target=_done_worker)
    _wait_for_status(job_id, "done")

    time.sleep(1.5)
    jobs._reap_once()
    assert jobs.get_status(job_id) is None, "expired job record should have been swept"


def test_running_jobs_are_not_expired(monkeypatch):
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 1)
    monkeypatch.setattr(settings, "JOB_RETENTION_SECONDS", 1)
    job_id = _submit("CC(C)=O", target=_slow_worker)
    time.sleep(1.5)
    jobs._reap_once()
    assert jobs.get_status(job_id)["status"] == "running"


def test_capacity_counts_running_and_queued_jobs_together(monkeypatch, store):
    # Without promotion in between (a dispatcher that hasn't polled yet),
    # submissions still can't exceed MAX_CONCURRENT_JOBS + MAX_QUEUE_DEPTH.
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 1)
    monkeypatch.setattr(settings, "MAX_QUEUE_DEPTH", 1)
    for _ in range(2):
        assert store.try_enqueue(str(time.time()), {"job_id": "x", "smiles": "C", "target": "m:f"}, 2)
    assert not store.try_enqueue("over", {"job_id": "over", "smiles": "C", "target": "m:f"}, 2)
    assert store.get("over") is None


def test_orphaned_running_job_is_reported_as_error_after_restart(store):
    # A dispatcher that died mid-job leaves the job "running" in the store
    # with no process behind it; the next dispatcher must not leave
    # clients polling that forever.
    store.try_enqueue("orphan", {"job_id": "orphan", "smiles": "C", "target": "m:f"}, 10)
    store.pop_pending()
    store.mark_running("orphan")

    jobs.recover_orphans()

    job = jobs.get_status("orphan")
    assert job["status"] == "error"
    assert "restarted" in job["error"]
    assert store.running_ids() == []


def test_web_worker_submit_does_not_start_processes_with_shared_store(monkeypatch, store):
    # With Redis, gunicorn workers must never own job processes -- only
    # the dispatcher does (otherwise a job would die with whichever worker
    # gunicorn happened to recycle).
    if not isinstance(store, jobs.RedisStore):
        pytest.skip("MemoryStore runs the dispatcher in-process by design")
    monkeypatch.setattr(settings, "MAX_CONCURRENT_JOBS", 2)
    job_id = jobs.submit("CC(C)=O", target=_done_worker)
    assert jobs.get_status(job_id)["status"] == "queued"
    assert jobs._PROCESSES == {}
    _wait_for_status(job_id, "done")


def test_concurrent_submits_never_overfill_capacity(store):
    # Many gunicorn workers hitting /submit at once: the capacity check and
    # the enqueue must be atomic, or two racers can both take the last slot.
    import threading

    accepted = []
    barrier = threading.Barrier(20)

    def _racer(i):
        barrier.wait()
        if store.try_enqueue(f"job{i}", {"job_id": f"job{i}", "smiles": "C", "target": "m:f"}, 5):
            accepted.append(i)

    threads = [threading.Thread(target=_racer, args=(i,)) for i in range(20)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert len(accepted) == 5
    assert store.job_count() == 5
