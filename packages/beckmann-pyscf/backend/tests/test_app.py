"""
Tests for app.py's routes: the health endpoint, synchronous input
validation on /submit (no job created for bad input), and rate limiting.
Does not exercise the real predict() pipeline -- see test_jobs.py for the
job-store mechanics and test_pipeline.py for the real (slow) end-to-end run.
"""
import pytest

import app as appmod
import jobs


@pytest.fixture(autouse=True)
def reset_state():
    jobs.JOBS.clear()
    jobs._PROCESSES.clear()
    jobs._STARTED_AT.clear()
    jobs._FINISHED_AT.clear()
    jobs._PENDING.clear()
    appmod.limiter.reset()
    yield
    for proc in jobs._PROCESSES.values():
        if proc.is_alive():
            proc.terminate()
            proc.join(timeout=5)


@pytest.fixture
def client():
    appmod.app.config["TESTING"] = True
    return appmod.app.test_client()


def test_health_returns_ok(client):
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.get_json() == {"status": "ok"}


def test_submit_missing_smiles_returns_400_no_job(client):
    resp = client.post("/submit", json={})
    assert resp.status_code == 400
    assert jobs.JOBS == {}


def test_submit_invalid_smiles_fails_synchronously_no_job_created(client):
    resp = client.post("/submit", json={"smiles": "not a real smiles!!!"})
    assert resp.status_code == 400
    assert "error" in resp.get_json()
    # The whole point of the synchronous pre-check: no job is created at
    # all for input that was never going to work.
    assert jobs.JOBS == {}
    assert jobs._PROCESSES == {}


def test_submit_no_ketone_fails_synchronously(client):
    resp = client.post("/submit", json={"smiles": "CCO"})  # ethanol, no ketone
    assert resp.status_code == 400
    assert jobs.JOBS == {}


def test_rate_limit_returns_429_after_threshold(client):
    # @limiter.limit(settings.RATE_LIMIT) reads settings.RATE_LIMIT once,
    # at decoration time (app.py's own import) -- monkeypatching the
    # setting afterwards wouldn't affect the already-bound decorator, so
    # this test exercises the real configured default directly rather than
    # pretending to reconfigure it.
    from beckmann_pyscf import settings
    limit_count = int(settings.RATE_LIMIT.split()[0])

    for _ in range(limit_count):
        resp = client.post("/submit", json={"smiles": "CCO"})
        assert resp.status_code == 400  # rejected by input validation, but still counts against the rate limit

    resp = client.post("/submit", json={"smiles": "CCO"})
    assert resp.status_code == 429
