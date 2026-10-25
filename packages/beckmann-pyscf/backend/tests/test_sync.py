"""
Tests for sync mode (BECKMANN_API_MODE=sync, the Cloud Run deployment):
jobs.run_blocking()'s heartbeat/timeout/abandon behavior, and the /predict
and /debug/client-ip routes. Like test_jobs.py, uses fast module-level
dummy workers (picklable for the spawn context) instead of the real
multi-minute predict() pipeline.
"""
import json
import multiprocessing
import os
import time

import pytest
from werkzeug.middleware.proxy_fix import ProxyFix

import app as appmod
import jobs
from beckmann_pyscf import settings

TETRALONE = "O=C1CCCC2=CC=CC=C21"


def _done_worker(job_id, smiles, result_queue):
    result_queue.put({"job_id": job_id, "status": "done", "result": {"smiles": smiles}})


def _one_second_worker(job_id, smiles, result_queue):
    import time as _time
    _time.sleep(1)
    result_queue.put({"job_id": job_id, "status": "done", "result": {}})


def _slow_worker(job_id, smiles, result_queue):
    import time as _time
    _time.sleep(30)
    result_queue.put({"job_id": job_id, "status": "done", "result": {}})


def _pid_gone(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return True
    return False


def _only_child_pid():
    children = multiprocessing.active_children()
    assert len(children) == 1, children
    return children[0].pid


# ---- jobs.run_blocking ---------------------------------------------------

def test_run_blocking_returns_the_final_record():
    records = list(jobs.run_blocking("CC(C)=O", heartbeat_s=0.1, target=_done_worker))
    assert records[-1] == {"status": "done", "result": {"smiles": "CC(C)=O"}}
    assert all(r is None for r in records[:-1])


def test_run_blocking_heartbeats_while_running():
    records = list(jobs.run_blocking("CC(C)=O", heartbeat_s=0.2, target=_one_second_worker))
    assert records.count(None) >= 2
    assert records[-1]["status"] == "done"


def test_run_blocking_times_out_and_kills_the_process(monkeypatch):
    monkeypatch.setattr(settings, "JOB_TIMEOUT_SECONDS", 1)
    gen = jobs.run_blocking("CC(C)=O", heartbeat_s=0.2, target=_slow_worker)
    assert next(gen) is None
    pid = _only_child_pid()
    records = list(gen)
    assert records[-1] == {"status": "timeout"}
    assert _pid_gone(pid)


def test_closing_run_blocking_early_kills_the_scan():
    # What happens when a client disconnects mid-scan: the heartbeat write
    # fails and the server closes the generator. The scan must actually
    # stop, not keep running (CPU-throttled, on Cloud Run) in the background.
    gen = jobs.run_blocking("CC(C)=O", heartbeat_s=0.2, target=_slow_worker)
    assert next(gen) is None
    pid = _only_child_pid()
    gen.close()
    assert _pid_gone(pid)


# ---- /predict and mode gating --------------------------------------------

@pytest.fixture(autouse=True)
def reset_state():
    jobs.reset_for_tests()
    appmod.limiter.reset()
    yield
    jobs.reset_for_tests()
    if appmod._SYNC_SLOT.locked():
        appmod._SYNC_SLOT.release()


@pytest.fixture
def client():
    appmod.app.config["TESTING"] = True
    return appmod.app.test_client()


@pytest.fixture
def sync_mode(monkeypatch):
    monkeypatch.setattr(settings, "API_MODE", "sync")


def _fake_run_blocking(records, closed=None):
    def _gen(smiles, heartbeat_s, target=None):
        try:
            yield from records
        finally:
            if closed is not None:
                closed.append(True)
    return _gen


def test_predict_is_disabled_in_async_mode(client):
    assert client.post("/predict", json={"smiles": TETRALONE}).status_code == 404


def test_submit_and_status_are_disabled_in_sync_mode(client, sync_mode):
    assert client.post("/submit", json={"smiles": TETRALONE}).status_code == 404
    assert client.get("/status/whatever").status_code == 404


def test_index_tells_the_frontend_which_mode(client, sync_mode):
    assert b'"sync"' in client.get("/").data


def test_predict_rejects_bad_input_before_streaming(client, sync_mode):
    resp = client.post("/predict", json={"smiles": "CCO"})
    assert resp.status_code == 400
    assert "error" in resp.get_json()
    assert not appmod._SYNC_SLOT.locked()


def test_predict_streams_heartbeats_then_the_record(client, sync_mode, monkeypatch):
    monkeypatch.setattr(jobs, "run_blocking", _fake_run_blocking(
        [None, None, {"status": "done", "result": {"prediction": "R"}}]))
    resp = client.post("/predict", json={"smiles": TETRALONE})
    assert resp.status_code == 200
    body = resp.get_data(as_text=True)
    assert body.startswith("  ")
    assert json.loads(body) == {"status": "done", "result": {"prediction": "R"}}
    assert not appmod._SYNC_SLOT.locked(), "slot must be released after the scan"


def test_predict_returns_503_while_another_scan_runs(client, sync_mode):
    appmod._SYNC_SLOT.acquire()
    resp = client.post("/predict", json={"smiles": TETRALONE})
    assert resp.status_code == 503
    assert "busy" in resp.get_json()["error"]


def test_client_disconnect_closes_the_scan_and_frees_the_slot(client, sync_mode, monkeypatch):
    closed = []
    monkeypatch.setattr(jobs, "run_blocking", _fake_run_blocking(
        [None] * 1000 + [{"status": "done", "result": {}}], closed))
    resp = client.post("/predict", json={"smiles": TETRALONE}, buffered=False)
    first = next(iter(resp.response))
    assert first.strip() == b""
    resp.close()
    assert closed == [True]
    assert not appmod._SYNC_SLOT.locked()


def test_slot_is_freed_even_if_the_stream_never_starts(client, sync_mode, monkeypatch):
    monkeypatch.setattr(jobs, "run_blocking", _fake_run_blocking([{"status": "done", "result": {}}]))
    resp = client.post("/predict", json={"smiles": TETRALONE}, buffered=False)
    resp.close()  # closed before a single chunk was read
    assert not appmod._SYNC_SLOT.locked()


def test_predict_is_rate_limited(client, sync_mode):
    limit_count = int(settings.RATE_LIMIT.split()[0])
    for _ in range(limit_count):
        assert client.post("/predict", json={"smiles": "CCO"}).status_code == 400
    assert client.post("/predict", json={"smiles": "CCO"}).status_code == 429


# ---- /debug/client-ip ------------------------------------------------------

def test_debug_client_ip_is_off_by_default(client):
    assert client.get("/debug/client-ip").status_code == 404


def test_debug_client_ip_keys_on_the_rightmost_forwarded_address(client, monkeypatch):
    # Cloud Run's front end appends the real client IP to whatever
    # X-Forwarded-For the client sent; with TRUSTED_PROXIES=1 the limiter
    # must key on that right-most entry, not a client-forged one.
    monkeypatch.setattr(settings, "DEBUG_CLIENT_IP", True)
    monkeypatch.setattr(settings, "TRUSTED_PROXIES", 1)
    monkeypatch.setattr(appmod.app, "wsgi_app", ProxyFix(appmod.app.wsgi_app, x_for=1))
    resp = client.get(
        "/debug/client-ip",
        headers={"X-Forwarded-For": "1.2.3.4, 203.0.113.7"},
        environ_base={"REMOTE_ADDR": "169.254.1.1"},
    )
    data = resp.get_json()
    assert data["rate_limit_key"] == "203.0.113.7"
    assert data["raw_x_forwarded_for"] == "1.2.3.4, 203.0.113.7"
    assert data["raw_peer_addr"] == "169.254.1.1"
