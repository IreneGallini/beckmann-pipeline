"""
Tests for the $0 safeguards (DEPLOY.md, "Staying at $0"): the free-tier
math and Monitoring-API reader in beckmann_pyscf/usage.py, the app's
per-scan usage gate on /predict, the GitHub Pages redirect, and the
hourly guard script's open/close decisions. No network: the metadata
server, Monitoring API and gcloud are all faked.
"""
import datetime as dt
import importlib.util
import io
import json
from pathlib import Path

import pytest

import app as appmod
import jobs
from beckmann_pyscf import settings, usage

TETRALONE = "O=C1CCCC2=CC=CC=C21"


# ---- free-tier math --------------------------------------------------------

def test_allowance_is_bound_by_whichever_free_tier_runs_out_first():
    # 1 vCPU / 4 GiB: memory binds (360000/4 = 90000 s < 180000 s).
    assert usage.monthly_limit_seconds(1, 4, 0.75) == 67500
    # 2 vCPU / 4 GiB: CPU and memory both bind at 90000 s.
    assert usage.monthly_limit_seconds(2, 4, 1.0) == 90000
    # 2 vCPU / 2 GiB: CPU binds.
    assert usage.monthly_limit_seconds(2, 2, 1.0) == 90000


def test_month_start_is_utc_first_of_month():
    now = dt.datetime(2026, 11, 17, 13, 5, 9, tzinfo=dt.timezone.utc)
    assert usage.month_start(now) == dt.datetime(2026, 11, 1, tzinfo=dt.timezone.utc)


def test_check_reserves_the_worst_case_scan():
    ok = {"billable_seconds": 60000, "egress_bytes": 0}
    assert usage.check(ok, 1, 4, 0.75) is None
    assert usage.check(ok, 1, 4, 0.75, reserve_seconds=3420) is None   # 63420 <= 67500
    assert usage.check(ok, 1, 4, 0.75, reserve_seconds=7600) is not None  # 67600 > 67500


def test_check_also_caps_outbound_data():
    heavy = {"billable_seconds": 0, "egress_bytes": 0.8 * usage.FREE_EGRESS_BYTES}
    assert "egress" in usage.check(heavy, 1, 4, 0.75)


def test_sum_metric_adds_all_points_across_pages(monkeypatch):
    pages = [
        {"timeSeries": [{"points": [{"value": {"doubleValue": 100.5}},
                                    {"value": {"doubleValue": 20}}]}],
         "nextPageToken": "p2"},
        {"timeSeries": [{"points": [{"value": {"int64Value": "300"}}]}]},
    ]
    seen = []

    def fake_urlopen(req, timeout):
        seen.append(req.full_url)
        assert req.headers["Authorization"] == "Bearer tok"
        return io.BytesIO(json.dumps(pages[len(seen) - 1]).encode())

    monkeypatch.setattr(usage.urllib.request, "urlopen", fake_urlopen)
    start = dt.datetime(2026, 11, 1, tzinfo=dt.timezone.utc)
    total = usage._sum_metric("tok", "proj", "run.googleapis.com/x", "svc", start, start)
    assert total == 420.5
    assert len(seen) == 2 and "pageToken=p2" in seen[1]
    assert "service_name+%3D+%22svc%22" in seen[0]


# ---- the app's per-scan gate -------------------------------------------------

@pytest.fixture(autouse=True)
def reset_state():
    jobs.reset_for_tests()
    appmod.limiter.reset()
    yield
    if appmod._SYNC_SLOT.locked():
        appmod._SYNC_SLOT.release()


@pytest.fixture
def client():
    appmod.app.config["TESTING"] = True
    return appmod.app.test_client()


@pytest.fixture
def gated(monkeypatch):
    monkeypatch.setattr(settings, "API_MODE", "sync")
    monkeypatch.setattr(settings, "USAGE_GATE", True)
    monkeypatch.setattr(usage, "metadata_token", lambda: ("tok", "proj"))
    started = []

    def fake_run_blocking(smiles, heartbeat_s, target=None):
        started.append(smiles)
        yield {"status": "done", "result": {}}

    monkeypatch.setattr(jobs, "run_blocking", fake_run_blocking)
    return started


def _usage(monkeypatch, billable_seconds):
    monkeypatch.setattr(usage, "month_usage",
                        lambda token, project, service: {"billable_seconds": billable_seconds,
                                                         "egress_bytes": 0})


def test_gate_allows_a_scan_that_fits(client, gated, monkeypatch):
    _usage(monkeypatch, 1000)
    resp = client.post("/predict", json={"smiles": TETRALONE})
    assert resp.status_code == 200
    assert gated == [TETRALONE]


def test_gate_refuses_a_scan_that_could_exceed_the_allowance(client, gated, monkeypatch):
    # 67500 s allowance at 1 vCPU / 4 GiB / 0.75; a worst-case scan
    # reserves JOB_TIMEOUT + startup allowance on top of usage so far.
    _usage(monkeypatch, 67500 - settings.JOB_TIMEOUT_SECONDS)
    resp = client.post("/predict", json={"smiles": TETRALONE})
    assert resp.status_code == 503
    assert "free compute allowance" in resp.get_json()["error"]
    assert gated == [], "no scan may start"
    assert not appmod._SYNC_SLOT.locked()


def test_gate_fails_closed_when_usage_is_unreadable(client, gated, monkeypatch):
    def boom(*a):
        raise OSError("metadata server unreachable")
    monkeypatch.setattr(usage, "month_usage", boom)
    resp = client.post("/predict", json={"smiles": TETRALONE})
    assert resp.status_code == 503
    assert "temporarily unavailable" in resp.get_json()["error"]
    assert gated == []


def test_bad_input_is_rejected_before_the_gate_spends_an_api_call(client, gated, monkeypatch):
    calls = []
    monkeypatch.setattr(usage, "month_usage", lambda *a: calls.append(a))
    assert client.post("/predict", json={"smiles": "CCO"}).status_code == 400
    assert calls == []


def test_predict_accepts_json_sent_as_text_plain(client, gated, monkeypatch):
    # The GitHub Pages frontend sends text/plain to avoid a CORS preflight.
    _usage(monkeypatch, 0)
    resp = client.post("/predict", data=json.dumps({"smiles": TETRALONE}),
                       headers={"Content-Type": "text/plain"})
    assert resp.status_code == 200
    assert gated == [TETRALONE]


def test_index_redirects_to_github_pages_when_configured(client, monkeypatch):
    monkeypatch.setattr(settings, "FRONTEND_URL", "https://example.github.io/beckmann-pipeline/")
    resp = client.get("/")
    assert resp.status_code == 302
    assert resp.headers["Location"] == "https://example.github.io/beckmann-pipeline/"


# ---- the hourly guard ------------------------------------------------------------

def _load_guard():
    path = Path(__file__).resolve().parents[4] / ".github/scripts/free_tier_guard.py"
    spec = importlib.util.spec_from_file_location("free_tier_guard", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _FakeGcloud:
    def __init__(self, public: bool, exists: bool = True):
        self.public, self.exists, self.calls = public, exists, []

    def __call__(self, *args, check=True):
        self.calls.append(args)

        class R:
            returncode, stdout = 0, ""
        r = R()
        if args[:3] == ("run", "services", "describe"):
            r.returncode = 0 if self.exists else 1
        elif args[:2] == ("auth", "print-access-token"):
            r.stdout = "tok\n"
        elif args[:3] == ("run", "services", "get-iam-policy"):
            members = ["allUsers"] if self.public else []
            r.stdout = json.dumps({"bindings": [{"role": "roles/run.invoker", "members": members}]})
        return r

    def changed(self):
        return [c[2] for c in self.calls if len(c) > 2 and c[2].endswith("iam-policy-binding")]


@pytest.fixture
def guard(monkeypatch):
    mod = _load_guard()
    for k, v in {"PROJECT_ID": "proj", "REGION": "us-central1", "SERVICE": "beckmann-pyscf",
                 "CLOUD_RUN_CPU": "1", "CLOUD_RUN_MEMORY_GIB": "4", "FREE_TIER_FRACTION": "0.75"}.items():
        monkeypatch.setenv(k, v)
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    return mod


def _run_guard(guard, monkeypatch, public, billable_seconds, exists=True):
    fake = _FakeGcloud(public, exists)
    monkeypatch.setattr(guard, "gcloud", fake)
    monkeypatch.setattr(guard.usage, "month_usage",
                        lambda t, p, s: {"billable_seconds": billable_seconds, "egress_bytes": 0})
    assert guard.main() == 0
    return fake.changed()


def test_guard_closes_public_access_near_the_allowance(guard, monkeypatch):
    # within one hour (the time until the next check) of 67500 s
    assert _run_guard(guard, monkeypatch, public=True, billable_seconds=64000) == ["remove-iam-policy-binding"]


def test_guard_leaves_a_healthy_public_service_alone(guard, monkeypatch):
    assert _run_guard(guard, monkeypatch, public=True, billable_seconds=10000) == []


def test_guard_reopens_once_a_new_month_resets_usage(guard, monkeypatch):
    assert _run_guard(guard, monkeypatch, public=False, billable_seconds=0) == ["add-iam-policy-binding"]


def test_guard_keeps_a_closed_service_closed_while_over(guard, monkeypatch):
    assert _run_guard(guard, monkeypatch, public=False, billable_seconds=70000) == []


def test_guard_is_a_no_op_before_the_first_deploy(guard, monkeypatch):
    assert _run_guard(guard, monkeypatch, public=False, billable_seconds=0, exists=False) == []
