"""
Flask entrypoint. templates/ and static/ are served from the sibling
frontend/ directory (explicit folders below), not from a backend/templates/
-- this repo keeps backend and frontend as separate top-level directories.

Production note: run this behind gunicorn via wsgi.py, not `python app.py`
directly (see wsgi.py's own docstring for the single-worker constraint).
`python app.py` is for local development only.
"""
import json
import logging
import os
import threading
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from flask import Flask, Response, abort, jsonify, redirect, render_template, request, stream_with_context
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from werkzeug.middleware.proxy_fix import ProxyFix

import jobs
from beckmann_pyscf import settings, usage
from beckmann_pyscf.pipeline import validate_smiles

logging.basicConfig(
    level=getattr(logging, settings.LOG_LEVEL, logging.INFO),
    format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    force=True,  # one of flask/flask_limiter's imports above installs a
    # default-format root handler first; without force=True, basicConfig()
    # silently no-ops and every log line falls back to the bare
    # "LEVEL:name:message" format instead of this one.
)
logger = logging.getLogger("beckmann_pyscf.app")

FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"

app = Flask(
    __name__,
    template_folder=str(FRONTEND_DIR / "templates"),
    # Not served at all when the page lives elsewhere (BECKMANN_FRONTEND_URL,
    # i.e. GitHub Pages): Cloud Run's outbound data is billed past 1 GiB a
    # month, and the Ketcher bundle alone is 24 MB per visit.
    static_folder=None if settings.FRONTEND_URL else str(FRONTEND_DIR / "static"),
)

# Same-origin only unless BECKMANN_CORS_ORIGINS is set -- correct today
# since frontend/ and backend/ are served from this one Flask app.
if settings.ALLOWED_ORIGINS:
    CORS(app, origins=settings.ALLOWED_ORIGINS)

# Behind a reverse proxy (Cloud Run's front end, nginx), the client's real IP
# only arrives via X-Forwarded-For -- without this, get_remote_address
# would see the proxy's IP and rate-limit every user as one.
if settings.TRUSTED_PROXIES:
    app.wsgi_app = ProxyFix(app.wsgi_app, x_for=settings.TRUSTED_PROXIES, x_proto=settings.TRUSTED_PROXIES)

# Counts live in Redis when it's configured, so the per-IP limit holds
# across every gunicorn worker rather than per worker process.
limiter = Limiter(
    get_remote_address,
    app=app,
    default_limits=[],
    storage_uri=settings.REDIS_URL or "memory://",
)

# Background jobs only exist in async mode; sync mode runs each scan
# inside its own /predict request instead.
if settings.API_MODE == "async":
    jobs.start_reaper()

# One scan at a time per process in sync mode (run gunicorn with -w 1, so
# this is one per container). Non-blocking: a second /predict gets a 503
# instead of queueing behind a ~20-30 min scan.
_SYNC_SLOT = threading.Lock()


@app.route("/")
def index():
    if settings.FRONTEND_URL:
        return redirect(settings.FRONTEND_URL, code=302)
    # api_base "": the page is served by this same app, so it calls the
    # API with relative URLs. (frontend/build_static.py renders the same
    # template for GitHub Pages with the Cloud Run URL instead.)
    return render_template("index.html", api_mode=settings.API_MODE, api_base="")


@app.route("/health")
def health():
    """Liveness/readiness probe target for whatever hosting this ends up
    behind. Deliberately cheap (no compute, no rate limit) and separate
    from the prediction API. Reports 503 if the shared job store is
    configured but unreachable, since /submit can't work then either."""
    if not jobs.get_store().ping():
        return jsonify({"status": "degraded", "error": "job store unreachable"}), 503
    return jsonify({"status": "ok"})


def _validated_smiles():
    """(smiles, None) or (None, error_response). Fast synchronous
    pre-check -- fail in milliseconds on obviously-bad input instead of
    spending minutes of compute first (no job is even started for a
    rejected SMILES)."""
    # force=True: the GitHub Pages frontend posts its JSON as text/plain, a
    # CORS "simple request", so the browser skips the OPTIONS preflight
    # round trip (one fewer billed Cloud Run request per prediction).
    smiles = (request.get_json(force=True, silent=True) or {}).get("smiles", "").strip()
    if not smiles:
        return None, (jsonify({"error": "No SMILES provided"}), 400)
    validation_error = validate_smiles(smiles)
    if validation_error is not None:
        logger.info("rejected invalid SMILES %r: %s", smiles, validation_error)
        return None, (jsonify({"error": validation_error}), 400)
    return smiles, None


@app.route("/submit", methods=["POST"])
@limiter.limit(settings.RATE_LIMIT, methods=["POST"])
def submit():
    if settings.API_MODE != "async":
        abort(404)
    smiles, error = _validated_smiles()
    if error is not None:
        return error

    try:
        job_id = jobs.submit(smiles)
    except RuntimeError as e:
        logger.warning("submission rejected, server at capacity: %s", e)
        return jsonify({"error": str(e)}), 503

    logger.info("job %s submitted (smiles=%r)", job_id, smiles)
    return jsonify({"job_id": job_id})


@app.route("/status/<job_id>")
def status(job_id: str):
    if settings.API_MODE != "async":
        abort(404)
    job = jobs.get_status(job_id)
    if job is None:
        return jsonify({"error": "Unknown job"}), 404
    return jsonify(job)


def _usage_gate_error() -> str | None:
    """With BECKMANN_USAGE_GATE=1: None if one more worst-case scan fits
    in this month's free-tier allowance, else the message to show. Fails
    closed -- if usage can't be read, no scan runs."""
    if not settings.USAGE_GATE:
        return None
    try:
        token, project = usage.metadata_token()
        used = usage.month_usage(token, project, settings.SERVICE_NAME)
    except Exception:
        logger.exception("usage gate: could not read this month's usage, refusing the scan")
        return "Predictions are temporarily unavailable, please try again later."
    reason = usage.check(
        used, settings.CLOUD_RUN_CPU, settings.CLOUD_RUN_MEMORY_GIB, settings.FREE_TIER_FRACTION,
        reserve_seconds=settings.JOB_TIMEOUT_SECONDS + settings.STARTUP_ALLOWANCE_SECONDS,
    )
    if reason is not None:
        logger.warning("usage gate: refusing scan, %s", reason)
        return ("This month's free compute allowance for this site is used up. "
                "Predictions resume on the 1st of next month (UTC).")
    logger.info("usage gate: ok (%.0f billable s used this month)", used["billable_seconds"])
    return None


@app.route("/predict", methods=["POST"])
@limiter.limit(settings.RATE_LIMIT, methods=["POST"])
def predict():
    """Sync mode: runs the whole scan inside this request. The response
    streams one space every SYNC_HEARTBEAT_SECONDS while the scan runs,
    then the final record ({"status": "done", "result": ...} or
    error/invalid_input/timeout) as JSON -- leading whitespace is valid
    JSON, so the client just parses the whole body. The status code is
    committed (200) before the scan starts, so scan failures are reported
    in the body's "status", not the HTTP code; only pre-scan rejections
    (400 bad input, 429 rate limit, 503 busy) use HTTP codes.

    The heartbeats double as disconnect detection: once the client is gone
    a write fails, the generator is closed, and jobs.run_blocking() kills
    the scan's process group."""
    if settings.API_MODE != "sync":
        abort(404)
    smiles, error = _validated_smiles()
    if error is not None:
        return error
    gate_error = _usage_gate_error()
    if gate_error is not None:
        return jsonify({"error": gate_error}), 503
    if not _SYNC_SLOT.acquire(blocking=False):
        return jsonify({"error": "The server is busy with another prediction, please try again later."}), 503

    logger.info("sync prediction started (smiles=%r)", smiles)

    released = threading.Event()

    def _release_slot():
        # Idempotent: runs from the stream's finally (after the scan is
        # dead) and again from call_on_close, which also covers a response
        # closed before the stream ever started -- an unstarted generator's
        # finally never runs, which would otherwise leak the slot forever.
        if not released.is_set():
            released.set()
            _SYNC_SLOT.release()

    def _stream():
        scan = jobs.run_blocking(smiles, heartbeat_s=settings.SYNC_HEARTBEAT_SECONDS)
        try:
            for record in scan:
                yield " " if record is None else json.dumps(record)
        finally:
            scan.close()  # kills the scan's process group if it's still running
            _release_slot()

    response = Response(
        stream_with_context(_stream()),
        mimetype="text/plain",
        # Ask any buffering proxy to pass the heartbeats through as they come.
        headers={"X-Accel-Buffering": "no", "Cache-Control": "no-store"},
    )
    response.call_on_close(_release_slot)
    return response


@app.route("/debug/client-ip")
def debug_client_ip():
    """Only when BECKMANN_DEBUG_CLIENT_IP=1: shows what the rate limiter
    keys on, to check BECKMANN_TRUSTED_PROXIES after a deploy (DEPLOY.md).
    Reveals nothing beyond the caller's own request."""
    if not settings.DEBUG_CLIENT_IP:
        abort(404)
    orig = request.environ.get("werkzeug.proxy_fix.orig", {})
    return jsonify({
        "rate_limit_key": get_remote_address(),
        "raw_x_forwarded_for": orig.get("HTTP_X_FORWARDED_FOR", request.headers.get("X-Forwarded-For")),
        "raw_peer_addr": orig.get("REMOTE_ADDR", request.environ.get("REMOTE_ADDR")),
        "trusted_proxies": settings.TRUSTED_PROXIES,
    })


if __name__ == "__main__":
    # debug=True is a real security hole in production (Werkzeug's
    # interactive debugger allows code execution from the browser on an
    # unhandled exception) -- never the default. Opt into it locally with
    # FLASK_DEBUG=1; production should use wsgi.py + gunicorn regardless.
    debug = os.environ.get("FLASK_DEBUG", "") == "1"
    app.run(debug=debug, port=5001)
