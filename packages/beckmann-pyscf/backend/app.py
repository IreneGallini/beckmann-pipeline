"""
Flask entrypoint. templates/ and static/ are served from the sibling
frontend/ directory (explicit folders below), not from a backend/templates/
-- this repo keeps backend and frontend as separate top-level directories.

Production note: run this behind gunicorn via wsgi.py, not `python app.py`
directly (see wsgi.py's own docstring for the single-worker constraint).
`python app.py` is for local development only.
"""
import logging
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from flask import Flask, jsonify, render_template, request
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

import jobs
from beckmann_pyscf import settings
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
    static_folder=str(FRONTEND_DIR / "static"),
)

# Same-origin only unless BECKMANN_CORS_ORIGINS is set -- correct today
# since frontend/ and backend/ are served from this one Flask app.
if settings.ALLOWED_ORIGINS:
    CORS(app, origins=settings.ALLOWED_ORIGINS)

limiter = Limiter(get_remote_address, app=app, default_limits=[])

jobs.start_reaper()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/health")
def health():
    """Liveness/readiness probe target for whatever hosting this ends up
    behind. Deliberately cheap (no compute, no rate limit) and separate
    from the prediction API."""
    return jsonify({"status": "ok"})


@app.route("/submit", methods=["POST"])
@limiter.limit(settings.RATE_LIMIT)
def submit():
    smiles = (request.get_json() or {}).get("smiles", "").strip()
    if not smiles:
        return jsonify({"error": "No SMILES provided"}), 400

    # Fast synchronous pre-check -- fail in milliseconds on obviously-bad
    # input instead of spending minutes of compute first (a job is not
    # even created for a rejected SMILES).
    validation_error = validate_smiles(smiles)
    if validation_error is not None:
        logger.info("rejected invalid SMILES %r: %s", smiles, validation_error)
        return jsonify({"error": validation_error}), 400

    try:
        job_id = jobs.submit(smiles)
    except RuntimeError as e:
        logger.warning("submission rejected, server at capacity: %s", e)
        return jsonify({"error": str(e)}), 503

    logger.info("job %s submitted (smiles=%r)", job_id, smiles)
    return jsonify({"job_id": job_id})


@app.route("/status/<job_id>")
def status(job_id: str):
    job = jobs.get_status(job_id)
    if job is None:
        return jsonify({"error": "Unknown job"}), 404
    return jsonify(job)


if __name__ == "__main__":
    # debug=True is a real security hole in production (Werkzeug's
    # interactive debugger allows code execution from the browser on an
    # unhandled exception) -- never the default. Opt into it locally with
    # FLASK_DEBUG=1; production should use wsgi.py + gunicorn regardless.
    debug = os.environ.get("FLASK_DEBUG", "") == "1"
    app.run(debug=debug, port=5001)
