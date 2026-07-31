"""
Flask entrypoint -- adapted from the original beckmann-pipeline app.py's
thread-spawn + JOBS-dict + /submit,/status polling pattern (proven design:
the pipeline takes ~1-2 min, so a synchronous request isn't viable). The only
change is what a job runs (pipeline.predict(), the HPC-free PySCF pipeline,
instead of beckmann.run_prediction()) and that the JOBS/threading logic now
lives in jobs.py rather than inline.

templates/ and static/ are served from the sibling frontend/ directory
(explicit folders below), not from a backend/templates/ -- this repo keeps
backend and frontend as separate top-level directories.
"""
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from flask import Flask, jsonify, render_template, request

import jobs

FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"

app = Flask(
    __name__,
    template_folder=str(FRONTEND_DIR / "templates"),
    static_folder=str(FRONTEND_DIR / "static"),
)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/submit", methods=["POST"])
def submit():
    smiles = (request.get_json() or {}).get("smiles", "").strip()
    if not smiles:
        return jsonify({"error": "No SMILES provided"}), 400
    job_id = jobs.submit(smiles)
    return jsonify({"job_id": job_id})


@app.route("/status/<job_id>")
def status(job_id: str):
    job = jobs.get_status(job_id)
    if job is None:
        return jsonify({"error": "Unknown job"}), 404
    return jsonify(job)


if __name__ == "__main__":
    app.run(debug=True, port=5001)
