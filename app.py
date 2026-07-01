import os
import uuid
import threading
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from flask import Flask, request, jsonify, render_template

from beckmann import run_prediction

app = Flask(__name__)

# job_id → {"status": "running"|"done"|"error", "result": ..., "error": ...}
JOBS: dict[str, dict] = {}


def _run_job(job_id: str, smiles: str) -> None:
    try:
        result = run_prediction(smiles)
        JOBS[job_id] = {"status": "done", "result": result}
    except Exception as e:
        JOBS[job_id] = {"status": "error", "error": str(e)}


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/submit", methods=["POST"])
def submit():
    smiles = (request.get_json() or {}).get("smiles", "").strip()
    if not smiles:
        return jsonify({"error": "No SMILES provided"}), 400
    job_id = str(uuid.uuid4())
    JOBS[job_id] = {"status": "running"}
    threading.Thread(target=_run_job, args=(job_id, smiles), daemon=True).start()
    return jsonify({"job_id": job_id})


@app.route("/status/<job_id>")
def status(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return jsonify({"error": "Unknown job"}), 404
    return jsonify(job)


if __name__ == "__main__":
    app.run(debug=True, port=5001)
