import os
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from flask import Flask, request, jsonify, render_template

from beckmann import run_prediction

app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/predict", methods=["POST"])
def predict():
    data = request.get_json()
    smiles = (data or {}).get("smiles", "").strip()
    if not smiles:
        return jsonify({"error": "No SMILES provided"}), 400

    try:
        result = run_prediction(smiles)
        return jsonify(result)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    app.run(debug=True, port=5001)
