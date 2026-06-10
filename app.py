import os
import shutil
import subprocess
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")

from pathlib import Path
from flask import Flask, request, jsonify, render_template

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

app = Flask(__name__)
PROJECT_ROOT = Path(__file__).parent


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/predict", methods=["POST"])
def predict():
    data = request.get_json()
    smiles = data.get("smiles", "").strip()
    if not smiles:
        return jsonify({"error": "No SMILES provided"}), 400

    # Write input file
    input_dir = PROJECT_ROOT / "data" / "inputs"
    output_dir = PROJECT_ROOT / "data" / "outputs" / "conformers"
    output_dir.mkdir(parents=True, exist_ok=True)

    smi_path = input_dir / "query.smi"
    with open(smi_path, "w") as f:
        f.write(f"{smiles} query\n")

    # Copy to output dir (Auto3D writes next to input)
    working_copy = output_dir / "query.smi"
    shutil.copy(smi_path, working_copy)

    # Run Auto3D as subprocess so multiprocessing works cleanly
    result = subprocess.run(
        ["python", "scripts/01_smiles_to_conformers.py", "--smiles_file", str(working_copy)],
        capture_output=True, text=True, cwd=str(PROJECT_ROOT)
    )

    if result.returncode != 0:
        return jsonify({"error": result.stderr}), 500

    # Find the output SDF
    sdf_files = sorted(output_dir.glob("query_*/query_out.sdf"))
    if not sdf_files:
        return jsonify({"error": "No output SDF found"}), 500

    sdf_path = sdf_files[-1]  # most recent run

    # Parse SDF and return conformer data
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    conformers = []
    for mol in suppl:
        if mol is None:
            continue
        energy = float(mol.GetProp("E_tot")) if mol.HasProp("E_tot") else None
        conformers.append({
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "conformer",
            "energy_hartree": energy,
            "energy_kcal": round(energy * 627.509, 3) if energy else None,
            "sdf_block": Chem.MolToMolBlock(mol),
        })

    return jsonify({"conformers": conformers})


if __name__ == "__main__":
    app.run(debug=True, port=5001)