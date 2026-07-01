"""
Beckmann rearrangement prediction pipeline.

High-level API for the web app:

    from beckmann import run_prediction
    result = run_prediction("O=C1CCC2=CC=CC=C21")
    # result["conformers"] — list of dicts with name, energy, sdf_block
"""
import os
import tempfile
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from beckmann.oximes import ketone_to_protonated_oximes, enumerate_ez
from beckmann.conformers import generate_conformers
from beckmann.optimize import select_and_optimize
from beckmann.analysis.classical import get_oxime_atoms, predict as _classical_predict


def run_prediction(smiles: str, workdir: Path | None = None) -> dict:
    """Full pipeline: ketone SMILES → optimized conformers ready for 3DMol.js.

    Returns {"conformers": [{"name", "energy_ev", "energy_kcal", "sdf_block"}, ...]}.
    If workdir is None, a temporary directory is created and cleaned up automatically.
    """
    _cleanup = workdir is None
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="beckmann_"))

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"RDKit cannot parse SMILES: {smiles!r}")

        oximes = ketone_to_protonated_oximes(mol)
        if not oximes:
            raise ValueError("No oxime products found — is this a ketone SMILES?")

        smi_path = workdir / "query.smi"
        with open(smi_path, 'w') as f:
            for ox in oximes:
                for iso, ez_label in enumerate_ez(ox):
                    f.write(f"{Chem.MolToSmiles(iso)} query_{ez_label}\n")

        conformers_dir = workdir / "conformers"
        conformers_sdf = generate_conformers(smi_path, conformers_dir)

        opt_dir = workdir / "optimized"
        best_sdf, _ = select_and_optimize(conformers_sdf, opt_dir)

        suppl = Chem.SDMolSupplier(str(best_sdf), removeHs=False)
        conformers = []
        for m in suppl:
            if m is None:
                continue
            energy = float(m.GetProp("E_aimnet2_eV")) if m.HasProp("E_aimnet2_eV") else None

            atom_ids = get_oxime_atoms(m)
            if atom_ids is not None:
                cox, nox, oox, c_aryl, c_allyl = atom_ids
                conf = m.GetConformer()
                d_aryl  = round(abs(rdMolTransforms.GetDihedralDeg(conf, oox, nox, cox, c_aryl)), 1)
                d_allyl = round(abs(rdMolTransforms.GetDihedralDeg(conf, oox, nox, cox, c_allyl)), 1)
                classical_pred = _classical_predict(d_aryl, d_allyl)
            else:
                d_aryl = d_allyl = None
                classical_pred = "inspect"

            conformers.append({
                "name":          m.GetProp("_Name") if m.HasProp("_Name") else "conformer",
                "energy_ev":     energy,
                "energy_kcal":   round(energy * 23.0605, 3) if energy is not None else None,
                "sdf_block":     Chem.MolToMolBlock(m),
                "classical_pred": classical_pred,
                "d_aryl":        d_aryl,
                "d_allyl":       d_allyl,
            })
        return {"conformers": conformers}

    finally:
        if _cleanup:
            import shutil
            shutil.rmtree(workdir, ignore_errors=True)