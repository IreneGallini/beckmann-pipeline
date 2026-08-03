"""
predict(smiles) -> dict: SMILES -> protonated oxime -> Auto3D conformers ->
AIMNet2 optimization -> PySCF-native wCNmax rigid scan -> R/F prediction.

No Gaussian/NBO7/Citadel anywhere in this call path -- see
backend/tests/test_no_hpc_dependency.py, which asserts that as a tested
property rather than just an intention.

Depends on beckmann-core (pip-installed, pinned version) for the shared
oxime/conformer/AIMNet2/geometry primitives and the wCNmax-minimum rule, and
on beckmann_pyscf.engine (this package's own PySCF wCNmax computation,
relocated unchanged from beckmann_alt/pair_nbo.py).
"""
import os
import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"  # required before Auto3D/AIMNet2 import on macOS

from rdkit import Chem

from beckmann_core.classical import get_oxime_atoms

from beckmann_pyscf.pipeline.conformers import parse_and_check_ketone, run_conformers, smiles_to_oxime_smi
from beckmann_pyscf.pipeline.optimize import run_optimize
from beckmann_pyscf.pipeline.predict import predict_outcome
from beckmann_pyscf.pipeline.wcnmax_pyscf import run_scan_series


def validate_smiles(smiles: str) -> str | None:
    """Fast synchronous pre-check: does this SMILES parse and contain a
    ketone/oxime-convertible group? Returns an error message string if
    not, None if it's fine to proceed. Backed by the exact same check
    predict() itself runs first (via parse_and_check_ketone()) -- calling
    this before submitting a job lets a bad SMILES fail in milliseconds
    instead of after minutes of Auto3D/AIMNet2/PySCF compute."""
    try:
        parse_and_check_ketone(smiles)
        return None
    except ValueError as e:
        return str(e)


def predict(smiles: str, workdir: Path | None = None) -> dict:
    """Full HPC-free pipeline for one SMILES. Returns:
        {
          "smiles": str,
          "sdf_block": str,             # AIMNet2-optimized winning geometry, MOL-block text
          "energy_ev": float | None,
          "wcnmax_series": [{"stage", "R_NO", "weight", "MO_index"}, ...],
          "prediction": "R" | "F",
          "minimum": dict | None,       # find_wcnmax_minimum()'s result, None if no interior minimum
        }
    Raises ValueError for an unparseable SMILES or a SMILES with no
    ketone/oxime-convertible group. Any other failure (Auto3D/AIMNet2/PySCF)
    propagates unwrapped -- caught generically by jobs.py's background job.
    """
    _cleanup = workdir is None
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="beckmann_pyscf_"))

    try:
        smi_path = smiles_to_oxime_smi(smiles, workdir)

        conformers_dir = workdir / "conformers"
        conformers_sdf = run_conformers(smi_path, conformers_dir)

        opt_dir = workdir / "optimized"
        best_mol = run_optimize(conformers_sdf, opt_dir)

        atom_ids = get_oxime_atoms(best_mol)
        if atom_ids is None:
            raise ValueError(
                "Could not resolve the oxime C=N-O / aryl / alkyl atom map on the optimized structure"
            )
        cox, nox, oox, c_aryl, c_allyl = atom_ids  # 0-based RDKit indices

        rows = run_scan_series(
            best_mol,
            ci=cox + 1, ni=nox + 1, oi=oox + 1, c_aryl=c_aryl + 1, c_alkyl=c_allyl + 1,
            name="query",
        )

        prediction, minimum = predict_outcome("query", rows)
        energy = float(best_mol.GetProp("E_aimnet2_eV")) if best_mol.HasProp("E_aimnet2_eV") else None

        return {
            "smiles": smiles,
            "sdf_block": Chem.MolToMolBlock(best_mol),
            "energy_ev": energy,
            "wcnmax_series": [
                {"stage": r["stage"], "R_NO": r["R_NO"], "weight": r["weight"], "MO_index": r["MO_index"]}
                for r in rows
            ],
            "prediction": prediction,
            "minimum": minimum,
        }
    finally:
        if _cleanup:
            import shutil
            shutil.rmtree(workdir, ignore_errors=True)
