"""
Real, slow end-to-end smoke test for generate_conformers()/select_and_optimize()
-- actually runs Auto3D + AIMNet2 on one small molecule. (The original
tests/test_step1_conformers.py and test_step2_aimnet.py validated the
34-molecule benchmark batch's file output specifically; that batch and its
tests now live in research/benchmark_pipeline/.)
"""
import math

import pytest
from rdkit import Chem

from beckmann_core.conformers import generate_conformers
from beckmann_core.optimize import select_and_optimize


@pytest.mark.slow
def test_conformers_then_optimize_end_to_end(tmp_path):
    smi_path = tmp_path / "query.smi"
    smi_path.write_text("CC(C)=N[OH2+] query_noEZ\n")  # protonated acetoxime cation

    conformers_sdf = generate_conformers(smi_path, tmp_path / "conformers", k=3)
    assert conformers_sdf.exists()

    conf_mols = [m for m in Chem.SDMolSupplier(str(conformers_sdf), removeHs=False) if m is not None]
    assert conf_mols
    assert 1 <= len(conf_mols) <= 3

    best_sdf, best_per_substrate_sdf = select_and_optimize(conformers_sdf, tmp_path / "optimized")
    assert best_sdf.exists() and best_per_substrate_sdf.exists()

    (opt_mol,) = [m for m in Chem.SDMolSupplier(str(best_per_substrate_sdf), removeHs=False) if m is not None]
    energy = float(opt_mol.GetProp("E_aimnet2_eV"))
    assert math.isfinite(energy)
