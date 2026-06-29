"""
Tests for script 02_select_and_optimize.py output.

Validates: data/output/aimnet_optimized/best_aimnet_optimized.sdf
Checks energy properties, one-structure-per-molecule, and that geometry
was not drastically distorted vs the Auto3D starting conformer.
"""
import math
from rdkit import Chem
from rdkit.Chem import AllChem


def test_aimnet_sdf_exists(aimnet_sdf_path):
    assert aimnet_sdf_path.exists()


def test_one_structure_per_molecule(aimnet_sdf_path):
    """Script 02 should write exactly one optimized structure per molecule."""
    names: list[str] = []
    for mol in Chem.SDMolSupplier(str(aimnet_sdf_path), removeHs=False):
        if mol is not None:
            names.append(mol.GetProp("_Name"))
    duplicates = [n for n in set(names) if names.count(n) > 1]
    assert not duplicates, f"Duplicate molecules in AIMNet SDF: {duplicates}"


def test_energy_properties_present_and_finite(aimnet_sdf_path):
    """E_aimnet2_eV and E_aimnet2_kcal must be present and numerically valid."""
    for mol in Chem.SDMolSupplier(str(aimnet_sdf_path), removeHs=False):
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        for prop in ("E_aimnet2_eV", "E_aimnet2_kcal"):
            assert mol.HasProp(prop), f"{name}: missing property {prop}"
            val = float(mol.GetProp(prop))
            assert math.isfinite(val), f"{name}: {prop}={val} is not finite"


def test_aimnet_covers_all_step1_molecules(aimnet_sdf_path, conformers_sdf_path):
    """Every molecule that came through conformer generation must appear here."""
    step1_names: set[str] = set()
    for mol in Chem.SDMolSupplier(str(conformers_sdf_path), removeHs=False):
        if mol is not None:
            step1_names.add(mol.GetProp("_Name"))

    step2_names: set[str] = set()
    for mol in Chem.SDMolSupplier(str(aimnet_sdf_path), removeHs=False):
        if mol is not None:
            step2_names.add(mol.GetProp("_Name"))

    missing = step1_names - step2_names
    assert not missing, (
        f"Molecules from step 1 missing in AIMNet output: {sorted(missing)}"
    )


def test_rmsd_step1_to_step2(aimnet_sdf_path, conformers_sdf_path):
    """AIMNet2 opt should fine-tune geometry, not rebuild it (RMSD < 2.0 Å)."""
    step1_best: dict[str, tuple[float, Chem.Mol]] = {}
    for mol in Chem.SDMolSupplier(str(conformers_sdf_path), removeHs=False):
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        e = float(mol.GetProp("E_tot"))
        if name not in step1_best or e < step1_best[name][0]:
            step1_best[name] = (e, mol)

    for mol2 in Chem.SDMolSupplier(str(aimnet_sdf_path), removeHs=False):
        if mol2 is None:
            continue
        name = mol2.GetProp("_Name")
        if name not in step1_best:
            continue
        _, mol1 = step1_best[name]
        rmsd = AllChem.GetBestRMS(Chem.Mol(mol1), Chem.Mol(mol2))
        assert rmsd < 2.0, (
            f"{name}: RMSD step1→step2 = {rmsd:.3f} Å — "
            "unusually large, check for optimization failure"
        )