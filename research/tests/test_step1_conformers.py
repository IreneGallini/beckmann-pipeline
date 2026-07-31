"""
Tests for script 01_smiles_to_conformers.py output.

Validates: data/output/conformers/molecules_*/molecules_out.sdf
Checks that all molecules from step 0 are present with valid conformers.
"""
from rdkit import Chem


def test_conformers_sdf_exists(conformers_sdf_path):
    assert conformers_sdf_path.exists()


def test_conformer_count_per_molecule(conformers_sdf_path):
    """Each molecule should have 1–5 conformers (k=5 in Auto3D)."""
    counts: dict[str, int] = {}
    for mol in Chem.SDMolSupplier(str(conformers_sdf_path), removeHs=False):
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        counts[name] = counts.get(name, 0) + 1
    assert counts, "No molecules found in conformers SDF"
    for name, count in counts.items():
        assert 1 <= count <= 5, f"{name}: {count} conformers (expected 1–5)"


def test_conformers_have_energy_property(conformers_sdf_path):
    """Every conformer must carry the E_tot property written by Auto3D."""
    for mol in Chem.SDMolSupplier(str(conformers_sdf_path), removeHs=False):
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        assert mol.HasProp("E_tot"), f"{name}: missing E_tot property"
        e = float(mol.GetProp("E_tot"))
        assert e < 0, f"{name}: E_tot={e} Hartree — expected negative (physical energy)"


def test_conformers_cover_all_step0_molecules(conformers_sdf_path, molecules_smi_path):
    """Every molecule in molecules.smi should appear in the conformers SDF."""
    smi_names: set[str] = set()
    for line in molecules_smi_path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            smi_names.add(parts[1])

    sdf_names: set[str] = set()
    for mol in Chem.SDMolSupplier(str(conformers_sdf_path), removeHs=False):
        if mol is not None:
            sdf_names.add(mol.GetProp("_Name"))

    missing = smi_names - sdf_names
    assert not missing, (
        f"Molecules in molecules.smi missing from conformers SDF: {sorted(missing)}"
    )