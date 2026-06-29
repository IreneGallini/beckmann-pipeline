"""
Tests for script 00_benchmark_to_oximes.py output.

Validates: data/input/molecules.smi and data/input/benchmark_meta.json
"""
import re
from rdkit import Chem

OXIME_PAT = Chem.MolFromSmarts("[C]=[N]-[O+]")
MOL_NAME_RE = re.compile(r"^mol_\d{3}_(E|Z)$")
EXPECTED_MOL_COUNT = 34


def _smi_lines(molecules_smi_path):
    return [
        line.split()
        for line in molecules_smi_path.read_text().splitlines()
        if line.strip()
    ]


def test_molecules_smi_nonempty(molecules_smi_path):
    assert len(_smi_lines(molecules_smi_path)) > 0


def test_molecules_smi_name_format(molecules_smi_path):
    """Every name must follow mol_XXX_E or mol_XXX_Z."""
    for parts in _smi_lines(molecules_smi_path):
        name = parts[1]
        assert MOL_NAME_RE.match(name), f"Unexpected name format: {name!r}"


def test_molecules_smi_valid_smiles(molecules_smi_path):
    for parts in _smi_lines(molecules_smi_path):
        smi, name = parts[0], parts[1]
        assert Chem.MolFromSmiles(smi) is not None, f"{name}: RDKit cannot parse SMILES {smi!r}"


def test_molecules_smi_protonated_oxime(molecules_smi_path):
    """Every SMILES must contain the protonated C=N-[O+] motif."""
    for parts in _smi_lines(molecules_smi_path):
        smi, name = parts[0], parts[1]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        assert mol.HasSubstructMatch(OXIME_PAT), (
            f"{name}: protonated oxime C=N-[O+] not found in {smi!r}"
        )


def test_molecules_smi_both_isomers(molecules_smi_path):
    """Each base mol_id must have both E and Z isomers."""
    names = {parts[1] for parts in _smi_lines(molecules_smi_path)}
    base_ids = {n.rsplit("_", 1)[0] for n in names}
    for base in sorted(base_ids):
        assert f"{base}_E" in names, f"{base}: E isomer missing from molecules.smi"
        assert f"{base}_Z" in names, f"{base}: Z isomer missing from molecules.smi"


def test_benchmark_meta_count(benchmark_meta):
    assert len(benchmark_meta) == EXPECTED_MOL_COUNT, (
        f"Expected {EXPECTED_MOL_COUNT} entries in benchmark_meta.json, got {len(benchmark_meta)}"
    )


def test_benchmark_meta_keys(benchmark_meta):
    required = {"smiles", "pct_A", "pct_B", "exp_outcome"}
    for mol_id, entry in benchmark_meta.items():
        missing = required - entry.keys()
        assert not missing, f"{mol_id}: missing keys {missing}"


def test_benchmark_meta_outcome_values(benchmark_meta):
    for mol_id, entry in benchmark_meta.items():
        assert entry["exp_outcome"] in ("R", "F"), (
            f"{mol_id}: exp_outcome={entry['exp_outcome']!r}, expected 'R' or 'F'"
        )


def test_benchmark_meta_ids_match_smi(benchmark_meta, molecules_smi_path):
    """mol_ids in benchmark_meta should appear as base names in molecules.smi."""
    smi_bases = {
        parts[1].rsplit("_", 1)[0]
        for parts in _smi_lines(molecules_smi_path)
    }
    for mol_id in benchmark_meta:
        assert mol_id in smi_bases, (
            f"{mol_id} in benchmark_meta.json has no corresponding entry in molecules.smi"
        )
