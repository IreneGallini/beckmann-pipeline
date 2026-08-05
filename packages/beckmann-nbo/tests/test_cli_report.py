"""Tests for cli_report.py's classical-vs-wCNmax comparison line -- built on
a synthetic, embedded RDKit Mol (not real data/output/ state) so this runs
in a fresh clone. Covers the "inspect" outcome explicitly, since
beckmann_core.classical.predict() is a genuine 3-way R/F/inspect result the
CLI must handle, not fold into R/F.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

import beckmann_nbo.cli_report as cli_report


def _embedded_oxime_mol(name: str) -> Chem.Mol:
    """A small protonated-oxime Mol (aryl ketone -> oxime, as
    beckmann_core.oximes.ketone_to_protonated_oximes produces) with a real
    3D conformer, so get_oxime_atoms()/GetDihedralDeg() have real geometry
    to work with -- exact R/F outcome isn't asserted (that's
    beckmann_core.classical's own test surface), only that the CLI-layer
    plumbing (atom lookup -> dihedral -> predict -> formatted line) runs
    end to end without crashing and reports one of the three valid
    outcomes."""
    mol = Chem.MolFromSmiles("C1=CC=C(C=C1)/C(=[NH+]\\O)/CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol.SetProp("_Name", name)
    return mol


def test_classical_vs_wcnmax_line_runs_end_to_end(monkeypatch):
    mol_name = "mol_test_E"
    rdkit_mol = _embedded_oxime_mol(mol_name)
    monkeypatch.setattr(cli_report, "_load_mols", lambda: {mol_name: rdkit_mol})

    line = cli_report.classical_vs_wcnmax_line(mol_name, c_map={}, wcnmax_pred="R")
    assert line.startswith(f"{mol_name}: classical=")
    assert "wcnmax=R" in line
    classical_outcome = line.split("classical=")[1].split()[0]
    assert classical_outcome in ("R", "F", "inspect")


def test_classical_vs_wcnmax_line_missing_molecule(monkeypatch):
    monkeypatch.setattr(cli_report, "_load_mols", lambda: {})

    line = cli_report.classical_vs_wcnmax_line("mol_missing_E", c_map={}, wcnmax_pred="F")
    assert "classical=unavailable" in line
    assert "wcnmax=F" in line


def test_classical_vs_wcnmax_line_no_oxime_atoms(monkeypatch):
    mol_name = "mol_no_oxime_E"
    plain_mol = Chem.MolFromSmiles("CCO")  # ethanol, no oxime at all
    monkeypatch.setattr(cli_report, "_load_mols", lambda: {mol_name: plain_mol})

    line = cli_report.classical_vs_wcnmax_line(mol_name, c_map={}, wcnmax_pred="R")
    assert "classical=inspect" in line
    assert "agreement=unclear" in line
