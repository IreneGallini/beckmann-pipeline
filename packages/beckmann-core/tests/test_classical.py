"""Unit tests for beckmann_core.classical -- atom-mapping + the classical
anti-periplanar baseline rule."""
from rdkit import Chem

from beckmann_core.classical import ANTI_THRESH, get_oxime_atoms, predict


def test_get_oxime_atoms_identifies_aryl_and_allyl():
    # alpha-tetralone-derived oxime: aromatic ring carbon vs. aliphatic ring carbon
    from beckmann_core.oximes import ketone_to_protonated_oximes
    mol = Chem.MolFromSmiles("O=C1CCC2=CC=CC=C21")
    (ox,) = ketone_to_protonated_oximes(mol)
    atom_ids = get_oxime_atoms(ox)
    assert atom_ids is not None
    cox_idx, nox_idx, oox_idx, c_aryl_idx, c_allyl_idx = atom_ids
    assert ox.GetAtomWithIdx(c_aryl_idx).GetIsAromatic()
    assert not ox.GetAtomWithIdx(c_allyl_idx).GetIsAromatic()


def test_get_oxime_atoms_no_match_returns_none():
    mol = Chem.MolFromSmiles("CCO")
    assert get_oxime_atoms(mol) is None


def test_predict_anti_aryl_gives_rearrangement():
    assert predict(d_aryl=170.0, d_allyl=20.0) == "R"


def test_predict_anti_allyl_gives_fragmentation():
    assert predict(d_aryl=20.0, d_allyl=170.0) == "F"


def test_predict_neither_anti_gives_inspect():
    assert predict(d_aryl=90.0, d_allyl=90.0) == "inspect"


def test_predict_both_anti_picks_larger():
    assert predict(d_aryl=ANTI_THRESH, d_allyl=179.0) == "F"
    assert predict(d_aryl=179.0, d_allyl=ANTI_THRESH) == "R"
