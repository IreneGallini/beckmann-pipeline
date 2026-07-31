"""
Unit tests for beckmann_core.oximes -- pure Mol -> Mol functions, no file I/O.
(The original tests/test_step0_oximes.py validated the benchmark-CSV batch
output of scripts/00_benchmark_to_oximes.py; that test and the script it
covers both now live in research/benchmark_pipeline/, since neither is part
of beckmann-core's own API.)
"""
from rdkit import Chem

from beckmann_core.oximes import enumerate_ez, get_cn_stereo, ketone_to_protonated_oximes

OXIME_CATION_PAT = Chem.MolFromSmarts("[C]=[N]-[O+]")

ALPHA_TETRALONE = "O=C1CCC2=CC=CC=C21"  # symmetric ring ketone -- yields E/Z isomers
ACETONE = "CC(C)=O"  # no aryl/allyl asymmetry, but still a valid ketone


def test_ketone_to_protonated_oximes_produces_cationic_oxime():
    mol = Chem.MolFromSmiles(ALPHA_TETRALONE)
    oximes = ketone_to_protonated_oximes(mol)
    assert oximes, "expected at least one oxime product"
    for ox in oximes:
        assert ox.HasSubstructMatch(OXIME_CATION_PAT)


def test_ketone_to_protonated_oximes_no_ketone_returns_empty():
    mol = Chem.MolFromSmiles("CCO")  # ethanol, no ketone
    assert ketone_to_protonated_oximes(mol) == []


def test_enumerate_ez_labels_every_isomer():
    mol = Chem.MolFromSmiles(ALPHA_TETRALONE)
    (ox,) = ketone_to_protonated_oximes(mol)
    isomers = enumerate_ez(ox)
    assert isomers, "expected at least one E/Z isomer"
    labels = {label for _, label in isomers}
    assert labels <= {"E", "Z", "noEZ"}
    for iso, label in isomers:
        assert get_cn_stereo(iso) == label


def test_enumerate_ez_deduplicates():
    """acetone's C=N has no distinguishable E/Z (both substituents are
    methyl) -- enumerate_ez should fall back to the single unmodified mol,
    not fabricate two identical isomers."""
    mol = Chem.MolFromSmiles(ACETONE)
    (ox,) = ketone_to_protonated_oximes(mol)
    isomers = enumerate_ez(ox)
    assert len(isomers) == 1
    assert isomers[0][1] == "noEZ"
