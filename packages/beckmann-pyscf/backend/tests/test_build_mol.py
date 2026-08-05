"""
Covers build_mol()'s per-element basis-patch dispatch (pyscf_livvo.py) --
see data/output/analysis/heteroatom_basis_coverage.csv for the full
heteroatom audit this guards. The property under test: vendoring a patch
for one element must never change basis construction for any other
element (additive-only, per the module's own docstring).
"""
import pytest
from pyscf import gto

from beckmann_pyscf.engine.pyscf_livvo import build_mol


def _case(atom_spec, charge=0, spin=0):
    return {"atom_spec": atom_spec, "charge": charge, "spin": spin}


def test_unpatched_element_keeps_plain_string_basis():
    mol = build_mol(_case([("C", (0, 0, 0)), ("H", (0, 0, 1.09))], spin=1))
    assert isinstance(mol.basis, str)


def test_br_builds_with_vendored_patch():
    mol = build_mol(_case([("Br", (0, 0, 0))], spin=1))
    assert mol.nao > 0


def test_as_builds_with_vendored_patch():
    mol = build_mol(_case([("As", (0, 0, 0))], spin=3))
    assert mol.nao > 0


def test_mixed_br_as_case_is_additive():
    br_only = build_mol(_case([("Br", (0, 0, 0))], spin=1))
    as_only = build_mol(_case([("As", (0, 0, 0))], spin=3))
    mixed = build_mol(_case([("Br", (0, 0, 0)), ("As", (0, 0, 3.0))], spin=0))
    assert mixed.nao == br_only.nao + as_only.nao


def test_iodine_has_no_vendored_ground_truth():
    with pytest.raises(gto.basis.BasisNotFoundError):
        build_mol(_case([("I", (0, 0, 0))], spin=1))
