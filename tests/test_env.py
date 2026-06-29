"""Environment sanity checks — run these first to verify the conda env is correct."""


def test_rdkit_import():
    from rdkit import Chem
    assert Chem.MolFromSmiles("c1ccccc1") is not None


def test_auto3d_import():
    from Auto3D.auto3D import options, main


def test_aimnet_import():
    from aimnet.calculators import AIMNet2Calculator, AIMNet2ASE


def test_ase_import():
    from ase.optimize import LBFGS
    from ase import Atoms
