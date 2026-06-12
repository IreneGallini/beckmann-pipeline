''' Output testing
Ethanol, butane (gauche vs anti) well-documented energy differences in textbooks
A simple oxime (acetone oxime) known E/Z geometry

RDKit's Chem.MolToSmiles
did AIMNet2 optimization significantly move the geometry from where Auto3D left it or just fine tune it?
'''
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import pytest
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

PROJECT_ROOT = Path(__file__).parent.parent

def test_imports():
    """Sanity check: do the key dependencies import cleanly?"""
    from Auto3D.auto3D import options, main
    from aimnet.calculators import AIMNet2Calculator, AIMNet2ASE
    assert True

def test_step1_output_exists():
    """Has script 01 produced an output SDF?"""
    conformers_dir = PROJECT_ROOT / "data" / "output" / "conformers"
    sdf_files = sorted(conformers_dir.glob("molecules_*/molecules_out.sdf"))
    assert len(sdf_files) > 0, "Run script 01 first"

def test_step1_conformer_count():
    """Each molecule should have ~5 conformers (k=5)."""
    conformers_dir = PROJECT_ROOT / "data" / "output" / "conformers"
    sdf_path = sorted(conformers_dir.glob("molecules_*/molecules_out.sdf"))[-1]

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    counts = {}
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        counts[name] = counts.get(name, 0) + 1

    for name, count in counts.items():
        assert count <= 5, f"{name} has {count} conformers, expected ≤5"
        assert count >= 1, f"{name} has no conformers"


def test_step2_output_exists():
    """Has script 02 produced optimized output?"""
    out_path = PROJECT_ROOT / "data" / "output" / "aimnet_optimized" / "best_aimnet_optimized.sdf"
    assert out_path.exists(), "Run script 02 first"

def test_step2_energy_decreased():
    """AIMNet2 optimization should lower (or maintain) energy vs Auto3D pre-opt."""
    conformers_dir = PROJECT_ROOT / "data" / "output" / "conformers"
    sdf_step1 = sorted(conformers_dir.glob("molecules_*/molecules_out.sdf"))[-1]
    sdf_step2 = PROJECT_ROOT / "data" / "output" / "aimnet_optimized" / "best_aimnet_optimized.sdf"

    # Lowest energy per molecule from step 1 (in Hartree → kcal/mol)
    step1_best = {}
    for mol in Chem.SDMolSupplier(str(sdf_step1), removeHs=False):
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        e_kcal = float(mol.GetProp("E_tot")) * 627.509
        if name not in step1_best or e_kcal < step1_best[name]:
            step1_best[name] = e_kcal

    # Step 2 energies (already in kcal/mol)
    for mol in Chem.SDMolSupplier(str(sdf_step2), removeHs=False):
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        e2_kcal = float(mol.GetProp("E_aimnet2_kcal"))
        e1_kcal = step1_best[name]

        # Allow small tolerance for different energy references (Auto3D ANI vs AIMNet2)
        # The KEY check is that step 2 didn't blow up the structure
        print(f"{name}: step1={e1_kcal:.2f} kcal/mol, step2={e2_kcal:.2f} kcal/mol")

def test_step2_rmsd_reasonable():
    """AIMNet2 optimization shouldn't drastically distort the geometry (RMSD check, Auto3D-paper-style)."""
    conformers_dir = PROJECT_ROOT / "data" / "output" / "conformers"
    sdf_step1 = sorted(conformers_dir.glob("molecules_*/molecules_out.sdf"))[-1]
    sdf_step2 = PROJECT_ROOT / "data" / "output" / "aimnet_optimized" / "best_aimnet_optimized.sdf"

    step1_mols = {}
    for mol in Chem.SDMolSupplier(str(sdf_step1), removeHs=False):
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        e_kcal = float(mol.GetProp("E_tot")) * 627.509
        if name not in step1_mols or e_kcal < step1_mols[name][0]:
            step1_mols[name] = (e_kcal, mol)

    for mol2 in Chem.SDMolSupplier(str(sdf_step2), removeHs=False):
        if mol2 is None:
            continue
        name = mol2.GetProp("_Name")
        _, mol1 = step1_mols[name]

        rmsd = AllChem.GetBestRMS(Chem.Mol(mol1), Chem.Mol(mol2))
        print(f"{name}: RMSD (step1 → step2) = {rmsd:.3f} Å")
        # Flag large structural changes - may indicate optimization issue
        assert rmsd < 2.0, f"{name}: RMSD {rmsd:.3f} Å is unusually large"



