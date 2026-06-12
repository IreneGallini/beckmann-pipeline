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

# RMSD testing built into RDKit
from rdkit.Chem import AllChem
rmsd = AllChem.GetBestRMS(mol_aimnet_optimized, mol_auto3d_conformer)

# compare energies -> check if energy decreased after optimization 



