"""
Step 2: select lowest-energy conformer per molecule, run AIMNet2 optimization
Input: Reads result SDF from 01_smiles_to_conformers.py
Output: one geometry per molecule in data/output/aimnet_optimized/best_aimnet_optimized.sdf
"""
import os
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import torch
import aimnet2calc

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE" 


PROJECT_ROOT = Path(__file__).parent.parent 

conformer_dir  = PROJECT_ROOT / "data" / "output"  / "conformers"
output_dir = PROJECT_ROOT / "data" / "output" / "aimnet_optimized"
output_dir.mkdir(parents=True, exist_ok=True)