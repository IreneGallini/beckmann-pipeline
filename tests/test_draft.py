''' Output testing
Ethanol, butane (gauche vs anti) well-documented energy differences in textbooks
A simple oxime (acetone oxime) known E/Z geometry

RDKit's Chem.MolToSmiles
did AIMNet2 optimization significantly move the geometry from where Auto3D left it or just fine tune it?
'''

from aimnet.calculators import AIMNet2Calculator, AIMNet2ASE
print("AIMNet2Calculator OK")

# RMSD testing built into RDKit
from rdkit.Chem import AllChem
rmsd = AllChem.GetBestRMS(mol_aimnet_optimized, mol_auto3d_conformer)

# compare energies -> check if energy decreased after optimization 



