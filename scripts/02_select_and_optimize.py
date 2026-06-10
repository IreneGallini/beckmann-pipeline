"""
Step 2: select lowest-energy conformer per molecule, run AIMNet2 optimization
Input: Reads result SDF from 01_smiles_to_conformers.py
Output: one geometry per molecule in data/output/aimnet_optimized/best_aimnet_optimized.sdf
"""
import os
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE" 

from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from aimnet.calculators import AIMNet2Calculator, AIMNet2ASE
from ase.optimize import LBFGS
from ase import Atoms




PROJECT_ROOT = Path(__file__).parent.parent 

conformer_dir  = PROJECT_ROOT / "data" / "output"  / "conformers"
output_dir = PROJECT_ROOT / "data" / "output" / "aimnet_optimized"
output_dir.mkdir(parents=True, exist_ok=True)

# Find most recent Auto3D output .sdf
sdf_files = sorted(conformers_dir.glob("molecules_*/molecules_out.sdf"))
if not sdf_files:
    raise FileNotFoundError("No Auto3D output SDF found. Run 01_smiles_to_conformers.py first.")

sdf_path = sdf_files[-1]
print(f"Reading conformers from: {sdf_path}")

# Load conformers, group by molecule name
suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
molecules = {}  # name -> list of (energy, mol)

for mol in suppl:
    if mol is None:
        continue
    name   = mol.GetProp("_Name") if mol.HasProp("_Name") else "unknown"
    energy = float(mol.GetProp("E_tot")) if mol.HasProp("E_tot") else float("inf")
    molecules.setdefault(name, []).append((energy, mol))


# Pick lowest energy conformer per molecule
best_conformers = {}
for name, confs in molecules.items():
    confs.sort(key=lambda x: x[0])  # sort by energy ascending
    best_energy, best_mol = confs[0]
    best_conformers[name] = best_mol
    print(f"{name}: selected conformer with E = {best_energy:.6f} Hartree "
          f"({best_energy * 627.509:.2f} kcal/mol) from {len(confs)} candidates")

# AIMNet2 optimization
calc = aimnet2calc.AIMNet2Calculator("aimnet2_b973c_0.jpt")

writer = Chem.SDWriter(str(output_dir / "best_aimnet_optimized.sdf"))

for name, mol in best_conformers.items():
    print(f"  Optimizing {name}...")
    try:
        opt_mol = calc.optimize(mol, fmax=0.05)  # convergence threshold in eV/Å
        energy  = calc.get_energy(opt_mol)
        opt_mol.SetProp("_Name", name)
        opt_mol.SetProp("E_aimnet_eV",    f"{energy:.6f}")
        opt_mol.SetProp("E_aimnet_kcal",  f"{energy * 23.0605:.4f}")
        writer.write(opt_mol)
        print(f"    Done. AIMNet2 E = {energy:.6f} eV ({energy * 23.0605:.2f} kcal/mol)")
    except Exception as e:
        print(f"    WARNING: AIMNet2 optimization failed for {name}: {e}")

writer.close()
print(f"\nStep 2 done. Optimized structures: {output_dir / 'best_aimnet_optimized.sdf'}")


'''
if __name__ == "__main__":
    path = os.path.join(root, "example/files/smiles.smi")  # You can specify the path to your file here
    config = Auto3DOptions(path=path, k=1, use_gpu=False)  # Configure Auto3D parameters
    out = main(config)  # Run Auto3D and get output path
    print(out)
'''

