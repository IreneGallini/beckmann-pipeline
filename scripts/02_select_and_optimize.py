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

conformers_dir  = PROJECT_ROOT / "data" / "output"  / "conformers"
output_dir = PROJECT_ROOT / "data" / "output" / "aimnet_optimized"
output_dir.mkdir(parents=True, exist_ok=True)

# Find most recent Auto3D output .sdfx
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
print("Selecting lowest energy conformer per molecule:")
best_conformers = {}
for name, confs in molecules.items():
    confs.sort(key=lambda x: x[0])  # sort by energy ascending
    best_energy, best_mol = confs[0]
    best_conformers[name] = best_mol
    print(f"{name}: selected conformer with E = {best_energy:.6f} Hartree "
          f"({best_energy * 627.509:.2f} kcal/mol) from {len(confs)} candidates")

# AIMNet2 optimization
print(f"\nRunning AIMNet2 optimization on {len(best_conformers)} structure(s)...")

base_calc = AIMNet2Calculator("aimnet2_2025")
writer = Chem.SDWriter(str(output_dir / "best_aimnet_optimized.sdf"))

for name, mol in best_conformers.items():
    print(f"  Optimizing {name}...")
    try:
        # Convert RDKit mol to ASE Atoms
        conf    = mol.GetConformer()
        coords  = conf.GetPositions()
        numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        atoms   = Atoms(numbers=numbers, positions=coords)

        # Attach AIMNet2 calculator; charge is auto-detected from formal charges
        # (0 for neutral oximes, +1 for protonated C=N-[OH2+] benchmark structures)
        mol_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        atoms.calc = AIMNet2ASE(base_calc, charge=mol_charge)

        # Run optimization
        opt = LBFGS(atoms, logfile=None)
        opt.run(fmax=0.05)  # convergence: max force < 0.05 eV/Å

        # Get final energy
        energy_ev   = atoms.get_potential_energy()
        energy_kcal = energy_ev * 23.0605

        # Write back optimized coords to RDKit mol
        new_conf = mol.GetConformer()
        for i, pos in enumerate(atoms.get_positions()):
            new_conf.SetAtomPosition(i, pos.tolist())

        mol.SetProp("_Name",          name)
        mol.SetProp("E_aimnet2_eV",   f"{energy_ev:.6f}")
        mol.SetProp("E_aimnet2_kcal", f"{energy_kcal:.4f}")
        writer.write(mol)

        print(f"    AIMNet2 E = {energy_ev:.6f} eV ({energy_kcal:.2f} kcal/mol)")
    
        
    except Exception as e:
        print(f"    WARNING: AIMNet2 optimization failed for {name}: {e}")

writer.close()
print(f"\nDone. Optimized structures saved to:")
print(f"  {output_dir / 'best_aimnet_optimized.sdf'}")


