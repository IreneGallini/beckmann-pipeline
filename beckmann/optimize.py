"""
Step 2: select lowest-energy conformer per molecule, run AIMNet2/ASE optimization.

Writes two output files:
  best_aimnet_optimized.sdf  — one optimized structure per isomer name (mol_XXX_E / mol_XXX_Z)
  best_per_substrate.sdf     — one structure per substrate (lowest-energy isomer wins)

Downstream DFT scripts read from best_per_substrate.sdf.
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

from beckmann.config import DATA_OUTPUT


def select_and_optimize(
    conformers_sdf: Path,
    output_dir: Path,
) -> tuple[Path, Path]:
    """Pick best conformer per isomer, run AIMNet2/ASE optimization.

    Returns (best_aimnet_optimized.sdf, best_per_substrate.sdf).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading conformers from: {conformers_sdf}")
    suppl = Chem.SDMolSupplier(str(conformers_sdf), removeHs=False)
    molecules: dict[str, list[tuple[float, Chem.Mol]]] = {}

    for mol in suppl:
        if mol is None:
            continue
        name   = mol.GetProp("_Name") if mol.HasProp("_Name") else "unknown"
        energy = float(mol.GetProp("E_tot")) if mol.HasProp("E_tot") else float("inf")
        molecules.setdefault(name, []).append((energy, mol))

    print("Selecting lowest energy conformer per molecule:")
    best_conformers: dict[str, Chem.Mol] = {}
    for name, confs in molecules.items():
        confs.sort(key=lambda x: x[0])
        best_energy, best_mol = confs[0]
        best_conformers[name] = best_mol
        print(f"  {name}: selected E = {best_energy:.6f} Hartree "
              f"({best_energy * 627.509:.2f} kcal/mol) from {len(confs)} candidates")

    print(f"\nRunning AIMNet2 optimization on {len(best_conformers)} structure(s)...")
    base_calc = AIMNet2Calculator("aimnet2_2025")

    best_sdf_path = output_dir / "best_aimnet_optimized.sdf"
    writer = Chem.SDWriter(str(best_sdf_path))
    optimized: dict[str, tuple[float, Chem.Mol]] = {}

    for name, mol in best_conformers.items():
        print(f"  Optimizing {name}...")
        try:
            conf    = mol.GetConformer()
            coords  = conf.GetPositions()
            numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
            atoms   = Atoms(numbers=numbers, positions=coords)

            mol_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
            atoms.calc = AIMNet2ASE(base_calc, charge=mol_charge)

            opt = LBFGS(atoms, logfile=None)
            opt.run(fmax=0.05)

            energy_ev   = atoms.get_potential_energy()
            energy_kcal = energy_ev * 23.0605

            new_conf = mol.GetConformer()
            for i, pos in enumerate(atoms.get_positions()):
                new_conf.SetAtomPosition(i, pos.tolist())

            mol.SetProp("_Name",          name)
            mol.SetProp("E_aimnet2_eV",   f"{energy_ev:.6f}")
            mol.SetProp("E_aimnet2_kcal", f"{energy_kcal:.4f}")
            writer.write(mol)
            optimized[name] = (energy_ev, mol)
            print(f"    AIMNet2 E = {energy_ev:.6f} eV ({energy_kcal:.2f} kcal/mol)")

        except Exception as e:
            print(f"    WARNING: AIMNet2 optimization failed for {name}: {e}")

    writer.close()

    substrate_best: dict[str, tuple[float, Chem.Mol]] = {}
    for name, (energy_ev, mol) in optimized.items():
        base = name.rsplit("_", 1)[0]
        if base not in substrate_best or energy_ev < substrate_best[base][0]:
            substrate_best[base] = (energy_ev, mol)

    sub_sdf_path = output_dir / "best_per_substrate.sdf"
    writer2 = Chem.SDWriter(str(sub_sdf_path))
    for base in sorted(substrate_best):
        _, mol = substrate_best[base]
        writer2.write(mol)
    writer2.close()

    print(f"\nDone. Optimized structures saved to:")
    print(f"  {best_sdf_path}  ({len(optimized)} isomers)")
    print(f"  {sub_sdf_path}  ({len(substrate_best)} substrates, lowest-energy isomer only)")
    return best_sdf_path, sub_sdf_path


def main() -> None:
    conformers_dir = DATA_OUTPUT / "conformers"
    output_dir     = DATA_OUTPUT / "aimnet_optimized"

    sdf_files = sorted(conformers_dir.glob("molecules_*/molecules_out.sdf"))
    if not sdf_files:
        raise FileNotFoundError(
            "No Auto3D output SDF found. Run beckmann.conformers.main() first."
        )
    select_and_optimize(sdf_files[-1], output_dir)


if __name__ == '__main__':
    main()
