"""Stage 0a: SMILES -> protonated oxime isomers (E/Z) -> 3D conformers (Auto3D).

Edit MOL_NAME/SMILES below for a new molecule, then:
    python 00_smiles_to_conformers.py
"""
import sys
from pathlib import Path

from rdkit import Chem

from _common import QUERY_PREFIX, sanitize_id, workdir_for

from beckmann_core.oximes import KETONE_PAT, enumerate_ez, ketone_to_protonated_oximes
from beckmann_core.conformers import generate_conformers

MOL_NAME = "test1"
SMILES = "O=C1CCC2=C1C=CC=C2"


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    workdir = workdir_for(mol_id)
    workdir.mkdir(parents=True, exist_ok=True)

    mol = Chem.MolFromSmiles(SMILES)
    if mol is None:
        print(f"ERROR: could not parse SMILES: {SMILES!r}", file=sys.stderr)
        sys.exit(1)
    if not mol.HasSubstructMatch(KETONE_PAT):
        print(f"ERROR: no ketone group found in: {SMILES!r}", file=sys.stderr)
        sys.exit(1)

    oximes = ketone_to_protonated_oximes(mol)
    if not oximes:
        print(f"ERROR: ketone-to-oxime conversion produced no products for: {SMILES!r}", file=sys.stderr)
        sys.exit(1)

    smi_path = workdir / f"{QUERY_PREFIX}_{mol_id}.smi"
    with open(smi_path, "w") as f:
        for ox in oximes:
            for iso, ez_label in enumerate_ez(ox):
                f.write(f"{Chem.MolToSmiles(iso)} {QUERY_PREFIX}_{mol_id}_{ez_label}\n")
    print(f"Wrote {smi_path}")

    print(f"\nGenerating conformers for {mol_id}...")
    sdf_path = generate_conformers(smi_path, workdir / "conformers")
    print(f"Conformers written to: {sdf_path}")
    print(f"\nNext: edit MOL_NAME/SMILES in 01_optimize_aimnet2.py to '{MOL_NAME}' and run it.")


if __name__ == "__main__":
    main()
