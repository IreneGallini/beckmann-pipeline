"""
Convert the 34-molecule benchmark CSV to molecules.smi + benchmark_meta.json.
Moved here (from beckmann/oximes.py's benchmark_to_oximes()/main()) since
it's CSV-batch-specific glue on top of beckmann_core.oximes's pure
per-molecule functions, not part of either product's own API.
"""
import csv
import json

from rdkit import Chem

from beckmann_core.oximes import KETONE_PAT, enumerate_ez, ketone_to_protonated_oximes
from beckmann_nbo.config import DATA_INPUT


def benchmark_to_oximes(csv_path, out_smi, out_meta) -> list[tuple[str, str]]:
    """Convert benchmark CSV to molecules.smi + benchmark_meta.json.

    Returns list of (smiles, name) pairs written to out_smi.
    """
    meta: dict = {}
    entries: list[tuple[str, str]] = []

    print(f"\n{'mol_id':<14}  {'E/Z out':>7}  SMILES")
    print('-' * 72)

    with open(csv_path, newline='') as f:
        for row in csv.DictReader(f):
            mol_id = f"mol_{int(row['id']):03d}"
            smiles = row['SMILES']
            meta[mol_id] = {
                'smiles':      smiles,
                'pct_A':       row['% product A'],
                'pct_B':       row['% product B'],
                'exp_outcome': row['exp_outcome'],
            }

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"  [ERROR] {mol_id}: invalid SMILES")
                continue

            if not mol.HasSubstructMatch(KETONE_PAT):
                print(f"  [WARN]  {mol_id}: no ketone found, skipping")
                continue

            oximes = ketone_to_protonated_oximes(mol)
            if not oximes:
                print(f"  [WARN]  {mol_id}: reaction gave no products")
                continue

            mol_entries: list[tuple[str, str]] = []
            for ox in oximes:
                for iso, ez_label in enumerate_ez(ox):
                    smi_out = Chem.MolToSmiles(iso)
                    name    = f"{mol_id}_{ez_label}"
                    mol_entries.append((smi_out, name))

            for smi_out, name in mol_entries:
                print(f"  {name:<14}  {len(mol_entries):>7}  {smi_out}")
            entries.extend(mol_entries)

    with open(out_smi, 'w') as f:
        for smi, name in entries:
            f.write(f"{smi} {name}\n")

    with open(out_meta, 'w') as f:
        json.dump(meta, f, indent=2)

    print(f"\nTotal: {len(entries)} entries -> {out_smi}")
    print(f"Benchmark metadata  -> {out_meta}")
    return entries


def main() -> None:
    benchmark_to_oximes(
        csv_path = DATA_INPUT / 'benchmark.csv',
        out_smi  = DATA_INPUT / 'molecules.smi',
        out_meta = DATA_INPUT / 'benchmark_meta.json',
    )


if __name__ == '__main__':
    main()
