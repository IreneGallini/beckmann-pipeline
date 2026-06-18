"""
Read benchmark.csv, convert each ketone SMILES to protonated activated oxime
(C=N-[OH2+]), generate E and Z isomers, write to data/input/molecules.smi.
Also writes data/input/benchmark_meta.json for downstream comparison in script 04.

Two-step conversion:
  1. C=O  →  C=N-OH   (neutral oxime, via reaction SMARTS)
  2. N-OH  →  N-[OH2+] (set formal charge +1 on O; sanitization gives 2 implicit H)
"""
import csv
import json
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdchem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers, StereoEnumerationOptions,
)

KETONE_PAT  = Chem.MolFromSmarts('[C;$([C](=O)([c,C])[c,C])]')
NEUTRAL_RXN = rdChemReactions.ReactionFromSmarts(
    '[#6:2]-[C:1](=[O])-[#6:3]>>[#6:2]-[C:1](=NO)-[#6:3]'
)
OXIME_O_PAT = Chem.MolFromSmarts('[C]=[N]-[OH1]')


def get_cn_stereo(mol: Chem.Mol) -> str:
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for bond in mol.GetBonds():
        if bond.GetBondType() != rdchem.BondType.DOUBLE:
            continue
        if {bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()} == {'C', 'N'}:
            stereo = bond.GetStereo()
            if stereo == rdchem.BondStereo.STEREOE:
                return 'E'
            if stereo == rdchem.BondStereo.STEREOZ:
                return 'Z'
    return 'noEZ'


def ketone_to_protonated_oximes(mol: Chem.Mol) -> list[Chem.Mol]:
    """Two-step: C=O → C=N-OH → C=N-[OH2+]."""
    seen: set[str] = set()
    result: list[Chem.Mol] = []

    for prod_tuple in NEUTRAL_RXN.RunReactants((mol,)):
        p = prod_tuple[0]
        try:
            Chem.SanitizeMol(p)
        except Exception:
            continue

        match = p.GetSubstructMatch(OXIME_O_PAT)
        if not match:
            continue

        o_idx = match[2]
        rw = Chem.RWMol(p)
        rw.GetAtomWithIdx(o_idx).SetFormalCharge(1)
        try:
            Chem.SanitizeMol(rw)
            protonated = rw.GetMol()
            smi = Chem.MolToSmiles(protonated)
            if smi not in seen:
                seen.add(smi)
                result.append(protonated)
        except Exception:
            pass

    return result


def enumerate_ez(mol: Chem.Mol) -> list[tuple[Chem.Mol, str]]:
    opts = StereoEnumerationOptions(unique=True, onlyUnassigned=True)
    isomers = list(EnumerateStereoisomers(mol, options=opts))
    if not isomers:
        isomers = [mol]
    return [(iso, get_cn_stereo(iso)) for iso in isomers]


def main() -> None:
    root      = Path(__file__).parent.parent
    csv_path  = root / 'data' / 'input' / 'benchmark.csv'
    smi_path  = root / 'data' / 'input' / 'molecules.smi'
    meta_path = root / 'data' / 'input' / 'benchmark_meta.json'

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

    with open(smi_path, 'w') as f:
        for smi, name in entries:
            f.write(f"{smi} {name}\n")

    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)

    print(f"\nTotal: {len(entries)} entries → {smi_path}")
    print(f"Benchmark metadata  → {meta_path}")


if __name__ == '__main__':
    main()
