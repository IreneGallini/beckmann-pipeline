"""
Step 0: Convert ketone SMILES (chemdraw.txt) → oxime SMILES (molecules.smi)

For each ketone:
  1. Detect the C=O group (excludes esters, aldehydes, amides)
  2. Convert C=O → C=N-OH using a reaction SMARTS
  3. Enumerate RDKit tautomers; keep only true oxime forms (C=N-OH)
  4. Generate both E and Z isomers of the C=N bond
  5. Write SMILES + name to data/input/molecules.smi

Output feeds directly into script 01 (Auto3D conformer generation).

Atom tracking note: across tautomers the oxime N and O may land at
different atom indices. For downstream DFT/NBO7 work, stamp atom map
numbers onto the oxime N (:3) and O (:4) here before writing.
"""

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers, StereoEnumerationOptions,
)

# Ketone: sp3/sp2 C with =O and exactly two carbon neighbours (no O, N, S attached)
KETONE_PAT = Chem.MolFromSmarts('[C;$([C](=O)([c,C])[c,C])]')
# True oxime: C double-bonded to N, N single-bonded to OH
OXIME_PAT = Chem.MolFromSmarts('[C]=[N]-[OH1]')
# Reaction: replace the ketone oxygen with N-OH
# Unmapped O in reactant is consumed; new N and O are added in product
RXN = rdChemReactions.ReactionFromSmarts(
    '[#6:2]-[C:1](=[O])-[#6:3]>>[#6:2]-[C:1](=NO)-[#6:3]'
)


def _fix_smiles(smiles: str) -> str:
    """Replace [O] (ChemDraw artifact for phenol O with missing H) with O."""
    return smiles.replace('[O]', 'O')


def parse_chemdraw(path: Path) -> list[tuple[str, str]]:
    """Parse alternating-line format: ID, SMILES, blank. Return (mol_id, smiles) list."""
    lines = [l.strip() for l in path.read_text().splitlines() if l.strip()]
    return [(lines[i], lines[i + 1]) for i in range(0, len(lines) - 1, 2)]


def mol_from_smiles(smiles: str, mol_id: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        fixed = _fix_smiles(smiles)
        mol = Chem.MolFromSmiles(fixed)
        if mol is not None:
            print(f"  [WARN] {mol_id}: [O] replaced with O — using fixed SMILES")
        else:
            print(f"  [ERROR] {mol_id}: unparseable SMILES, skipping")
    return mol


def ketone_to_oximes(mol: Chem.Mol) -> list[Chem.Mol]:
    """Run C=O → C=N-OH reaction; return unique sanitised products."""
    seen: set[str] = set()
    result: list[Chem.Mol] = []
    for prod_tuple in RXN.RunReactants((mol,)):
        p = prod_tuple[0]
        try:
            Chem.SanitizeMol(p)
            smi = Chem.MolToSmiles(p)
            if smi not in seen:
                seen.add(smi)
                result.append(p)
        except Exception:
            pass
    return result


def enumerate_tautomers(mol: Chem.Mol) -> list[Chem.Mol]:
    """Return tautomers that still contain the C=N-OH (oxime) motif."""
    tautomers = rdMolStandardize.TautomerEnumerator().Enumerate(mol)
    valid = [t for t in tautomers if t.HasSubstructMatch(OXIME_PAT)]
    return valid if valid else [mol]


def get_cn_stereo(mol: Chem.Mol) -> str:
    """Return 'E', 'Z', or 'noEZ' for the C=N double bond in mol."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for bond in mol.GetBonds():
        if bond.GetBondType() != rdchem.BondType.DOUBLE:
            continue
        syms = {bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()}
        if syms == {'C', 'N'}:
            stereo = bond.GetStereo()
            if stereo == rdchem.BondStereo.STEREOE:
                return 'E'
            if stereo == rdchem.BondStereo.STEREOZ:
                return 'Z'
    return 'noEZ'


def enumerate_ez(mol: Chem.Mol) -> list[tuple[Chem.Mol, str]]:
    """Enumerate E/Z isomers of unassigned stereocenters; label each by C=N stereo."""
    opts = StereoEnumerationOptions(unique=True, onlyUnassigned=True)
    isomers = list(EnumerateStereoisomers(mol, options=opts))
    if not isomers:
        isomers = [mol]
    return [(iso, get_cn_stereo(iso)) for iso in isomers]


def main() -> None:
    root = Path(__file__).parent.parent
    chemdraw_path = root / 'data' / 'input' / 'chemdraw.txt'
    output_path   = root / 'data' / 'input' / 'molecules.smi'

    pairs = parse_chemdraw(chemdraw_path)

    seen_smiles: set[str] = set()
    entries: list[tuple[str, str]] = []

    print(f"\n{'mol_id':<14} {'ketones':>7} {'tautomers':>9} {'E/Z out':>8}")
    print('-' * 42)

    for mol_id, smiles in pairs:
        mol = mol_from_smiles(smiles, mol_id)
        if mol is None:
            print(f"  {mol_id:<14} {'—':>7} {'—':>9} {'—':>8}")
            continue

        if not mol.HasSubstructMatch(KETONE_PAT):
            print(f"  [WARN] {mol_id}: no ketone C=O found, skipping")
            print(f"  {mol_id:<14} {'0':>7} {'—':>9} {'—':>8}")
            continue

        oximes = ketone_to_oximes(mol)
        n_ketones = len(oximes)
        mol_entries: list[tuple[str, str]] = []

        for ox_idx, oxime in enumerate(oximes, 1):
            tautomers = enumerate_tautomers(oxime)
            for tau_idx, tau in enumerate(tautomers, 1):
                for iso, ez_label in enumerate_ez(tau):
                    smi = Chem.MolToSmiles(iso)
                    # Canonicalise for deduplication across molecules
                    canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
                    if canon in seen_smiles:
                        continue
                    seen_smiles.add(canon)
                    ox_part = f"_ox{ox_idx}" if n_ketones > 1 else ""
                    name = f"{mol_id}{ox_part}_tau{tau_idx}_{ez_label}"
                    mol_entries.append((smi, name))

        n_tau = len({n.split('_tau')[1].split('_')[0] for _, n in mol_entries}) if mol_entries else 0
        print(f"  {mol_id:<14} {n_ketones:>7} {n_tau:>9} {len(mol_entries):>8}")
        entries.extend(mol_entries)

    with open(output_path, 'w') as f:
        for smi, name in entries:
            f.write(f"{smi} {name}\n")

    print(f"\nTotal: {len(entries)} entries → {output_path}")


if __name__ == '__main__':
    main()