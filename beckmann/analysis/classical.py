"""
Baseline benchmark: does the classical anti-periplanar rule predict Beckmann outcomes?

Classical rule: the substituent ANTI (|dihedral| ≥150°) to the N-O leaving group
migrates. Aryl anti → rearrangement (R); alkyl anti → fragmentation (F).

Result: ~59% accuracy on 34 molecules — motivates the DFT/NBO approach.
"""
import csv
import json
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from beckmann.config import DATA_INPUT, DATA_OUTPUT

# Uses neutral [O;!R] (not [O+]) because best_aimnet_optimized.sdf has all isomers
# including neutral-like geometries; the protonation state matters less for dihedral analysis.
OXIME_PAT   = Chem.MolFromSmarts('[C:1]=[N:2]-[O;!R:3]')
ANTI_THRESH = 150.0

FIELDS = [
    'mol_id',
    'c_aryl_idx', 'c_ox_idx', 'n_ox_idx', 'o_ox_idx', 'c_allyl_idx',
    'Emin_E_eV', 'Emin_Z_eV', 'delta_E_Z_E_kcal',
    'lowest_isomer',
    'dihedral_O_N_C_aryl', 'dihedral_O_N_C_allyl',
    'anti_group', 'beckmann_pred',
    'exp_pct_A', 'exp_pct_B', 'exp_outcome',
    'agreement', 'notes',
]


def abs_dihedral(conf, i: int, j: int, k: int, l: int) -> float:
    return abs(rdMolTransforms.GetDihedralDeg(conf, i, j, k, l))


def get_oxime_atoms(mol: Chem.Mol):
    match = mol.GetSubstructMatch(OXIME_PAT)
    if not match:
        return None
    cox_idx, nox_idx, oox_idx = match
    c_aryl_idx = c_allyl_idx = None
    for nbr in mol.GetAtomWithIdx(cox_idx).GetNeighbors():
        if nbr.GetIdx() == nox_idx:
            continue
        if nbr.GetIsAromatic():
            c_aryl_idx = nbr.GetIdx()
        else:
            c_allyl_idx = nbr.GetIdx()
    if c_aryl_idx is None or c_allyl_idx is None:
        return None
    return cox_idx, nox_idx, oox_idx, c_aryl_idx, c_allyl_idx


def predict(d_aryl: float, d_allyl: float) -> str:
    aryl_anti  = d_aryl  >= ANTI_THRESH
    allyl_anti = d_allyl >= ANTI_THRESH
    if aryl_anti and not allyl_anti:
        return 'R'
    if allyl_anti and not aryl_anti:
        return 'F'
    if aryl_anti and allyl_anti:
        return 'R' if d_aryl >= d_allyl else 'F'
    return 'inspect'


def mol_to_xyz(mol: Chem.Mol, title: str) -> str:
    conf  = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), title]
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<3}  {p.x:>12.6f}  {p.y:>12.6f}  {p.z:>12.6f}")
    return '\n'.join(lines) + '\n'


def run_benchmark(
    sdf_path: Path,
    meta_path: Path,
    out_dir: Path,
) -> list[dict]:
    """Run classical rule benchmark. Returns list of result row dicts."""
    csv_out  = out_dir / 'classical_rule_results.csv'
    xyz_root = out_dir / 'classical_rule_xyz'

    meta: dict = json.load(open(meta_path))
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    structures: dict[str, Chem.Mol] = {}
    for mol in suppl:
        if mol is None:
            continue
        structures[mol.GetProp('_Name')] = mol

    parents: dict[str, dict[str, str]] = {}
    for name in structures:
        parts  = name.rsplit('_', 1)
        parent = parts[0]
        ez     = parts[1] if len(parts) > 1 else 'noEZ'
        parents.setdefault(parent, {})[ez] = name

    rows = []
    for parent, ez_map in sorted(parents.items()):
        m   = meta.get(parent, {})
        row = {f: '' for f in FIELDS}
        row['mol_id']      = parent
        row['exp_pct_A']   = m.get('pct_A', '')
        row['exp_pct_B']   = m.get('pct_B', '')
        row['exp_outcome'] = m.get('exp_outcome', '')

        e_mols = {ez: structures[name] for ez, name in ez_map.items() if name in structures}
        energies = {
            ez: float(mol.GetProp('E_aimnet2_eV'))
            for ez, mol in e_mols.items()
            if mol.HasProp('E_aimnet2_eV')
        }

        if 'E' in energies:
            row['Emin_E_eV'] = f"{energies['E']:.6f}"
        if 'Z' in energies:
            row['Emin_Z_eV'] = f"{energies['Z']:.6f}"
        if 'E' in energies and 'Z' in energies:
            delta_kcal = (energies['Z'] - energies['E']) * 23.0605
            row['delta_E_Z_E_kcal'] = f"{delta_kcal:.2f}"

        if not energies:
            row['notes'] = 'no optimized structures'
            rows.append(row)
            continue

        lowest_ez = min(energies, key=energies.__getitem__)
        row['lowest_isomer'] = lowest_ez
        best_mol = e_mols[lowest_ez]
        conf     = best_mol.GetConformer()

        mol_dir = xyz_root / parent
        mol_dir.mkdir(parents=True, exist_ok=True)
        for ez, mol in e_mols.items():
            (mol_dir / f"{parent}_{ez}_opt.xyz").write_text(
                mol_to_xyz(mol, f"{parent}_{ez}_opt")
            )

        atom_ids = get_oxime_atoms(best_mol)
        if atom_ids is None:
            row['notes']         = 'oxime atoms not identified'
            row['beckmann_pred'] = 'inspect'
            row['agreement']     = 'unclear'
            rows.append(row)
            continue

        cox_idx, nox_idx, oox_idx, c_aryl_idx, c_allyl_idx = atom_ids
        row['c_ox_idx']    = cox_idx    + 1
        row['n_ox_idx']    = nox_idx    + 1
        row['o_ox_idx']    = oox_idx    + 1
        row['c_aryl_idx']  = c_aryl_idx + 1
        row['c_allyl_idx'] = c_allyl_idx + 1

        d_aryl  = abs_dihedral(conf, oox_idx, nox_idx, cox_idx, c_aryl_idx)
        d_allyl = abs_dihedral(conf, oox_idx, nox_idx, cox_idx, c_allyl_idx)
        row['dihedral_O_N_C_aryl']  = f"{d_aryl:.1f}"
        row['dihedral_O_N_C_allyl'] = f"{d_allyl:.1f}"

        pred = predict(d_aryl, d_allyl)
        row['beckmann_pred'] = pred
        row['anti_group']    = ('aryl' if pred == 'R'
                                else 'allyl' if pred == 'F'
                                else 'unclear')

        exp = row['exp_outcome']
        row['agreement'] = ('yes'     if pred == exp
                            else 'inspect' if pred == 'inspect'
                            else 'no')
        rows.append(row)
        print(f"  {parent}  lowest={lowest_ez}  "
              f"d_aryl={d_aryl:.1f}°  d_allyl={d_allyl:.1f}°  "
              f"pred={pred}  exp={exp}  → {row['agreement']}")

    out_dir.mkdir(parents=True, exist_ok=True)
    with open(csv_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    total   = len(rows)
    agree   = sum(1 for r in rows if r['agreement'] == 'yes')
    disagree= sum(1 for r in rows if r['agreement'] == 'no')
    inspect = sum(1 for r in rows if r['agreement'] in ('inspect', 'unclear'))
    print(f"\nResults  → {csv_out}")
    print(f"Summary: {total} molecules | agree {agree} | disagree {disagree} | inspect {inspect}")
    print(f"Classical rule accuracy: {agree}/{total} ({100*agree/total:.0f}%)")
    return rows


def main() -> None:
    run_benchmark(
        sdf_path  = DATA_OUTPUT / 'aimnet_optimized' / 'best_aimnet_optimized.sdf',
        meta_path = DATA_INPUT  / 'benchmark_meta.json',
        out_dir   = DATA_OUTPUT / 'analysis',
    )


if __name__ == '__main__':
    main()
