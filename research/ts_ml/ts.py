"""
Locate transition states via Gaussian QST2/QST3, verify with frequency analysis,
and confirm connectivity with IRC.

Parallel to inputs.py/scan.py in structure. Endpoints (reactant + product/
intermediate) come from ts_ml.ts_products, which guarantees identical
atom ordering between them -- required for QST2/QST3, which interpolates
atom-by-atom between the supplied Cartesian structures.

Route lines follow the same FUNCTIONAL/BASIS/SOLVENT/CHARGE/MULTIPLICITY
constants as the rest of the DFT pipeline (beckmann/config.py) -- nothing here
duplicates or overrides them.

Output: data/output/dft_opt/{mol}/{mol}_{label}.gjf
"""
import re
from pathlib import Path

from rdkit import Chem

from beckmann_core.constants import CHARGE, MULTIPLICITY
from beckmann_nbo.config import (
    DATA_OUTPUT,
    FUNCTIONAL, BASIS, NPROC, MEM_GB,
    NBO_KEYWORDS, SOLVENT,
)

TS_KEYWORDS = "calcfc"  # appended inside opt=(qstN,...); noeigentest not needed for
                          # QST2/QST3 (unlike a bare opt=ts guess) since QSTn's own
                          # linear-synchronous-transit step supplies a reasonable
                          # starting Hessian direction, not just a raw guess geometry.

IMAG_FREQ_RE = re.compile(r"(\d+)\s+imaginary frequenc")
FREQ_BLOCK_RE = re.compile(r"^\s*Frequencies\s+--\s+(-?\d+\.\d+)")


def mol_to_atom_tuples(mol: Chem.Mol) -> list[tuple[str, float, float, float]]:
    """RDKit Mol (with a conformer) -> [(symbol, x, y, z), ...], 1 tuple per atom
    in RDKit atom-index order (0-based index i = atom i+1 in Gaussian's numbering).
    """
    conf = mol.GetConformer()
    return [
        (atom.GetSymbol(), *conf.GetAtomPosition(atom.GetIdx()))
        for atom in mol.GetAtoms()
    ]


def _coord_block(atoms: list[tuple[str, float, float, float]]) -> str:
    return "\n".join(
        f"{sym:<3}  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
        for sym, x, y, z in atoms
    )


def build_qst_gjf(
    mol: str,
    label: str,
    endpoints: list[tuple[str, list[tuple[str, float, float, float]]]],
) -> str:
    """QST2 (2 endpoints) or QST3 (3, third is a TS guess) input.

    `endpoints` = [(title, atom_tuples), ...], in order: reactant, product
    (, ts_guess). Each must have identical atom count/order (guaranteed by
    ts_ml.ts_products, which builds every endpoint from a copy of the
    same reactant Mol). Route includes freq (combined opt+freq job, standard
    practice) and pop=nbo7read so the converged TS geometry gets NBO/CMO
    descriptors directly comparable to the rest of the pipeline's data.
    """
    n = len(endpoints)
    if n not in (2, 3):
        raise ValueError(f"QST needs 2 or 3 endpoints, got {n}")
    qst_kw = f"qst{n}"

    header = (
        f"%chk={mol}_{label}.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} opt=({qst_kw},{TS_KEYWORDS}) freq "
        f"pop=nbo7read {SOLVENT} scf=(tight,xqc)\n"
    )

    body = ""
    for title, atoms in endpoints:
        body += (
            f"\n{mol} {label} -- {title}\n"
            f"\n"
            f"{CHARGE} {MULTIPLICITY}\n"
            f"{_coord_block(atoms)}\n"
        )

    return header + body + f"\n$NBO {NBO_KEYWORDS} $END\n\n\n"


def build_irc_gjf(mol: str, label: str, ts_chk_stem: str) -> str:
    """IRC job chained off a verified TS's .chk, both directions combined."""
    return (
        f"%chk={mol}_{label}.chk\n"
        f"%oldchk={ts_chk_stem}.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} irc=(calcfc,forward,reverse) "
        f"geom=checkpoint guess=read {SOLVENT}\n"
        f"\n"
        f"{mol} {label} IRC\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        f"\n\n"
    )


def verify_ts(log_path: Path, atoms_of_interest: dict[str, int]) -> dict:
    """Human-review report for a completed TS .gjf log: termination status,
    imaginary-frequency count, and the imaginary mode's per-atom displacement
    ranking. Does NOT auto-decide pass/fail -- prints for a human to judge
    against the expected reaction coordinate, same precedent as JOB_ISSUES.md.

    `atoms_of_interest` = {label: 1-based_atom_index}, e.g.
    {"ni": 12, "oi": 13, "ci": 11, "c_aryl": 6} for the rearrangement channel.
    """
    text = log_path.read_text()
    lines = text.splitlines()

    normal_termination = "Normal termination of Gaussian 16" in lines[-1] if lines else False

    imag_match = IMAG_FREQ_RE.search(text)
    n_imaginary = int(imag_match.group(1)) if imag_match else None

    # Locate the first "Frequencies --" line with a negative value (the imaginary
    # mode). Gaussian sorts modes ascending and prints 3 side-by-side columns per
    # block, so a real (single) imaginary mode is always the leftmost column of
    # the first block -- each displacement row is "Atom AN X Y Z X Y Z X Y Z"
    # (up to 11 fields), and we only want the first mode's X/Y/Z (columns 3-5).
    displacement_by_atom: dict[int, float] = {}
    for i, line in enumerate(lines):
        m = FREQ_BLOCK_RE.match(line)
        if m and float(m.group(1)) < 0:
            # displacement rows start after "Atom  AN  X  Y  Z ..." header, a few lines down
            for j in range(i, min(i + 200, len(lines))):
                row = lines[j].split()
                if len(row) >= 5 and row[0].isdigit() and row[1].isdigit():
                    atom_idx = int(row[0])
                    dx, dy, dz = (float(v) for v in row[2:5])
                    displacement_by_atom[atom_idx] = (dx**2 + dy**2 + dz**2) ** 0.5
                elif displacement_by_atom:
                    break  # ran past the displacement block
            break

    ranked = sorted(displacement_by_atom.items(), key=lambda t: -t[1])
    top_atoms = {idx for idx, _ in ranked[:6]}
    of_interest_ranks = {
        label: (idx, displacement_by_atom.get(idx))
        for label, idx in atoms_of_interest.items()
    }

    report = {
        "normal_termination": normal_termination,
        "n_imaginary": n_imaginary,
        "top_displaced_atoms": ranked[:6],
        "atoms_of_interest_displacement": of_interest_ranks,
        "atoms_of_interest_in_top6": {
            label: idx in top_atoms for label, (idx, _) in of_interest_ranks.items()
        },
    }
    return report


def main() -> None:
    """Pilot scope: mol_002_E rearrangement TS (TS1_A1) QST2 input only."""
    from beckmann_nbo.descriptors import _load_mols

    mol = "mol_002_E"
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    aimnet_dir = DATA_OUTPUT / "aimnet_optimized"

    reactant = _load_mols()[mol]
    product_sdf = aimnet_dir / f"{mol}_product_rearr.sdf"
    product = next(Chem.SDMolSupplier(str(product_sdf), removeHs=False))

    endpoints = [
        ("reactant (R_A1)", mol_to_atom_tuples(reactant)),
        ("rearrangement product (P_A1, AIMNet2 guess)", mol_to_atom_tuples(product)),
    ]
    gjf_text = build_qst_gjf(mol, "ts1_a1", endpoints)

    out_path = dft_opt_dir / mol / f"{mol}_ts1_a1.gjf"
    out_path.write_text(gjf_text)
    print(f"-- {mol}: TS1_A1 QST2 input -> {out_path}")


def print_verification_report(mol: str, label: str, report: dict) -> None:
    print(f"-- {mol} {label} TS verification (human review required, not auto-decided) --")
    print(f"   Normal termination: {report['normal_termination']}")
    print(f"   Imaginary frequencies: {report['n_imaginary']} (expect exactly 1)")
    print(f"   Top displaced atoms (1-based, |displacement|): {report['top_displaced_atoms']}")
    print(f"   Atoms of interest in top-6 displaced: {report['atoms_of_interest_in_top6']}")
    print("   -> Confirm the imaginary mode's dominant atoms match the expected "
          "reaction coordinate before treating this as a verified TS.")


if __name__ == "__main__":
    main()
