"""
Prepare Gaussian 16 single-point + NBO7 input files for all benchmark molecules.

Reads best_aimnet_optimized.sdf and writes {mol_name}.gjf for each structure.
Output: data/output/dft_sp/

Use this for a fast NBO analysis directly on the AIMNet2-optimised geometry,
without a DFT geometry re-optimisation. For the full opt→NBO two-stage
workflow use scripts/dft/prepare_opt.py instead.

Submission on Citadel:
  python scripts/dft/hpc_sync.py --dir data/output/dft_sp upload
  python scripts/dft/hpc_sync.py --dir data/output/dft_sp submit-sp
  python scripts/dft/hpc_sync.py --dir data/output/dft_sp download
"""

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem

# ── DFT settings ──────────────────────────────────────────────────────────────
FUNCTIONAL   = "wB97XD"
BASIS        = "6-311+G(d,p)"
NPROC        = 8
MEM_GB       = 16
CHARGE       = 1   # protonated activated oxime (C=N-[OH2+])
MULTIPLICITY = 1
# CMO keyword requires NBO7 (separately licensed). Citadel runs NBO 3.1 (bundled
# with Gaussian 16), which does not support CMO. Add CMO here only if NBO7 is
# installed and linked on the cluster.
NBO_KEYWORDS = "E2PERT BNDIDX NBOSUM"
# ──────────────────────────────────────────────────────────────────────────────

OXIME_PAT = Chem.MolFromSmarts('[C:1]=[N:2]-[O+:3]')


def _gjf(name: str, coords: list[tuple], oxime_label: str) -> str:
    header = (
        f"%chk={name}.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} sp pop=nboread\n"
        f"\n"
        f"{name}  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
    )
    coord_block = "\n".join(
        f"{sym:<3}  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
        for sym, x, y, z in coords
    )
    return header + coord_block + f"\n$NBO {NBO_KEYWORDS} $END\n\n\n"


def main() -> None:
    root   = Path(__file__).parent.parent.parent
    sdf    = root / "data" / "output" / "aimnet_optimized" / "best_per_substrate.sdf"
    outdir = root / "data" / "output" / "dft_sp"
    outdir.mkdir(parents=True, exist_ok=True)

    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    mols  = [m for m in suppl if m is not None]

    print(f"\n{'Name':<24} {'Atoms':>5}  {'Oxime (1-based)':>18}")
    print("-" * 52)

    for mol in mols:
        name = mol.GetProp("_Name")
        conf = mol.GetConformer()
        coords = [
            (atom.GetSymbol(), *conf.GetAtomPosition(i))
            for i, atom in enumerate(mol.GetAtoms())
        ]
        match = mol.GetSubstructMatch(OXIME_PAT)
        if match:
            ci, ni, oi = (idx + 1 for idx in match)
            oxime_label = f"[oxime: C{ci}=N{ni}-O{oi}]"
        else:
            oxime_label = "[oxime: not found]"

        mol_dir = outdir / name
        mol_dir.mkdir(exist_ok=True)
        (mol_dir / f"{name}.gjf").write_text(_gjf(name, coords, oxime_label))
        print(f"  {name:<24} {len(coords):>5}  {oxime_label:>18}")

    print(f"\n{len(mols)} .gjf files → {outdir}")
    print(
        "\nTo submit on Citadel:\n"
        "  python scripts/dft/hpc_sync.py --dir data/output/dft_sp upload\n"
        "  python scripts/dft/hpc_sync.py --dir data/output/dft_sp submit-sp"
    )


if __name__ == "__main__":
    main()
