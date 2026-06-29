"""
Step 3: Prepare Gaussian 16 input files for DFT/NBO7 analysis

Reads best_aimnet_optimized.sdf and writes, for each structure:
  - {mol_name}.gjf    Gaussian 16 input (single-point + NBO7 on AIMNet2 geometry)

Output: data/output/dft_inputs/

Submission on Citadel (shared server, no SLURM):
  python scripts/hpc_sync.py --dir data/output/dft_inputs upload
  python scripts/hpc_sync.py --dir data/output/dft_inputs submit-sp
  python scripts/hpc_sync.py --dir data/output/dft_inputs download

Note: for DFT geometry optimisation before NBO, use scripts/05_prepare_test_opt.py
instead. That generates a two-stage opt→NBO workflow in data/output/dft_opt_test/.
"""

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem

# ── DFT settings ──────────────────────────────────────────────────────────────
FUNCTIONAL   = "wB97XD"
BASIS        = "6-311+G(d,p)"
JOB_TYPE     = "sp"       # "sp" = single-point on AIMNet2 geometry; "opt" = DFT re-opt
NPROC        = 8
MEM_GB       = 16
CHARGE       = 1   # protonated activated oxime (C=N-[OH2+])
MULTIPLICITY = 1

# NBO7 keywords fed to the $NBO…$END section (pop=nboread activates this block).
# E2PERT: donor–acceptor perturbation table — key for σ*/CN-handoff analysis
# BNDIDX: Wiberg bond indices
# NBOSUM: NBO summary table
NBO_KEYWORDS = "E2PERT BNDIDX NBOSUM"
# ──────────────────────────────────────────────────────────────────────────────

OXIME_PAT = Chem.MolFromSmarts('[C:1]=[N:2]-[O+:3]')  # matches [OH2+] in protonated activated oxime


def _gjf(name: str, coords: list[tuple], oxime_label: str) -> str:
    route = f"#p {FUNCTIONAL}/{BASIS} {JOB_TYPE} pop=nboread"
    header = (
        f"%chk={name}.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"{route}\n"
        f"\n"
        f"{name}  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
    )
    coord_block = "\n".join(
        f"{sym:<3}  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
        for sym, x, y, z in coords
    )
    nbo_block = f"\n$NBO {NBO_KEYWORDS} $END\n\n\n"
    return header + coord_block + nbo_block



def main() -> None:
    root   = Path(__file__).parent.parent
    sdf    = root / "data" / "output" / "aimnet_optimized" / "best_aimnet_optimized.sdf"
    outdir = root / "data" / "output" / "dft_inputs"
    outdir.mkdir(parents=True, exist_ok=True)

    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    mols  = [m for m in suppl if m is not None]

    print(f"\n{'Name':<24} {'Atoms':>5}  {'Oxime (1-based)':>18}  {'Files'}")
    print("-" * 62)

    for mol in mols:
        name = mol.GetProp("_Name")
        conf = mol.GetConformer()

        coords = [
            (atom.GetSymbol(),
             *conf.GetAtomPosition(i))
            for i, atom in enumerate(mol.GetAtoms())
        ]

        # Annotate oxime atoms (1-based for Gaussian / NBO output reference)
        match = mol.GetSubstructMatch(OXIME_PAT)
        if match:
            ci, ni, oi = (idx + 1 for idx in match)
            oxime_label = f"[oxime: C{ci}=N{ni}-O{oi}]"
        else:
            oxime_label = "[oxime: not found]"

        mol_dir = outdir / name
        mol_dir.mkdir(exist_ok=True)

        (mol_dir / f"{name}.gjf").write_text(_gjf(name, coords, oxime_label))
        print(f"  {name:<24} {len(coords):>5}  {oxime_label:>18}  ✓")

    print(f"\n{len(mols)} .gjf files → {outdir}")
    print(
        "\nTo submit on Citadel:\n"
        "  python scripts/hpc_sync.py --dir data/output/dft_inputs upload\n"
        "  python scripts/hpc_sync.py --dir data/output/dft_inputs submit-sp"
    )


if __name__ == "__main__":
    main()