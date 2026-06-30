"""
Prepare two-stage Gaussian input files for DFT geometry optimisation + NBO analysis.

For the test set (molecules 002, 006, 020, 021, E and Z isomers), generates
two .gjf files per structure in data/output/dft_opt/{name}/:

  Stage 1 — {name}_opt.gjf
      wB97XD/6-311+G(d,p) opt
      Starting geometry: AIMNet2-optimised coordinates from best_aimnet_optimized.sdf
      Output: {name}_opt.chk  (contains DFT-optimised geometry)

  Stage 2 — {name}_nbo.gjf
      wB97XD/6-311+G(d,p) sp pop=nboread geom=checkpoint guess=read
      Reads geometry from {name}_opt.chk — run AFTER Stage 1 completes

These are separate jobs so a NBO7 failure does not destroy the optimised
geometry already stored in the checkpoint.

Submission on Citadel:
  python scripts/dft/hpc_sync.py --mol 002 upload
  python scripts/dft/hpc_sync.py --mol 002 submit-opt
  python scripts/dft/hpc_sync.py status
  python scripts/dft/hpc_sync.py --mol 002 submit-nbo   # after Stage 1 finishes
  python scripts/dft/hpc_sync.py --mol 002 download
"""

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem

# ── Settings ──────────────────────────────────────────────────────────────────
TEST_IDS     = {"002", "006", "020", "021"}
FUNCTIONAL   = "wB97XD"
BASIS        = "6-311+G(d,p)"
NPROC        = 8
MEM_GB       = 16
CHARGE       = 1
MULTIPLICITY = 1
# CMO requires NBO7 linked into the Gaussian build (NBOEXE mechanism).
# Citadel's g16 does not support NBOEXE — l607.exe (NBO 3.1) runs regardless.
# Remove CMO until the admin links NBO7 into g16 (or provides a patched l607.exe).
# With the current setup, E2PERT+BNDIDX+NBOSUM enable Psi and d/dR descriptors.
# Lambda and wCNmax (CMO-based) are deferred.
NBO_KEYWORDS = "E2PERT BNDIDX NBOSUM"
# ──────────────────────────────────────────────────────────────────────────────

OXIME_PAT = Chem.MolFromSmarts('[C:1]=[N:2]-[O+:3]')


def _opt_gjf(name: str, coords: list[tuple], oxime_label: str) -> str:
    """Stage 1: geometry optimisation — no NBO block."""
    return (
        f"%chk={name}_opt.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} opt\n"
        f"\n"
        f"{name} opt  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        + "\n".join(
            f"{sym:<3}  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
            for sym, x, y, z in coords
        )
        + "\n\n\n"
    )


def _nbo_gjf(name: str, oxime_label: str) -> str:
    """Stage 2: NBO single-point — reads DFT geometry from opt checkpoint."""
    return (
        f"%chk={name}_nbo.chk\n"
        f"%oldchk={name}_opt.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} sp pop=nboread geom=checkpoint guess=read\n"
        f"\n"
        f"{name} NBO  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        f"\n"
        f"$NBO {NBO_KEYWORDS} $END\n"
        f"\n\n"
    )


def _scan_gjf(name: str, ni: int, oi: int, oxime_label: str) -> str:
    """Stage 3: relaxed N-O bond scan — 5 points (R to R+0.4 A), NBO at each.

    Section order matters: ModRedundant specs come immediately after charge/mult,
    then the NBO block as a second additional input section.
    """
    return (
        f"%chk={name}_scan.chk\n"
        f"%oldchk={name}_opt.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} opt=(ModRedundant,MaxCycles=200) pop=nboread geom=checkpoint guess=read\n"
        f"\n"
        f"{name} scan  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        f"\n"
        f"B {ni} {oi} S 4 0.1\n"
        f"\n"
        f"$NBO {NBO_KEYWORDS} $END\n"
        f"\n\n"
    )


def main() -> None:
    root   = Path(__file__).parent.parent.parent
    sdf    = root / "data" / "output" / "aimnet_optimized" / "best_per_substrate.sdf"
    outdir = root / "data" / "output" / "dft_opt"
    outdir.mkdir(parents=True, exist_ok=True)

    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    mols  = [m for m in suppl if m is not None]
    test_mols = [m for m in mols if m.GetProp("_Name").split("_")[1] in TEST_IDS]

    print(f"\n{'Name':<24} {'Atoms':>5}  {'Oxime':>20}  Stage1  Stage2  Stage3")
    print("-" * 80)

    for mol in test_mols:
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
        (mol_dir / f"{name}_opt.gjf").write_text(_opt_gjf(name, coords, oxime_label))
        (mol_dir / f"{name}_nbo.gjf").write_text(_nbo_gjf(name, oxime_label))
        if match:
            (mol_dir / f"{name}_scan.gjf").write_text(_scan_gjf(name, ni, oi, oxime_label))
            scan_mark = "✓"
        else:
            print(f"  WARNING: {name} — oxime pattern not found, _scan.gjf not written")
            scan_mark = "✗"
        print(f"  {name:<24} {len(coords):>5}  {oxime_label:>20}   ✓      ✓      {scan_mark}")

    print(f"\n{len(test_mols)} structures written to {outdir}")
    print("\nSubmit on Citadel via hpc_sync.py:")
    print("  python scripts/dft/hpc_sync.py --mol 002 upload")
    print("  python scripts/dft/hpc_sync.py --mol 002 submit-opt")
    print("  python scripts/dft/hpc_sync.py status")
    print("  python scripts/dft/hpc_sync.py --mol 002 submit-nbo")
    print("  python scripts/dft/hpc_sync.py --mol 002 download")


if __name__ == "__main__":
    main()
