"""
Step 5 (test set): Prepare two-stage Gaussian input files for DFT geometry
optimisation + NBO analysis.

For molecules 002, 007, 020, 021 (E and Z isomers), generates two jobs per
structure:

  Stage 1 — {name}_opt.gjf / {name}_opt_submit.sh
      wB97XD/6-311+G(d,p) opt
      Starting geometry: AIMNet2-optimised coordinates from best_aimnet_optimized.sdf
      Output: {name}_opt.chk  (contains DFT-optimised geometry)

  Stage 2 — {name}_nbo.gjf / {name}_nbo_submit.sh
      wB97XD/6-311+G(d,p) sp pop=nboread geom=checkpoint guess=read
      Reads geometry from {name}_opt.chk — run AFTER Stage 1 completes
      Output: NBO7 donor/acceptor table, Wiberg bond indices, NBO summary

These are separate jobs so that a NBO7 failure (which is not uncommon) does
not destroy the optimised geometry already stored in the checkpoint.

Output: data/output/dft_opt_test/{name}/
"""

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem

# ── Settings ──────────────────────────────────────────────────────────────────
TEST_IDS     = {"002", "006", "020", "021"}   # mol IDs to process (zero-padded)
FUNCTIONAL   = "wB97XD"
BASIS        = "6-311+G(d,p)"
NPROC        = 8
MEM_GB       = 16
CHARGE       = 1    # protonated activated oxime (C=N-[OH2+])
MULTIPLICITY = 1
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


def _slurm(name: str, gjf: str, job_label: str) -> str:
    job_name = f"{name}_{job_label}"
    return (
        "#!/bin/bash\n"
        f"#SBATCH --job-name={job_name}\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task={NPROC}\n"
        f"#SBATCH --mem={MEM_GB}GB\n"
        "#SBATCH --time=24:00:00\n"
        f"#SBATCH --output={job_name}_%j.out\n"
        f"#SBATCH --error={job_name}_%j.err\n"
        "# Adjust module name and --partition / --account for your cluster:\n"
        "module load gaussian/16\n"
        "\n"
        f"g16 < {gjf} > {name}_{job_label}.log\n"
    )


def main() -> None:
    root   = Path(__file__).parent.parent
    sdf    = root / "data" / "output" / "aimnet_optimized" / "best_aimnet_optimized.sdf"
    outdir = root / "data" / "output" / "dft_opt_test"
    outdir.mkdir(parents=True, exist_ok=True)

    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    mols  = [m for m in suppl if m is not None]

    # Filter to test set
    test_mols = [m for m in mols if m.GetProp("_Name").split("_")[1] in TEST_IDS]

    print(f"\n{'Name':<24} {'Atoms':>5}  {'Oxime':>20}  Stage1  Stage2")
    print("-" * 72)

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

        opt_gjf  = f"{name}_opt.gjf"
        nbo_gjf  = f"{name}_nbo.gjf"

        (mol_dir / opt_gjf).write_text(_opt_gjf(name, coords, oxime_label))
        (mol_dir / nbo_gjf).write_text(_nbo_gjf(name, oxime_label))
        (mol_dir / f"{name}_opt_submit.sh").write_text(_slurm(name, opt_gjf, "opt"))
        (mol_dir / f"{name}_nbo_submit.sh").write_text(_slurm(name, nbo_gjf, "nbo"))

        print(f"  {name:<24} {len(coords):>5}  {oxime_label:>20}   ✓      ✓")

    # Two-stage cluster submission guide
    guide = (
        "#!/bin/bash\n"
        "# Two-stage submission for DFT test set.\n"
        "# Run from the dft_opt_test/ directory on the HPC cluster.\n"
        "#\n"
        "# Stage 1: submit all optimisations\n"
        "for dir in */; do\n"
        '    name="${dir%/}"\n'
        '    cd "$dir" && sbatch "${name}_opt_submit.sh" && cd ..\n'
        "done\n"
        "#\n"
        "# Stage 2: after ALL Stage 1 jobs finish, submit NBO single-points\n"
        "# (the _nbo.gjf reads the .chk from Stage 1, so Stage 1 must complete first)\n"
        "for dir in */; do\n"
        '    name="${dir%/}"\n'
        '    cd "$dir" && sbatch "${name}_nbo_submit.sh" && cd ..\n'
        "done\n"
    )
    (outdir / "submit_stages.sh").write_text(guide)

    print(f"\n{len(test_mols)} structures written to {outdir}")
    print("\nFiles per molecule:")
    print("  {name}_opt.gjf         Stage 1 — DFT geometry optimisation")
    print("  {name}_nbo.gjf         Stage 2 — NBO single-point (reads opt .chk)")
    print("  {name}_opt_submit.sh   SLURM script for Stage 1")
    print("  {name}_nbo_submit.sh   SLURM script for Stage 2")
    print("\nSubmission workflow (on cluster):")
    print("  scp -r data/output/dft_opt_test/ user@cluster:~/beckmann/")
    print("  cd ~/beckmann/dft_opt_test")
    print("  # Submit Stage 1 (opt) for all 8 structures:")
    print("  for dir in */; do name=${dir%/}; cd $dir && sbatch ${name}_opt_submit.sh && cd ..; done")
    print("  # After Stage 1 finishes, submit Stage 2 (NBO):")
    print("  for dir in */; do name=${dir%/}; cd $dir && sbatch ${name}_nbo_submit.sh && cd ..; done")
    print("\nFrom the NBO .log files you will want to parse:")
    print("  - 'NATURAL BOND ORBITALS' section (orbital character, e.g. CN-handoff)")
    print("  - 'E2PERT' table (donor→acceptor perturbation energies, key for σ* analysis)")
    print("  - 'BNDIDX' table (Wiberg bond indices for N-O, C=N, migrating C-C)")


if __name__ == "__main__":
    main()