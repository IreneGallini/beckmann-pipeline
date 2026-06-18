"""
Step 3: Prepare Gaussian 16 input files for DFT/NBO7 analysis

Reads best_aimnet_optimized.sdf and writes, for each structure:
  - {mol_name}.gjf         Gaussian 16 input (single-point + NBO7)
  - {mol_name}_submit.sh   SLURM job submission template

Also writes data/output/dft_inputs/submit_all.sh to batch-submit from the cluster.

HPC hand-off (manual):
  scp -r data/output/dft_inputs/ user@cluster:~/beckmann/
  cd ~/beckmann/dft_inputs && bash submit_all.sh
  scp user@cluster:~/beckmann/dft_inputs/**/*.log data/output/dft_inputs/
"""

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem

# ── DFT settings ──────────────────────────────────────────────────────────────
FUNCTIONAL   = "B3LYP"
BASIS        = "6-311+G(d,p)"
JOB_TYPE     = "sp"       # "sp" = single-point on AIMNet2 geometry; "opt" = DFT re-opt
NPROC        = 8
MEM_GB       = 16
CHARGE       = 0
MULTIPLICITY = 1

# NBO7 keywords fed to the $NBO…$END section (pop=nboread activates this block).
# E2PERT: donor–acceptor perturbation table — key for σ*/CN-handoff analysis
# BNDIDX: Wiberg bond indices
# NBOSUM: NBO summary table
NBO_KEYWORDS = "E2PERT BNDIDX NBOSUM"
# ──────────────────────────────────────────────────────────────────────────────

OXIME_PAT = Chem.MolFromSmarts('[C:1]=[N:2]-[OH1:3]')


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


def _slurm(name: str) -> str:
    return (
        "#!/bin/bash\n"
        f"#SBATCH --job-name={name}\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task={NPROC}\n"
        f"#SBATCH --mem={MEM_GB}GB\n"
        "#SBATCH --time=24:00:00\n"
        f"#SBATCH --output={name}_%j.out\n"
        f"#SBATCH --error={name}_%j.err\n"
        "# Adjust the module name and add --partition / --account as needed:\n"
        "module load gaussian/16\n"
        "\n"
        f"g16 < {name}.gjf > {name}.log\n"
    )


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
        (mol_dir / f"{name}_submit.sh").write_text(_slurm(name))

        print(f"  {name:<24} {len(coords):>5}  {oxime_label:>18}  ✓")

    # Master submission script (run from dft_inputs/ on the cluster)
    submit_all = (
        "#!/bin/bash\n"
        "# Run from the dft_inputs/ directory on the HPC cluster:\n"
        "#   bash submit_all.sh\n"
        "for dir in */; do\n"
        '    name="${dir%/}"\n'
        '    cd "$dir" && sbatch "${name}_submit.sh" && cd ..\n'
        "done\n"
    )
    (outdir / "submit_all.sh").write_text(submit_all)

    print(f"\n{len(mols)} input sets → {outdir}")
    print(f"Master submission script → {outdir / 'submit_all.sh'}")
    print(
        "\nTo transfer to HPC:\n"
        f"  scp -r {outdir}/ user@cluster:~/beckmann/dft_inputs/\n"
        "  cd ~/beckmann/dft_inputs && bash submit_all.sh"
    )


if __name__ == "__main__":
    main()