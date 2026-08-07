"""
Step 3a of the priority E/Z re-run: generate Stage 1 (_opt.gjf) + Stage 2
(_nbo.gjf) inputs for the *other* isomer of each of the 9 priority
molecules, without touching the already-run isomer's existing dft_opt/
directory.

beckmann_nbo.inputs.prepare_opt() writes every isomer present in the SDF it's
given for a molecule id, so each molecule's SDF is filtered down to a
single-row temp SDF (just the other isomer) before calling prepare_opt() --
otherwise it would try to rewrite the already-completed isomer's directory
too.

These 9 directories are written to their own local dir, data/output/dft_opt_ez_other/
(same pattern as the existing dft_opt_stepscan/ side-directory), NOT into
data/output/dft_opt/ alongside the already-run isomers. This matters for
hpc_sync.py: its --mol filter globs mol_{id}_*/ by NUMERIC ID ONLY, isomer-
agnostic (see beckmann_nbo.hpc.mol_dirs()) -- if the new mol_002_Z lived next
to the already-completed mol_002_E in the same --dir, "--mol 002 submit-opt"
would resubmit and overwrite mol_002_E's already-finished Stage 1 job too.
Keeping the 9 new isomers in a separate --dir sidesteps that entirely.
"""
import tempfile
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem

from beckmann_nbo.config import DATA_OUTPUT
from beckmann_nbo.inputs import prepare_opt

from ez_energy_check import energy_table

BEST_AIMNET_SDF = DATA_OUTPUT / "aimnet_optimized" / "best_aimnet_optimized.sdf"
DFT_OPT_DIR = DATA_OUTPUT / "dft_opt_ez_other"


def other_isomer_names(rows: list[dict]) -> list[str]:
    return [f"{r['mol_id']}_{r['other_isomer']}" for r in rows if r["other_isomer"]]


def main(mol_ids=None) -> None:
    from ez_energy_check import ALREADY_RUN
    rows = energy_table(mol_ids=mol_ids if mol_ids is not None else tuple(ALREADY_RUN))
    names = other_isomer_names(rows)
    print(f"Other isomers to prepare: {names}")

    suppl = Chem.SDMolSupplier(str(BEST_AIMNET_SDF), removeHs=False)
    mols_by_name = {m.GetProp("_Name"): m for m in suppl if m is not None}

    for name in names:
        mol_dir = DFT_OPT_DIR / name
        if mol_dir.exists():
            print(f"  SKIP {name}: {mol_dir} already exists -- not overwriting")
            continue
        mol = mols_by_name.get(name)
        if mol is None:
            print(f"  SKIP {name}: not found in {BEST_AIMNET_SDF}")
            continue

        mol_id = name.rsplit("_", 1)[0]
        test_id = mol_id.split("_")[1]

        with tempfile.NamedTemporaryFile(suffix=".sdf", mode="w", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        writer = Chem.SDWriter(str(tmp_path))
        writer.write(mol)
        writer.close()

        prepare_opt(tmp_path, DFT_OPT_DIR, test_ids={test_id})
        tmp_path.unlink()

    print("\nDone. Upload + submit-opt each molecule via hpc_sync.py, using --dir to keep")
    print("these isolated from the already-completed isomers in data/output/dft_opt/:")
    for name in names:
        mol_id = name.rsplit("_", 1)[0]
        test_id = mol_id.split("_")[1]
        print(f"  python packages/beckmann-nbo/scripts/hpc_sync.py --dir {DFT_OPT_DIR} --mol {test_id} upload")
        print(f"  python packages/beckmann-nbo/scripts/hpc_sync.py --dir {DFT_OPT_DIR} --mol {test_id} submit-opt")


if __name__ == "__main__":
    main()
