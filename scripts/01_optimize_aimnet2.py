"""Stage 0b: pick the lowest-energy conformer per isomer, AIMNet2/ASE geometry
optimization. Needs 00_smiles_to_conformers.py to have already run for this
MOL_NAME (reads its "*_out.sdf" from data/output/query_predictions/<id>/conformers/).

Edit MOL_NAME below to match what you used in 00, then:
    python 01_optimize_aimnet2.py
"""
import sys

from _common import sanitize_id, workdir_for

from beckmann_core.optimize import select_and_optimize

MOL_NAME = "test1"


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    workdir = workdir_for(mol_id)
    conformers_dir = workdir / "conformers"

    # Auto3D writes into a timestamped subdirectory it creates itself
    # (conformers_dir/mol_<id>_<timestamp>/mol_<id>_out.sdf), not directly
    # into conformers_dir.
    sdf_files = sorted(conformers_dir.glob("*/*_out.sdf"))
    if not sdf_files:
        print(
            f"ERROR: no '*/*_out.sdf' found in {conformers_dir} -- run "
            f"00_smiles_to_conformers.py for '{MOL_NAME}' first.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"AIMNet2 geometry optimization for {mol_id}...")
    best_sdf, best_per_substrate_sdf = select_and_optimize(sdf_files[-1], workdir / "aimnet_optimized")
    print(f"\nWrote {best_sdf}")
    print(f"Wrote {best_per_substrate_sdf}")
    print(f"\nNext: edit MOL_NAME in 02_prepare_stage12_inputs.py to '{MOL_NAME}' and run it.")


if __name__ == "__main__":
    main()
