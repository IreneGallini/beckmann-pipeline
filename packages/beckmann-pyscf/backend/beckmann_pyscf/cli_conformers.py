"""`beckmann-pyscf conformers` -- SMILES -> Auto3D conformers only, the
first pipeline stage in isolation. Thin dispatch onto
beckmann_pyscf.pipeline.conformers's validate_smiles()/smiles_to_oxime_smi()/
run_conformers(), unchanged.
"""
import sys
from pathlib import Path

from beckmann_pyscf.pipeline import validate_smiles
from beckmann_pyscf.pipeline.conformers import run_conformers, smiles_to_oxime_smi


def cmd_conformers(args) -> None:
    error = validate_smiles(args.smiles)
    if error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

    out = Path(args.out) if args.out else Path("beckmann_pyscf_runs") / args.name / "conformers"
    out.parent.mkdir(parents=True, exist_ok=True)

    print(f"Generating conformers for {args.name!r}...")
    smi_path = smiles_to_oxime_smi(args.smiles, out.parent, mol_name=args.name)
    sdf_path = run_conformers(smi_path, out)

    print(f"\nWrote: {sdf_path}")
