"""`beckmann-pyscf optimize` -- conformers SDF -> AIMNet2-optimized geometry
only, standalone (accepts any conformers SDF, not just `conformers`'
output). Thin dispatch onto beckmann_pyscf.pipeline.optimize.run_optimize()
plus beckmann_core.classical.get_oxime_atoms() for the atom-map sanity
check printed below -- this is the "did anything get missed" checkpoint: a
bad/ambiguous atom map here would otherwise only surface as a crash in
`scan`.
"""
import sys
from pathlib import Path

from rdkit import Chem

from beckmann_core.classical import get_oxime_atoms
from beckmann_pyscf.pipeline.optimize import run_optimize


def cmd_optimize(args) -> None:
    conformers_sdf = Path(args.conformers_sdf)
    if not conformers_sdf.exists():
        print(f"ERROR: {conformers_sdf} not found", file=sys.stderr)
        sys.exit(1)

    name = args.name or conformers_sdf.stem
    out = Path(args.out) if args.out else Path("beckmann_pyscf_runs") / name / "optimized"
    out.mkdir(parents=True, exist_ok=True)

    print(f"AIMNet2 geometry optimization for {name!r}...")
    try:
        mol = run_optimize(conformers_sdf, out)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    best_path = out / "best.sdf"
    with Chem.SDWriter(str(best_path)) as writer:
        writer.write(mol)

    energy_ev = float(mol.GetProp("E_aimnet2_eV")) if mol.HasProp("E_aimnet2_eV") else None
    atom_ids = get_oxime_atoms(mol)

    print(f"\nAIMNet2 energy: {energy_ev:.6f} eV" if energy_ev is not None else "\nAIMNet2 energy: unavailable")
    if atom_ids is None:
        print("Oxime atom map: NOT FOUND -- get_oxime_atoms() could not identify the C=N-O / aryl / alkyl atoms.")
        print("  A downstream `scan` call will need --ci/--ni/--oi/--c-aryl/--c-alkyl passed explicitly.")
    else:
        cox, nox, oox, c_aryl, c_alkyl = atom_ids
        print(
            "Oxime atom map (0-based RDKit indices): "
            f"C={cox} N={nox} O={oox} c_aryl={c_aryl} c_alkyl={c_alkyl}"
        )
        print(
            "  (1-based, for `scan`'s --ci/--ni/--oi/--c-aryl/--c-alkyl overrides): "
            f"--ci {cox + 1} --ni {nox + 1} --oi {oox + 1} --c-aryl {c_aryl + 1} --c-alkyl {c_alkyl + 1}"
        )

    print(f"\nWrote: {best_path}")
