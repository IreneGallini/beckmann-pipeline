"""`beckmann-pyscf predict` -- SMILES -> conformers -> AIMNet2 opt -> PySCF
wCNmax scan -> R/F prediction, the full pipeline in one call, with [1/4]..
[4/4] progress printed between stages (matching beckmann-nbo predict's UX).

Deliberately re-implements the same stage sequence
beckmann_pyscf.pipeline.predict() already runs, calling the exact same
underlying functions in the exact same order, rather than calling
pipeline.predict() as one opaque call -- this is what gives per-stage
progress lines for the conformers/optimize stages (the scan stage already
prints its own per-point progress internally, from wcnmax_pyscf.py, and
that's reused unchanged here). pipeline.predict() itself is untouched and
still covered by its own tests; this is a second, CLI-specific caller of
the same building blocks, not a replacement.
"""
import sys
from pathlib import Path

from rdkit import Chem

from beckmann_core.classical import get_oxime_atoms
from beckmann_pyscf.cli_common import write_series_csv, write_summary
from beckmann_pyscf.pipeline import validate_smiles
from beckmann_pyscf.pipeline.conformers import run_conformers, smiles_to_oxime_smi
from beckmann_pyscf.pipeline.optimize import run_optimize
from beckmann_pyscf.pipeline.plot import plot_wcnmax
from beckmann_pyscf.pipeline.predict import predict_outcome
from beckmann_pyscf.pipeline.wcnmax_pyscf import run_scan_series


def cmd_predict(args) -> None:
    error = validate_smiles(args.smiles)
    if error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

    name = args.name
    out = Path(args.out) if args.out else Path("beckmann_pyscf_runs") / name
    out.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] Generating conformers for {name!r}...")
    smi_path = smiles_to_oxime_smi(args.smiles, out, mol_name=name)
    conformers_sdf = run_conformers(smi_path, out / "conformers")

    print(f"[2/4] AIMNet2 geometry optimization for {name!r}...")
    mol = run_optimize(conformers_sdf, out / "optimized")
    energy_ev = float(mol.GetProp("E_aimnet2_eV")) if mol.HasProp("E_aimnet2_eV") else None

    atom_ids = get_oxime_atoms(mol)
    if atom_ids is None:
        print(
            f"ERROR: could not resolve the oxime C=N-O / aryl / alkyl atom map on the "
            f"optimized structure for {name!r}",
            file=sys.stderr,
        )
        sys.exit(1)
    cox, nox, oox, c_aryl, c_alkyl = atom_ids

    print(f"[3/4] Running PySCF wCNmax scan for {name!r}...")
    series = run_scan_series(
        mol, cox + 1, nox + 1, oox + 1, c_aryl + 1, c_alkyl + 1, name,
    )

    print(f"[4/4] Predicting outcome for {name!r}...")
    prediction, minimum = predict_outcome(name, series)

    optimized_path = out / "optimized.sdf"
    with Chem.SDWriter(str(optimized_path)) as writer:
        writer.write(mol)

    series_path = out / "wcnmax_series.csv"
    write_series_csv(series, series_path)
    summary_path = out / "summary.txt"
    summary_text = write_summary(
        summary_path, name=name, prediction=prediction, minimum=minimum, energy_ev=energy_ev,
    )

    print(f"\n{summary_text}")
    print("Wrote:")
    print(f"  {optimized_path}")
    print(f"  {series_path}")
    print(f"  {summary_path}")

    if args.plot:
        plot_path = out / "wcnmax_vs_rno.png"
        plot_wcnmax(series, minimum, plot_path, name=name)
        print(f"  {plot_path}")
