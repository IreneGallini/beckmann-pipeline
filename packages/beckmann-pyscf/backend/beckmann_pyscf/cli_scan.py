"""`beckmann-pyscf scan` -- optimized SDF -> PySCF wCNmax scan + R/F
prediction only, standalone. Thin dispatch onto
beckmann_pyscf.pipeline.wcnmax_pyscf.run_scan_series() and
beckmann_pyscf.pipeline.predict.predict_outcome(), plus
beckmann_core.classical.get_oxime_atoms() for auto-detecting the atom map
when --ci/--ni/--oi/--c-aryl/--c-alkyl aren't all given explicitly.
"""
import sys
from pathlib import Path

from rdkit import Chem

from beckmann_core.classical import get_oxime_atoms
from beckmann_pyscf.cli_common import write_series_csv, write_summary
from beckmann_pyscf.pipeline.plot import plot_wcnmax
from beckmann_pyscf.pipeline.predict import predict_outcome
from beckmann_pyscf.pipeline.wcnmax_pyscf import run_scan_series

OVERRIDE_FIELDS = ["ci", "ni", "oi", "c_aryl", "c_alkyl"]


def _resolve_atom_indices(args, mol) -> tuple[int, int, int, int, int]:
    overrides = {f: getattr(args, f) for f in OVERRIDE_FIELDS}
    given = {f: v for f, v in overrides.items() if v is not None}

    if len(given) == len(OVERRIDE_FIELDS):
        return tuple(overrides[f] for f in OVERRIDE_FIELDS)
    if given:
        missing = [f for f in OVERRIDE_FIELDS if f not in given]
        print(
            f"ERROR: pass all five of --ci/--ni/--oi/--c-aryl/--c-alkyl together, "
            f"or none to auto-detect (missing: {', '.join('--' + f.replace('_', '-') for f in missing)})",
            file=sys.stderr,
        )
        sys.exit(1)

    atom_ids = get_oxime_atoms(mol)
    if atom_ids is None:
        print(
            "ERROR: get_oxime_atoms() could not identify the C=N-O / aryl / alkyl atoms on this "
            "structure -- pass --ci/--ni/--oi/--c-aryl/--c-alkyl explicitly (see `optimize`'s printed atom map).",
            file=sys.stderr,
        )
        sys.exit(1)
    cox, nox, oox, c_aryl, c_alkyl = atom_ids
    return cox + 1, nox + 1, oox + 1, c_aryl + 1, c_alkyl + 1


def cmd_scan(args) -> None:
    sdf_path = Path(args.sdf)
    if not sdf_path.exists():
        print(f"ERROR: {sdf_path} not found", file=sys.stderr)
        sys.exit(1)

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        print(f"ERROR: no readable structure in {sdf_path}", file=sys.stderr)
        sys.exit(1)
    mol = mols[0]

    name = args.name or sdf_path.stem
    out = Path(args.out) if args.out else Path("beckmann_pyscf_runs") / name / "scan"
    out.mkdir(parents=True, exist_ok=True)

    ci, ni, oi, c_aryl, c_alkyl = _resolve_atom_indices(args, mol)
    print(
        f"Atom map (1-based) for {name!r}: ci={ci} ni={ni} oi={oi} "
        f"c_aryl={c_aryl} c_alkyl={c_alkyl}"
    )

    scan_kwargs = {}
    if args.r_min is not None:
        scan_kwargs["r_min"] = args.r_min
    if args.r_max is not None:
        scan_kwargs["r_max"] = args.r_max
    if args.r_step is not None:
        scan_kwargs["r_step"] = args.r_step

    print(f"Running PySCF wCNmax scan for {name!r}...")
    series = run_scan_series(mol, ci, ni, oi, c_aryl, c_alkyl, name, **scan_kwargs)
    prediction, minimum = predict_outcome(name, series)

    series_path = out / "wcnmax_series.csv"
    write_series_csv(series, series_path)
    summary_path = out / "summary.txt"
    summary_text = write_summary(summary_path, name=name, prediction=prediction, minimum=minimum)

    print(f"\n{summary_text}")
    print(f"Wrote: {series_path}")
    print(f"Wrote: {summary_path}")

    if args.plot:
        plot_path = out / "wcnmax_vs_rno.png"
        plot_wcnmax(series, minimum, plot_path, name=name)
        print(f"Wrote: {plot_path}")
