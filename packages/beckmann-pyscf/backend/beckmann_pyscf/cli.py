"""`beckmann-pyscf` -- terminal CLI for the open-source, HPC-free
AIMNet2+PySCF Beckmann rearrangement pipeline. The recommended way to run
this pipeline locally; the Flask web app in this same package (backend/
app.py) stays available separately as a hosted prototype for external
collaborators without a local Python environment -- see
packages/beckmann-pyscf/README.md.

Four subcommands: `predict` runs the full pipeline in one call; `conformers`/
`optimize`/`scan` run one stage at a time so a stage's output can be
inspected on its own before trusting a full predict() run. Every subcommand
is a thin dispatch onto existing beckmann_pyscf.pipeline functions -- see
each cli_*.py module's docstring for exactly which ones.

Mirrors the sibling Gaussian/NBO7 package's CLI argparse-with-subparsers
pattern. No global --dry-run/--mol/--dir here: those are that package's
HPC-job-targeting concepts, which this synchronous, no-HPC pipeline has no
equivalent of -- each subcommand owns its own flags instead.
"""
import argparse

from beckmann_pyscf.cli_conformers import cmd_conformers
from beckmann_pyscf.cli_optimize import cmd_optimize
from beckmann_pyscf.cli_predict import cmd_predict
from beckmann_pyscf.cli_scan import cmd_scan


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="beckmann-pyscf",
        description="Open-source, HPC-free Beckmann rearrangement pipeline -- SMILES to R/F prediction.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  beckmann-pyscf predict --smiles 'O=C1CCC2=C1C=CC=C2' --name test1 --plot\n"
            "  beckmann-pyscf predict --csv molecules.csv --plot\n"
            "  beckmann-pyscf conformers --smiles 'O=C1CCC2=C1C=CC=C2' --name test1\n"
            "  beckmann-pyscf optimize --conformers-sdf beckmann_pyscf_runs/test1/conformers/test1_out.sdf\n"
            "  beckmann-pyscf scan --sdf beckmann_pyscf_runs/test1_out/optimized/best.sdf --plot\n"
        ),
    )

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    p_predict = sub.add_parser("predict", help="SMILES (or --csv) -> conformers -> AIMNet2 opt -> PySCF wCNmax scan -> R/F prediction")
    p_predict.add_argument("--smiles", help="A single ketone SMILES string")
    p_predict.add_argument("--name", default="query", help="Short id for --smiles (default: 'query')")
    p_predict.add_argument("--csv", help="CSV with 'id'/'SMILES' columns (same shape as data/input/benchmark.csv) for a whole batch")
    p_predict.add_argument("--out", help="Output directory (default: ./beckmann_pyscf_runs/<name>/, or ./beckmann_pyscf_runs/<id>/ per --csv row)")
    p_predict.add_argument("--plot", action="store_true", help="Also write wcnmax_vs_rno.png")
    p_predict.set_defaults(func=cmd_predict)

    p_conformers = sub.add_parser("conformers", help="SMILES -> Auto3D conformers only")
    p_conformers.add_argument("--smiles", required=True, help="A single ketone SMILES string")
    p_conformers.add_argument("--name", default="query", help="Short id for this molecule (default: 'query')")
    p_conformers.add_argument("--out", help="Output directory (default: ./beckmann_pyscf_runs/<name>/conformers/)")
    p_conformers.set_defaults(func=cmd_conformers)

    p_optimize = sub.add_parser("optimize", help="conformers SDF -> AIMNet2-optimized geometry only")
    p_optimize.add_argument("--conformers-sdf", required=True, help="Path to a conformers SDF (any SDF, not necessarily from 'conformers')")
    p_optimize.add_argument("--name", help="Short id (default: derived from --conformers-sdf's filename)")
    p_optimize.add_argument("--out", help="Output directory (default: ./beckmann_pyscf_runs/<name>/optimized/)")
    p_optimize.set_defaults(func=cmd_optimize)

    p_scan = sub.add_parser("scan", help="optimized SDF -> PySCF wCNmax scan + R/F prediction only")
    p_scan.add_argument("--sdf", required=True, help="Path to an AIMNet2-optimized (or any 3D) single-molecule SDF")
    p_scan.add_argument("--name", help="Short id (default: derived from --sdf's filename)")
    p_scan.add_argument("--out", help="Output directory (default: ./beckmann_pyscf_runs/<name>/scan/)")
    p_scan.add_argument("--plot", action="store_true", help="Also write wcnmax_vs_rno.png")
    p_scan.add_argument("--ci", type=int, help="1-based oxime C atom index (override auto-detection)")
    p_scan.add_argument("--ni", type=int, help="1-based oxime N atom index (override auto-detection)")
    p_scan.add_argument("--oi", type=int, help="1-based oxime O atom index (override auto-detection)")
    p_scan.add_argument("--c-aryl", type=int, dest="c_aryl", help="1-based aryl migrating-group atom index (override)")
    p_scan.add_argument("--c-alkyl", type=int, dest="c_alkyl", help="1-based alkyl migrating-group atom index (override)")
    p_scan.add_argument("--r-min", type=float, dest="r_min", help="Scan window start, Angstrom (default: 1.50)")
    p_scan.add_argument("--r-max", type=float, dest="r_max", help="Scan window end, Angstrom (default: 1.80)")
    p_scan.add_argument("--r-step", type=float, dest="r_step", help="Scan step size, Angstrom (default: 0.05)")
    p_scan.set_defaults(func=cmd_scan)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
