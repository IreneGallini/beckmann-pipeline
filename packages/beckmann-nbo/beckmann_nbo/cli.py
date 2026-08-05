"""`beckmann-nbo` -- CLI wrapper around this package's existing Gaussian/
NBO7/HPC modules, so a new user with their own SSH-reachable HPC cluster
(Gaussian16 + NBO7) can go SMILES -> prediction/plots without reading the
full monorepo CLAUDE.md. Every subcommand below is a thin dispatch onto
existing beckmann_core/beckmann_nbo functions -- see each cli_*.py module's
docstring for exactly which ones.

Mirrors beckmann_nbo/hpc.py::main()'s argparse-with-subparsers pattern.
"""
import argparse

from beckmann_nbo.cli_init import cmd_init
from beckmann_nbo.cli_predict import cmd_predict
from beckmann_nbo.cli_recover import cmd_recover
from beckmann_nbo.cli_report import cmd_report
from beckmann_nbo.cli_status import cmd_status_command
from beckmann_nbo.cli_verify import cmd_verify


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="beckmann-nbo",
        description="Gaussian/NBO7 Beckmann rearrangement pipeline -- SMILES to R/F prediction.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  beckmann-nbo init\n"
            "  beckmann-nbo verify\n"
            "  beckmann-nbo predict --smiles 'O=C1CCC2=C1C=CC=C2' --name test1\n"
            "  beckmann-nbo predict --continue qtest1 --dir data/output/query_predictions/qtest1/dft_opt\n"
            "  beckmann-nbo status --mol qtest1 --dir data/output/query_predictions/qtest1/dft_opt\n"
            "  beckmann-nbo recover --mol 020\n"
            "  beckmann-nbo report --mol 020 --out /tmp/report_020 --advanced\n"
        ),
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Print every command that would run, but do not execute it.")
    parser.add_argument("--mol", metavar="ID",
                        help="Target a single molecule/query by id.")
    parser.add_argument("--dir", metavar="PATH", default=None,
                        help="Local job directory to operate on (default: data/output/dft_opt).")

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    p_init = sub.add_parser("init", help="Write .env with your cluster's SSH/Gaussian/NBO7 settings")
    p_init.set_defaults(func=cmd_init)

    p_verify = sub.add_parser("verify", help="Preflight-check SSH connectivity, G16_PATH, and NBO7 setup")
    p_verify.set_defaults(func=cmd_verify)

    p_predict = sub.add_parser("predict", help="SMILES (or --csv) -> generate + submit Stage 1+2 DFT jobs")
    p_predict.add_argument("--smiles", help="A single ketone SMILES string")
    p_predict.add_argument("--name", help="Short id for --smiles (default: 'query')")
    p_predict.add_argument("--csv", help="CSV with 'id'/'SMILES' columns (same shape as data/input/benchmark.csv)")
    p_predict.add_argument("--continue", dest="continue_name", metavar="ID",
                            help="Generate + submit Stage 3 for a molecule whose Stage 1 already completed")
    p_predict.add_argument("--workdir", help="Working directory for generated files (default: data/output/query_predictions)")
    p_predict.set_defaults(func=cmd_predict)

    p_status = sub.add_parser("status", help="Per-stage job status, recovery caveats, and a live wCNmax prediction")
    p_status.add_argument("--no-download", action="store_true",
                          help="Skip downloading fresh logs from the cluster before checking status")
    p_status.set_defaults(func=cmd_status_command)

    p_recover = sub.add_parser("recover", help="Run one pass of the automated oscillation-recovery ladder")
    p_recover.set_defaults(func=cmd_recover)

    p_report = sub.add_parser("report", help="Write wCNmax/bond-order/E2PERT plots + a classical-vs-wCNmax comparison")
    p_report.add_argument("--out", required=True, help="Output directory for plots and classical_vs_wcnmax.txt")
    p_report.add_argument("--advanced", action="store_true", help="Also write the E2PERT-vs-R(N-O) plot")
    p_report.set_defaults(func=cmd_report)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
