"""
Plot one molecule's wCNmax-vs-R(N-O) chart from already-parsed data.
Not part of the regular prediction pipeline -- a quick single-molecule
alternative to summarize_descriptors.py's full benchmark-set run, for
when you just want to look at one substrate's curve.

Requires `parse_cmo.py`/`descriptors.py` to have already been run, so
`channel_descriptors.csv` and `cmo_channel_extraction.csv` exist in
data/output/analysis/. Reuses the same reusable plotting helper
(viz.plot_wcnmax_single()) summarize_descriptors.py's per-molecule panels
are built from -- no separate plotting logic here.

Usage:
    PYTHONPATH=research python research/analysis_scripts/plot_single_wcnmax.py mol_002_E
    PYTHONPATH=research python research/analysis_scripts/plot_single_wcnmax.py mol_002_E --out mol_002_E_wcnmax.png
"""
import argparse
import csv
from pathlib import Path

from beckmann_core.wcnmax_rule import find_wcnmax_extremum
from beckmann_nbo.config import DATA_OUTPUT
from beckmann_nbo.descriptors import load_series

from viz import plot_wcnmax_single

ANALYSIS_DIR = DATA_OUTPUT / "analysis"


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mol", help="molecule name, e.g. mol_002_E")
    parser.add_argument("--out", type=Path, default=None,
                         help="output PNG path (default: <mol>_wcnmax.png in the current directory)")
    args = parser.parse_args()

    channel_rows = _read_csv(ANALYSIS_DIR / "channel_descriptors.csv")
    extraction_rows = _read_csv(ANALYSIS_DIR / "cmo_channel_extraction.csv")

    r_values, y_by_descriptor = load_series(args.mol, channel_rows)
    if not r_values:
        raise SystemExit(f"No rows found for '{args.mol}' in {ANALYSIS_DIR / 'channel_descriptors.csv'}")

    extremum = find_wcnmax_extremum(args.mol, extraction_rows)

    fig = plot_wcnmax_single(args.mol, r_values, y_by_descriptor["wcnmax"], extremum=extremum)
    out_path = args.out or Path(f"{args.mol}_wcnmax.png")
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
