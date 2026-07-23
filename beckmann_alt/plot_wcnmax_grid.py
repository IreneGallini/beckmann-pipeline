"""
Small-multiples grid of the open-source PySCF wCNmax(R) series for all 32 covered
benchmark molecules (mol_005_E/mol_032_E permanently skipped -- bromine, unsupported
by this PySCF install's basis files, see Notes_open_source_alt.md), in the same visual
format as the trusted pipeline's data/output/analysis/plots/wcnmax_grid.png
(scripts/analysis/summarize_descriptors.py -> beckmann.dft.viz.plot_wcnmax_grid()).

Reuses plot_wcnmax_grid()/find_wcnmax_extremum() as-is -- both are already generic,
dict-in, no NBO7-specific logic. The only new piece is build_per_mol_series(), since
there's no open-source equivalent of channel_descriptors.csv (the descriptor-summary
CSV load_series() expects) -- wcnmax_channel_extraction_opensource.csv is extraction-
row shaped instead, so the series is built directly from it.

Output: data/output/analysis/plots/wcnmax_grid_opensource.png
"""
import csv
import json

import matplotlib.pyplot as plt

from beckmann.config import DATA_INPUT, DATA_OUTPUT
from beckmann.dft.descriptors import find_wcnmax_extremum
from beckmann.dft.viz import plot_wcnmax_grid


def build_per_mol_series(rows: list[dict]) -> dict[str, tuple[list[float], list[float]]]:
    by_mol: dict[str, dict[float, float]] = {}
    for row in rows:
        if row["channel"] != "cn" or row["weight"] in (None, "", "None"):
            continue
        by_mol.setdefault(row["mol"], {})[round(float(row["R_NO"]), 4)] = float(row["weight"])
    return {mol: (sorted(by_r), [by_r[r] for r in sorted(by_r)]) for mol, by_r in by_mol.items()}


def main() -> None:
    analysis_dir = DATA_OUTPUT / "analysis"
    rows = list(csv.DictReader(open(analysis_dir / "wcnmax_channel_extraction_opensource.csv")))
    outcomes = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    per_mol_series = build_per_mol_series(rows)
    mols = sorted(per_mol_series)

    extrema = {mol: find_wcnmax_extremum(mol, rows) for mol in mols}
    outcome_by_mol = {}
    pct_by_mol = {}
    for mol in mols:
        meta = outcomes[f"mol_{mol.split('_')[1]}"]
        outcome = meta["exp_outcome"]
        outcome_by_mol[mol] = outcome
        key = "pct_A" if outcome == "R" else "pct_B"
        val = meta.get(key)
        pct_by_mol[mol] = float(val) if val not in (None, "", "None") else None

    fig = plot_wcnmax_grid(per_mol_series, extrema, outcome_by_mol, pct_by_mol)
    plots_dir = analysis_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    out_path = plots_dir / "wcnmax_grid_opensource.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"-- wrote {out_path} ({len(mols)} molecules)")


if __name__ == "__main__":
    main()
