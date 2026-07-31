"""
Plot Wiberg bond index (NAO basis) for the two migrating C-C bonds --
central-C to aryl-C (rearrangement channel) and central-C to alkyl-C
(fragmentation channel) -- across the N-O scan, for every substrate with
bond-order data available. Not part of the regular prediction pipeline --
for supervisor discussion, alongside the existing wCNmax/Lambda/Psi plots
(scripts/analysis/summarize_descriptors.py).

Produces:
  data/output/analysis/plots/bond_order_scan.png -- one small-multiples grid
      (sized to however many molecules bond_order_scan.csv actually
      contains -- parse_wiberg.py already scopes that to ALL_IDS, see
      beckmann/dft/inputs.py), one panel per molecule, R(N-O) on the x-axis,
      two lines per panel (aryl/alkyl bond order)
"""
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

from beckmann_nbo.config import DATA_OUTPUT

ANALYSIS_DIR = DATA_OUTPUT / "analysis"
PLOTS_DIR    = ANALYSIS_DIR / "plots"

# Color encodes bond identity (aryl vs alkyl channel), not experimental
# outcome -- there's no R/F axis here, unlike summarize_descriptors.py's
# per-descriptor plots, so the same 2 colors are used consistently in every
# panel and a single shared legend covers the whole figure.
SERIES = {
    "bond_order_aryl":  ("C-C(aryl) bond order",  "tab:blue"),
    "bond_order_alkyl": ("C-C(alkyl) bond order", "tab:orange"),
}


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def load_series(mol: str, rows: list[dict]) -> list[tuple[float, dict]]:
    """(R, row) pairs for one molecule, sorted by R(N-O)."""
    pts = [(float(r["R"]), r) for r in rows if r["mol"] == mol and r["R"] not in (None, "", "None")]
    pts.sort(key=lambda t: t[0])
    return pts


def main() -> None:
    rows = _read_csv(ANALYSIS_DIR / "bond_order_scan.csv")
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    mols = sorted({row["mol"] for row in rows})
    ncols = math.ceil(math.sqrt(len(mols))) if mols else 1
    nrows = math.ceil(len(mols) / ncols) if mols else 1

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.8 * nrows), squeeze=False)
    for ax, mol in zip(axes.flat, mols):
        pts = load_series(mol, rows)
        if not pts:
            ax.set_visible(False)
            continue
        r_values = [r for r, _ in pts]
        for field, (label, color) in SERIES.items():
            y_values = [float(row[field]) for _, row in pts]
            ax.plot(r_values, y_values, marker="o", markersize=4, linewidth=1.5, color=color, label=label)
        ax.set_xlabel("R(N-O)  (Å)", fontsize=8)
        ax.set_ylabel("Wiberg bond index", fontsize=8)
        ax.set_title(mol.split("_")[1], fontsize=9, color="dimgray")
        ax.tick_params(labelsize=7)

    for ax in axes.flat[len(mols):]:
        ax.axis("off")

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("C-C bond order (Wiberg, NAO basis) vs. N-O distance", y=1.05)
    fig.tight_layout()

    out_path = PLOTS_DIR / "bond_order_scan.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"-- wrote {out_path}")


if __name__ == "__main__":
    main()
