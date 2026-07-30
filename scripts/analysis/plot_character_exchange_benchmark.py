"""
Compare the diabatic character-exchange pattern across the full benchmark
set, colored by experimental R/F outcome. Not part of the regular prediction
pipeline -- for supervisor discussion, reading
data/output/analysis/character_exchange_benchmark.csv and
character_exchange_benchmark_summary.csv
(scripts/analysis/character_exchange_benchmark.py). Deliberately does NOT
compute or suggest a classification threshold -- these plots exist to make
the pattern visible, not to decide a cutoff. Reuses beckmann.dft.viz's
established R/F color convention (OUTCOME_COLOR) rather than introducing new
colors -- consistent with every other R/F-colored plot in this project.

Separate file from scripts/analysis/plot_character_exchange.py (the
reference-case script) -- does not import or touch it, or its outputs.

Produces:
  data/output/analysis/plots/character_exchange_benchmark_overlay.png -- one
      line per molecule, f_CN_CC vs delta_R, colored by R/F
  data/output/analysis/plots/character_exchange_benchmark_delta.png -- one
      point per molecule, delta on the y-axis, jittered by R/F on the x-axis
"""
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from beckmann.config import DATA_OUTPUT
from beckmann.dft.viz import OUTCOME_COLOR

ANALYSIS_DIR = DATA_OUTPUT / "analysis"
PLOTS_DIR = ANALYSIS_DIR / "plots"

# Fixed x-position per outcome for the delta strip plot, jittered around this.
OUTCOME_X = {"R": 0, "F": 1}


def _read_csv(path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def plot_overlay(detail_rows: list[dict]) -> None:
    by_mol: dict[str, list[dict]] = defaultdict(list)
    for row in detail_rows:
        by_mol[row["mol"]].append(row)

    fig, ax = plt.subplots(figsize=(8, 6))
    seen_outcomes = set()
    for mol, rows in sorted(by_mol.items()):
        rows.sort(key=lambda r: float(r["delta_R"]))
        outcome = rows[0]["exp_outcome"]
        x = [float(r["delta_R"]) for r in rows]
        y = [float(r["f_CN_CC"]) for r in rows]
        label = {"R": "Rearrangement (R)", "F": "Fragmentation (F)"}[outcome] if outcome not in seen_outcomes else None
        seen_outcomes.add(outcome)
        ax.plot(x, y, color=OUTCOME_COLOR[outcome], linewidth=1.3, alpha=0.7, marker="o", markersize=3, label=label)

    ax.set_xlabel("delta R(N-O)  (Å, relative to each molecule's own first scan point)")
    ax.set_ylabel("f$_{CN}$ of max-w$_{CC}$ MO")
    ax.set_title("Diabatic character exchange across the benchmark set, by experimental outcome")
    ax.grid(True, color="#e1e0d9")
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()

    out_path = PLOTS_DIR / "character_exchange_benchmark_overlay.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"-- wrote {out_path}")


def plot_delta(summary_rows: list[dict]) -> None:
    rng = np.random.default_rng(0)

    fig, ax = plt.subplots(figsize=(5, 6))
    for outcome, x0 in OUTCOME_X.items():
        deltas = [float(r["delta"]) for r in summary_rows if r["exp_outcome"] == outcome]
        if not deltas:
            continue
        jitter = rng.uniform(-0.12, 0.12, size=len(deltas))
        label = {"R": "Rearrangement (R)", "F": "Fragmentation (F)"}[outcome]
        ax.scatter(
            [x0] * len(deltas) + jitter, deltas, color=OUTCOME_COLOR[outcome],
            s=40, alpha=0.75, edgecolors="none", label=f"{label} (n={len(deltas)})",
        )

    ax.set_xticks(list(OUTCOME_X.values()))
    ax.set_xticklabels(["R", "F"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel(r"$\Delta$ = f$_{CN,start}$ - f$_{CN,end}$")
    ax.set_title("Character-exchange magnitude by experimental outcome")
    ax.grid(True, axis="y", color="#e1e0d9")
    ax.set_axisbelow(True)
    ax.legend(loc="upper center", frameon=False, bbox_to_anchor=(0.5, -0.08), ncol=2)
    fig.tight_layout()

    out_path = PLOTS_DIR / "character_exchange_benchmark_delta.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"-- wrote {out_path}")


def main() -> None:
    detail_rows = _read_csv(ANALYSIS_DIR / "character_exchange_benchmark.csv")
    summary_rows = _read_csv(ANALYSIS_DIR / "character_exchange_benchmark_summary.csv")
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    plot_overlay(detail_rows)
    plot_delta(summary_rows)


if __name__ == "__main__":
    main()
