"""
wCNmax-vs-R(N-O) plot for one query molecule's scan series
(wcnmax_pyscf.run_scan_series() output, or pipeline.predict()'s
"wcnmax_series"). Visually consistent with (not code-shared with) the
sibling Gaussian/NBO7 package's own report plotting -- same colors/axis
style -- but this package must never import from that package (see
backend/tests/test_no_hpc_dependency.py), so the shared constants/helper
are copied by value, not imported.

Unlike that package's version, which reconciles rows from multiple
Gaussian log sources via a series-resolution step, this pipeline's series
is already one flat, already-sorted list per molecule (see
backend/tests/test_pipeline.py), so no reconciliation step is needed here.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SERIES_1 = "#2a78d6"  # blue


def _style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="#d8d7d0", linewidth=0.6, alpha=0.6)
    ax.set_axisbelow(True)


def plot_wcnmax(series: list[dict], minimum: dict | None, out_path: Path, name: str = "query") -> None:
    """series: rows shaped like run_scan_series()'s output / predict()'s
    "wcnmax_series" (must have "R_NO"/"weight" keys). minimum: find_wcnmax_
    minimum()'s result (R_star/w_star/...) or None. No-ops on an empty
    series, matching the sibling package's own report-plotting guard."""
    pts = sorted(
        (float(r["R_NO"]), float(r["weight"]))
        for r in series if r.get("weight") not in (None, "", "None")
    )
    if not pts:
        return

    fig, ax = plt.subplots(figsize=(6, 4.2))
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=SERIES_1, linewidth=2, marker="o", markersize=6)
    if minimum is not None:
        ax.plot(minimum["R_star"], minimum["w_star"], marker="o", markersize=10,
                 markerfacecolor="none", markeredgecolor="#e34948", markeredgewidth=2,
                 label="interior minimum")
        ax.legend(frameon=False)
    ax.set_xlabel("R(N-O) / Å")
    ax.set_ylabel("w$_{CN}^{max}$")
    ax.set_title(f"{name}: wCNmax vs R(N-O)")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
