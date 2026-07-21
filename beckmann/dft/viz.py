"""
Reusable wCNmax visualization -- built as importable package code (not
script-local) so the per-molecule chart can eventually back the Flask web UI
as well as the analysis scripts, per CLAUDE.md's "logic lives in the package"
convention. See beckmann/dft/descriptors.py for the data side (load_series(),
find_wcnmax_minimum()) this module renders.

Past ~7 series, per-line color/legend stops scaling (dataviz skill's
series-count ladder) -- and the actual analytic question here, "is wCNmax
monotonic or does it have an interior minimum," is a per-molecule shape
question best answered on its own small axes, not by tangling many lines
into one shared plot. plot_wcnmax_grid() is the small-multiples answer to
"see all N substrates in one graph, as readable as possible."
"""
import math

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

OUTCOME_COLOR = {"R": "tab:green", "F": "tab:red"}


def _outcome_alpha(pct: float | None) -> float:
    """Map a 0-100 selectivity percentage (of whichever outcome the line's
    color already encodes) to a line alpha -- a decisive 95/5 split renders
    near-opaque, a borderline ~50/50 case renders faint, so two molecules
    sharing the same binary R/F color are still visually distinguishable by
    how decisive their experimental outcome actually was. None (no
    percentage available) renders fully opaque -- never guess a shade."""
    if pct is None:
        return 1.0
    # pct is "how much of the reported outcome" -- 50 is maximally ambiguous,
    # 100 is maximally decisive. Floor at 0.35 so faint lines stay visible.
    return 0.35 + 0.65 * (abs(pct - 50) / 50)


def plot_wcnmax_single(mol: str, r_values: list[float], y_values: list[float | None],
                        extremum: dict | None = None, outcome: str | None = None,
                        pct: float | None = None, ax: Axes | None = None) -> Figure:
    """One molecule's wCNmax(R) curve, with its interior minimum (if any)
    annotated as an open circle. Reusable: pass ax=None for a standalone
    figure (scripts, or a future Flask route rendering one substrate), or an
    existing ax to embed in a larger figure (plot_wcnmax_grid() below)."""
    owns_fig = ax is None
    if owns_fig:
        fig, ax = plt.subplots(figsize=(6, 4.5))
    else:
        fig = ax.figure

    color = OUTCOME_COLOR.get(outcome, "gray")
    alpha = _outcome_alpha(pct)
    pts = [(r, y) for r, y in zip(r_values, y_values) if y is not None]
    if pts:
        pts.sort()
        xs, ys = zip(*pts)
        ax.plot(xs, ys, marker="o", markersize=5 if owns_fig else 3,
                 linewidth=2 if owns_fig else 1.3, color=color, alpha=alpha)
    if extremum is not None:
        ax.scatter([extremum["R_star"]], [extremum["w_star"]],
                    s=80 if owns_fig else 30, facecolors="none",
                    edgecolors="black", linewidths=1.3, zorder=5)

    if owns_fig:
        ax.set_xlabel("R(N-O)  (Å)")
        ax.set_ylabel("wCNmax  (nitrilium/CN channel)")
        title = f"{mol}: wCNmax vs. N-O distance"
        if extremum is not None:
            title += f"  (interior minimum @ R={extremum['R_star']:.3f} Å)"
        ax.set_title(title, fontsize=10)
        fig.tight_layout()
    return fig


def plot_wcnmax_grid(per_mol_series: dict[str, tuple[list[float], list[float | None]]],
                      extrema: dict[str, dict | None],
                      outcomes: dict[str, str],
                      pcts: dict[str, float] | None = None) -> Figure:
    """All molecules in one figure, small-multiples grid: one tiny axes per
    molecule, shared y-axis scale for comparability, bordered by
    experimental R/F outcome, each built via plot_wcnmax_single(). This is
    the 'all N substrates in one graph' comparison view -- see module
    docstring for why small multiples rather than one shared axes."""
    pcts = pcts or {}
    mols = sorted(per_mol_series)
    n = len(mols)
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)

    all_y = [y for r_values, y_values in per_mol_series.values() for y in y_values if y is not None]
    y_min, y_max = (min(all_y), max(all_y)) if all_y else (0.0, 1.0)
    pad = (y_max - y_min) * 0.08 or 0.01

    fig, axes = plt.subplots(nrows, ncols, figsize=(2.3 * ncols, 2.1 * nrows), squeeze=False)
    for ax, mol in zip(axes.flat, mols):
        r_values, y_values = per_mol_series[mol]
        outcome = outcomes.get(mol)
        plot_wcnmax_single(mol, r_values, y_values, extrema.get(mol), outcome,
                            pcts.get(mol), ax=ax)
        ax.set_ylim(y_min - pad, y_max + pad)
        for spine in ax.spines.values():
            spine.set_edgecolor(OUTCOME_COLOR.get(outcome, "gray"))
            spine.set_linewidth(1.8)
        ax.set_xticks([])
        ax.set_yticks([])
        # Substrate ID direct-labeled per panel -- text wears a text token
        # (dimgray), never the series/border color, per the dataviz skill.
        ax.set_title(mol.split("_")[1], fontsize=8, color="dimgray", pad=2)

    for ax in axes.flat[n:]:
        ax.axis("off")

    fig.suptitle(
        "wCNmax(R) across all substrates  —  border color = experimental outcome "
        "(green = R, red = F), open circle = interior minimum",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig
