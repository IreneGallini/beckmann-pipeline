"""
Summarize and plot the channel-resolved descriptors (Psi, Lambda, wCNmax,
w17max, w78max) across the N-O scan for the benchmark set, for discussion
with a PI -- not part of the regular prediction pipeline.

Produces:
  data/output/analysis/plots/{descriptor}.png  -- one plot per descriptor,
      split into R-outcome / F-outcome side-by-side panels (not one shared
      axes -- keeps each panel's line count manageable as the substrate
      count grows toward 34), R(N-O) on the x-axis, direct-labeled by
      substrate ID.
  data/output/analysis/plots/wcnmax_grid.png    -- all substrates' wCNmax(R)
      in one small-multiples figure (see beckmann.dft.viz.plot_wcnmax_grid),
      the 'all N in one graph' comparison view.
  data/output/analysis/descriptor_summary.md   -- condensed table: d/dR for
      each descriptor per substrate, plus whether wCNmax shows an interior
      extremum (the paper's central experimental/computational signature).

All data-loading logic (load_series(), find_wcnmax_extremum()) and the
reusable per-molecule wCNmax figure (plot_wcnmax_single/grid) live in the
beckmann package (beckmann.dft.descriptors / beckmann.dft.viz) -- this script
is a thin caller, per CLAUDE.md's "logic lives in the package" convention.
"""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

from beckmann_core.wcnmax_rule import find_wcnmax_extremum
from beckmann_nbo.config import DATA_INPUT, DATA_OUTPUT
from beckmann_nbo.descriptors import DESCRIPTORS, load_series
from beckmann_nbo.parse_cmo import classify_crossing

from viz import OUTCOME_COLOR, plot_wcnmax_grid

ANALYSIS_DIR = DATA_OUTPUT / "analysis"
PLOTS_DIR    = ANALYSIS_DIR / "plots"

LABELS = {
    "psi": "Ψ (Hyperconjugative Competition)",
    "log_lambda": "log₁₀(Λ)  (Frontier Dominance)",
    "wcnmax": "wCNmax  (nitrilium/CN channel)",
    "w17max": "w17max  (rearrangement channel)",
    "w78max": "w78max  (fragmentation channel)",
}

# Direct endpoint-labeling reads fine up to about this many lines on one
# axes -- past it, the labels themselves become the crowding problem, so
# fall back to letting the (still 2-entry) legend carry identity instead.
MAX_DIRECT_LABELS = 10


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _pct_of_outcome(meta: dict, outcome: str) -> float | None:
    """% of whichever product the experimental outcome reports (pct_A for
    'R', pct_B for 'F') -- the 'how decisive was this call' number used to
    shade lines in plot_wcnmax_grid()."""
    key = "pct_A" if outcome == "R" else "pct_B"
    val = meta.get(key)
    return float(val) if val not in (None, "", "None") else None


def _plot_descriptor_rf_split(descriptor: str, mols: list[str], per_mol_series: dict, outcomes: dict) -> None:
    """One descriptor, R-outcome and F-outcome side by side -- halves the
    line count per axes versus one shared plot, and keeps the R/F contrast
    (the actual scientific question) as the organizing structure rather
    than an afterthought color split."""
    fig, (ax_r, ax_f) = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    mols_by_outcome = {"R": [], "F": []}
    for mol in mols:
        mol_id = mol.split("_")[1]
        outcome = outcomes[f"mol_{mol_id}"]["exp_outcome"]
        mols_by_outcome.setdefault(outcome, []).append(mol)

    for ax, outcome in [(ax_r, "R"), (ax_f, "F")]:
        group = mols_by_outcome.get(outcome, [])
        color = OUTCOME_COLOR.get(outcome, "gray")
        label_all = len(group) <= MAX_DIRECT_LABELS
        for mol in group:
            r_values, y_by_descriptor = per_mol_series[mol]
            pts = [(r, y) for r, y in zip(r_values, y_by_descriptor[descriptor]) if y is not None]
            if not pts:
                continue
            pts.sort()
            xs, ys = zip(*pts)
            mol_id = mol.split("_")[1]
            ax.plot(xs, ys, marker="o", markersize=6, linewidth=1.6, color=color)
            if label_all:
                ax.annotate(mol_id, (xs[-1], ys[-1]), xytext=(6, 0),
                            textcoords="offset points", fontsize=8, color="dimgray", va="center")
        if not label_all:
            # Too many lines for direct labels to stay legible -- fall back
            # to a plain count in the panel title instead of cramming text.
            ax.text(0.02, 0.98, f"{len(group)} substrates", transform=ax.transAxes,
                     fontsize=8, color="dimgray", va="top")
        ax.set_xlabel("R(N-O)  (Å)")
        ax.set_title(f"Outcome = {outcome}", fontsize=10)

    ax_r.set_ylabel(LABELS[descriptor])
    fig.suptitle(f"{LABELS[descriptor]} vs. N-O distance", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out_path = PLOTS_DIR / f"{descriptor}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"-- wrote {out_path}")


def main() -> None:
    channel_rows    = _read_csv(ANALYSIS_DIR / "channel_descriptors.csv")
    slopes_rows     = {row["mol"]: row for row in _read_csv(ANALYSIS_DIR / "descriptor_slopes.csv")}
    extraction_rows = _read_csv(ANALYSIS_DIR / "cmo_channel_extraction.csv")
    outcomes        = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    mols = sorted({row["mol"] for row in channel_rows})
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    per_mol_series = {}
    extrema = {}
    crossings = {}
    outcome_by_mol = {}
    pct_by_mol = {}
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    for mol in mols:
        per_mol_series[mol] = load_series(mol, channel_rows)
        extrema[mol] = find_wcnmax_extremum(mol, extraction_rows)
        crossings[mol] = classify_crossing(mol, extraction_rows, dft_opt_dir / mol)
        mol_id = mol.split("_")[1]
        meta = outcomes[f"mol_{mol_id}"]
        outcome_by_mol[mol] = meta["exp_outcome"]
        pct_by_mol[mol] = _pct_of_outcome(meta, meta["exp_outcome"])

    # ---- per-descriptor exploration plots (R/F split) ----
    for descriptor in DESCRIPTORS:
        _plot_descriptor_rf_split(descriptor, mols, per_mol_series, outcomes)

    # ---- all-substrate wCNmax comparison grid ----
    wcnmax_series = {mol: (per_mol_series[mol][0], per_mol_series[mol][1]["wcnmax"]) for mol in mols}
    grid_fig = plot_wcnmax_grid(wcnmax_series, extrema, outcome_by_mol, pct_by_mol)
    grid_path = PLOTS_DIR / "wcnmax_grid.png"
    grid_fig.savefig(grid_path, dpi=150)
    plt.close(grid_fig)
    print(f"-- wrote {grid_path}")

    # ---- summary table ----
    lines = [
        "| mol | exp | d(Ψ)/dR | d(log₁₀Λ)/dR | d(wCNmax)/dR | d(w17max)/dR | d(w78max)/dR | "
        "wCNmax extremum | dip depth | crossing classification |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for mol in mols:
        mol_id = mol.split("_")[1]
        outcome = outcomes[f"mol_{mol_id}"]["exp_outcome"]
        slopes = slopes_rows[mol]
        extremum_info = extrema.get(mol)
        if extremum_info is None:
            extremum, depth_str = "no", "n/a"
        else:
            extremum = (
                f"yes @ R={extremum_info['R_star']:.4f} (MO {extremum_info['MO_index']}, "
                f"epsilon={extremum_info['epsilon_i_star']:.4f} a.u.)"
            )
            depth_str = f"{extremum_info['depth']:.4f}"

        crossing = crossings.get(mol) or {}
        crossing_label = crossing.get("label", "n/a")
        if crossing.get("reason"):
            crossing_str = f"{crossing_label} ({crossing['reason']})"
        else:
            crossing_str = crossing_label

        def fmt(key):
            val = slopes.get(key)
            return f"{float(val):.3f}" if val not in (None, "", "None") else "n/a"

        lines.append(
            f"| {mol} | {outcome} | {fmt('d_psi_dR')} | {fmt('d_log_lambda_dR')} | "
            f"{fmt('d_wcnmax_dR')} | {fmt('d_w17max_dR')} | {fmt('d_w78max_dR')} | {extremum} | "
            f"{depth_str} | {crossing_str} |"
        )

    table = "\n".join(lines)
    summary_path = ANALYSIS_DIR / "descriptor_summary.md"
    mol_ids = ", ".join(f"mol_{mol.split('_')[1]}" for mol in mols)
    summary_path.write_text(
        f"# Descriptor summary ({mol_ids})\n\n"
        "d/dR = least-squares slope over each molecule's N-O scan series. "
        "'wCNmax extremum' = interior local min/max in the wCNmax(R) series "
        "(the paper's central signature -- Table 2 reports this only for the one "
        "rearranging reference compound, none of the three fragmenting ones). "
        "'crossing classification' = beckmann.dft.parse_cmo.classify_crossing()'s "
        "verdict on whether the wCNmax MO handoff is a CONFIRMED avoided crossing "
        "(bracketed narrow eigenvalue gap between the identity-tracked pre/post-handoff "
        "MO pair, AND roughly conserved CN weight across them) vs. an unconfirmed "
        "handoff vs. no handoff at all. Across all 34 molecules the crossing partner is "
        "consistently the N-O sigma*/sigma antibond, not the aryl-migrating C-C "
        "antibond (data/output/analysis/cn_crossing_report.csv), so no aryl-coefficient "
        "swap is checked.\n\n"
        + table + "\n"
    )
    print(f"\n-- wrote {summary_path}\n")
    print(table)


if __name__ == "__main__":
    main()
