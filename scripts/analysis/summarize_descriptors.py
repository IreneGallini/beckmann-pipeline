"""
Summarize and plot the channel-resolved descriptors (Psi, Lambda, wCNmax,
w17max, w78max) across the N-O scan for all 4 test molecules, for discussion
with a supervisor -- not part of the regular prediction pipeline.

Produces:
  data/output/analysis/plots/{descriptor}.png  -- one plot per descriptor,
      R(N-O) on the x-axis, one line per substrate, color-coded by
      experimental outcome (R = rearrangement, F = fragmentation)
  data/output/analysis/descriptor_summary.md   -- condensed table: d/dR for
      each descriptor per substrate, plus whether wCNmax shows an interior
      extremum (the paper's central experimental/computational signature)
"""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

from beckmann.config import DATA_INPUT, DATA_OUTPUT
from beckmann.dft.descriptors import resolve_series

ANALYSIS_DIR = DATA_OUTPUT / "analysis"
PLOTS_DIR    = ANALYSIS_DIR / "plots"

DESCRIPTORS = ["psi", "log_lambda", "wcnmax", "w17max", "w78max"]
LABELS = {
    "psi": "Ψ (Hyperconjugative Competition)",
    "log_lambda": "log₁₀(Λ)  (Frontier Dominance)",
    "wcnmax": "wCNmax  (nitrilium/CN channel)",
    "w17max": "w17max  (rearrangement channel)",
    "w78max": "w78max  (fragmentation channel)",
}
# Color still encodes experimental outcome (R vs F) rather than substrate identity --
# with dozens of substrates eventually, a distinct hue per line stops scaling long
# before a 2-color R/F split does. Individual lines are told apart by a direct label
# at the endpoint instead (see main()), not by color.
OUTCOME_COLOR = {"R": "tab:green", "F": "tab:red"}


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _float_or_none(v):
    if v in (None, "", "None"):
        return None
    return float(v)


def load_series(mol: str, channel_rows: list[dict]) -> tuple[list[float], dict[str, list[float | None]]]:
    """R(N-O) values and per-descriptor y-values for the 5-point series, in SERIES_STAGES order."""
    by_stage = {row["stage"]: row for row in channel_rows if row["mol"] == mol}
    series = resolve_series(by_stage)
    r_values = [float(row["r_no"]) for row in series]
    y_by_descriptor = {d: [_float_or_none(row[d]) for row in series] for d in DESCRIPTORS}
    return r_values, y_by_descriptor


def find_wcnmax_extremum(mol: str, extraction_rows: list[dict]) -> dict | None:
    """R_star/w_star/MO_index/epsilon_i_star at the interior wCNmax extremum, if any.

    MO_index/epsilon_i_star are backfilled from cmo_channel_extraction.csv's 'cn'
    channel rows (beckmann/dft/parse_cmo.py) rather than recomputed here -- that's
    the only place which virtual MO achieved the max weight, and its orbital energy,
    are actually recorded.
    """
    by_stage = {
        r["stage"]: r for r in extraction_rows
        if r["mol"] == mol and r["channel"] == "cn" and r["weight"] not in (None, "", "None")
    }
    rows = resolve_series(by_stage)
    pts = [(float(r["R_NO"]), float(r["weight"]), r["MO_index"], r["epsilon_i_star"]) for r in rows]
    if len(pts) < 3:
        return None
    pts.sort(key=lambda p: p[0])
    for i in range(1, len(pts) - 1):
        _, w_prev, _, _ = pts[i - 1]
        r_cur, w_cur, mo_cur, eps_cur = pts[i]
        _, w_next, _, _ = pts[i + 1]
        if (w_cur < w_prev and w_cur < w_next) or (w_cur > w_prev and w_cur > w_next):
            return {
                "R_star": r_cur, "w_star": w_cur, "MO_index": mo_cur,
                "epsilon_i_star": float(eps_cur) if eps_cur not in (None, "", "None") else None,
            }
    return None


def main() -> None:
    channel_rows    = _read_csv(ANALYSIS_DIR / "channel_descriptors.csv")
    slopes_rows     = {row["mol"]: row for row in _read_csv(ANALYSIS_DIR / "descriptor_slopes.csv")}
    extraction_rows = _read_csv(ANALYSIS_DIR / "cmo_channel_extraction.csv")
    outcomes        = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    mols = sorted({row["mol"] for row in channel_rows})
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    per_mol_series = {}
    for mol in mols:
        r_values, y_by_descriptor = load_series(mol, channel_rows)
        per_mol_series[mol] = (r_values, y_by_descriptor)

    # ---- plots ----
    for descriptor in DESCRIPTORS:
        fig, ax = plt.subplots(figsize=(6, 4.5))
        seen_outcomes = set()
        for mol in mols:
            r_values, y_by_descriptor = per_mol_series[mol]
            y_values = y_by_descriptor[descriptor]
            pts = [(r, y) for r, y in zip(r_values, y_values) if y is not None]
            if not pts:
                continue
            pts.sort()
            xs, ys = zip(*pts)
            mol_id = mol.split("_")[1]
            outcome = outcomes[f"mol_{mol_id}"]["exp_outcome"]
            color = OUTCOME_COLOR.get(outcome, "gray")
            # One legend entry per outcome (R/F), not per molecule -- individual
            # lines are identified by the direct label at their endpoint instead,
            # so the legend stays 2 entries regardless of substrate count.
            label = outcome if outcome not in seen_outcomes else None
            seen_outcomes.add(outcome)
            ax.plot(xs, ys, marker="o", markersize=7, linewidth=2, color=color, label=label)
            # Direct label at the line's endpoint -- text token color (not the
            # series color), per beckmann-dataviz: text never wears the data color.
            ax.annotate(
                mol_id, (xs[-1], ys[-1]), xytext=(6, 0), textcoords="offset points",
                fontsize=8, color="dimgray", va="center",
            )
        ax.set_xlabel("R(N-O)  (Å)")
        ax.set_ylabel(LABELS[descriptor])
        ax.set_title(f"{LABELS[descriptor]} vs. N-O distance")
        ax.legend(title="outcome")
        fig.tight_layout()
        out_path = PLOTS_DIR / f"{descriptor}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"-- wrote {out_path}")

    # ---- summary table ----
    lines = [
        "| mol | exp | d(Ψ)/dR | d(log₁₀Λ)/dR | d(wCNmax)/dR | d(w17max)/dR | d(w78max)/dR | wCNmax extremum |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for mol in mols:
        mol_id = mol.split("_")[1]
        outcome = outcomes[f"mol_{mol_id}"]["exp_outcome"]
        slopes = slopes_rows[mol]
        extremum_info = find_wcnmax_extremum(mol, extraction_rows)
        if extremum_info is None:
            extremum = "no"
        else:
            extremum = (
                f"yes @ R={extremum_info['R_star']:.4f} (MO {extremum_info['MO_index']}, "
                f"epsilon={extremum_info['epsilon_i_star']:.4f} a.u.)"
            )

        def fmt(key):
            val = slopes.get(key)
            return f"{float(val):.3f}" if val not in (None, "", "None") else "n/a"

        lines.append(
            f"| {mol} | {outcome} | {fmt('d_psi_dR')} | {fmt('d_log_lambda_dR')} | "
            f"{fmt('d_wcnmax_dR')} | {fmt('d_w17max_dR')} | {fmt('d_w78max_dR')} | {extremum} |"
        )

    table = "\n".join(lines)
    summary_path = ANALYSIS_DIR / "descriptor_summary.md"
    mol_ids = ", ".join(f"mol_{mol.split('_')[1]}" for mol in mols)
    summary_path.write_text(
        f"# Descriptor summary ({mol_ids})\n\n"
        "d/dR = least-squares slope over the 5-point N-O scan. "
        "'wCNmax extremum' = interior local min/max in the wCNmax(R) series "
        "(the paper's central signature -- Table 2 reports this only for the one "
        "rearranging reference compound, none of the three fragmenting ones).\n\n"
        + table + "\n"
    )
    print(f"\n-- wrote {summary_path}\n")
    print(table)


if __name__ == "__main__":
    main()
