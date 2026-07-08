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
from beckmann.dft.descriptors import SERIES_STAGES

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
    series = [by_stage[s] for s in SERIES_STAGES if s in by_stage]
    r_values = [float(row["r_no"]) for row in series]
    y_by_descriptor = {d: [_float_or_none(row[d]) for row in series] for d in DESCRIPTORS}
    return r_values, y_by_descriptor


def has_interior_extremum(r_values: list[float], y_values: list[float | None]) -> bool:
    """True if there's a local min or max strictly between the first and last point (paper's wCNmax signature)."""
    pts = [(r, y) for r, y in zip(r_values, y_values) if y is not None]
    if len(pts) < 3:
        return False
    pts.sort()
    for i in range(1, len(pts) - 1):
        y_prev, y_cur, y_next = pts[i - 1][1], pts[i][1], pts[i + 1][1]
        if (y_cur < y_prev and y_cur < y_next) or (y_cur > y_prev and y_cur > y_next):
            return True
    return False


def main() -> None:
    channel_rows = _read_csv(ANALYSIS_DIR / "channel_descriptors.csv")
    slopes_rows  = {row["mol"]: row for row in _read_csv(ANALYSIS_DIR / "descriptor_slopes.csv")}
    outcomes     = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    mols = sorted({row["mol"] for row in channel_rows})
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    per_mol_series = {}
    for mol in mols:
        r_values, y_by_descriptor = load_series(mol, channel_rows)
        per_mol_series[mol] = (r_values, y_by_descriptor)

    # ---- plots ----
    for descriptor in DESCRIPTORS:
        fig, ax = plt.subplots(figsize=(6, 4.5))
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
            ax.plot(xs, ys, marker="o", label=f"mol_{mol_id} ({outcome})",
                     color=OUTCOME_COLOR.get(outcome, "gray"))
        ax.set_xlabel("R(N-O)  (Å)")
        ax.set_ylabel(LABELS[descriptor])
        ax.set_title(f"{LABELS[descriptor]} vs. N-O distance")
        ax.legend()
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
        r_values, y_by_descriptor = per_mol_series[mol]
        extremum = "yes" if has_interior_extremum(r_values, y_by_descriptor["wcnmax"]) else "no"

        def fmt(key):
            val = slopes.get(key)
            return f"{float(val):.3f}" if val not in (None, "", "None") else "n/a"

        lines.append(
            f"| {mol} | {outcome} | {fmt('d_psi_dR')} | {fmt('d_log_lambda_dR')} | "
            f"{fmt('d_wcnmax_dR')} | {fmt('d_w17max_dR')} | {fmt('d_w78max_dR')} | {extremum} |"
        )

    table = "\n".join(lines)
    summary_path = ANALYSIS_DIR / "descriptor_summary.md"
    summary_path.write_text(
        "# Descriptor summary (mol_002, 006, 020, 021)\n\n"
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
