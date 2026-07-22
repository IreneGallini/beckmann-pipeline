"""
Plot mol_006_E's wCNmax(R) trend from both methods: the trusted NBO7 series
(data/output/analysis/cmo_channel_extraction.csv) and the open-source
per-atom-pair series (data/output/analysis/wcnmax_channel_extraction_opensource.csv,
beckmann_alt.pair_nbo) -- both at the standard 0.05 A scan resolution, 7
directly-computed points each (no interpolated geometries -- see
Notes_open_source_alt.md for why an earlier interpolated-point investigation
was tried and removed). mol_006_E is the one molecule in the 6-molecule
test-set comparison where the two methods disagree on whether a genuine
interior wCNmax minimum exists (see Notes_open_source_alt.md /
wcnmax_rule_results_opensource.csv) -- NBO7 finds one (depth 0.1010,
predicting rearrangement, matching experiment); the open-source series does
not.

Output: data/output/analysis/plots/mol006_opensource_vs_nbo7_wcnmax.png
"""
import csv

import matplotlib.pyplot as plt

from beckmann.config import DATA_OUTPUT

MOL = "mol_006_E"

# Okabe-Ito colorblind-safe pair.
COLOR_TRUSTED = "#0072B2"      # blue
COLOR_OPEN_SOURCE = "#E69F00"  # orange


def load_cn_series(path, mol: str) -> tuple[list[float], list[float]]:
    by_r = {}
    for row in csv.DictReader(open(path)):
        if row["mol"] == mol and row["channel"] == "cn" and row["weight"] not in (None, "", "None"):
            by_r[round(float(row["R_NO"]), 4)] = float(row["weight"])
    xs = sorted(by_r)
    ys = [by_r[x] for x in xs]
    return xs, ys


def main() -> None:
    analysis_dir = DATA_OUTPUT / "analysis"
    trusted_xs, trusted_ys = load_cn_series(analysis_dir / "cmo_channel_extraction.csv", MOL)
    open_xs, open_ys = load_cn_series(analysis_dir / "wcnmax_channel_extraction_opensource.csv", MOL)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        trusted_xs, trusted_ys, marker="o", markersize=7, linewidth=2,
        color=COLOR_TRUSTED, label="NBO7 (trusted)",
    )
    ax.plot(
        open_xs, open_ys, marker="o", markersize=7, linewidth=2,
        color=COLOR_OPEN_SOURCE, label="Open-source (per-atom-pair)",
    )

    ax.set_xlabel("R(N-O)  (Å)")
    ax.set_ylabel("wCNmax  (nitrilium/CN channel)")
    ax.set_title(f"{MOL}: wCNmax vs. N-O distance -- NBO7 vs. open-source method")
    ax.legend(frameon=False)
    fig.tight_layout()

    out_path = analysis_dir / "plots" / "mol006_opensource_vs_nbo7_wcnmax.png"
    fig.savefig(out_path, dpi=150)
    print(f"-- wrote {out_path}")


if __name__ == "__main__":
    main()
