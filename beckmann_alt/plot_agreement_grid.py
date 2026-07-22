"""
Visualize data/output/analysis/wcnmax_rule_results_opensource.csv -- for every
molecule tested so far (17: the original 6 test-set molecules + 11 more F-labeled
substrates run to check specificity, see Notes_open_source_alt.md), whether the
open-source wCNmax-minima rule and the trusted NBO7 rule each correctly predicted
the experimental R/F outcome.

Result this reproduces: both methods land at exactly 9/17 (53%) -- tied, but by
being wrong on different molecules (an even 2-2 split among the 4 molecules where
the two methods' predicted label actually differs), not by agreeing with each
other.

Output: data/output/analysis/plots/wcnmax_opensource_vs_nbo7_agreement.png
"""
import csv

import matplotlib.pyplot as plt

from beckmann.config import DATA_OUTPUT

# Okabe-Ito colorblind-safe pair. Same hues as the mol_006_E comparison plot, but
# here they encode CORRECT/INCORRECT (a different axis, not method identity) --
# marker shape is also varied so identity never depends on color alone.
COLOR_CORRECT = "#0072B2"    # blue
COLOR_INCORRECT = "#E69F00"  # orange


def load_rows(path):
    return list(csv.DictReader(open(path)))


def main() -> None:
    analysis_dir = DATA_OUTPUT / "analysis"
    rows = load_rows(analysis_dir / "wcnmax_rule_results_opensource.csv")

    # F-labeled first, then R-labeled; alphabetical within each group.
    rows.sort(key=lambda r: (r["exp_outcome"], r["mol"]))

    mols = [r["mol"] for r in rows]
    n = len(mols)
    x = list(range(n))

    os_correct = [r["pyscf_agreement"] == "yes" for r in rows]
    nbo7_correct = [r["nbo_agreement"] == "yes" for r in rows]

    fig, ax = plt.subplots(figsize=(max(9, n * 0.75), 3.6))

    # Shade the F/R groups so the class split is visible at a glance.
    n_f = sum(1 for r in rows if r["exp_outcome"] == "F")
    if n_f:
        ax.axvspan(-0.5, n_f - 0.5, color="0.94", zorder=0)
    if n_f < n:
        ax.axvspan(n_f - 0.5, n - 0.5, color="1.0", zorder=0)

    for row_y, label, correct in [(1, "NBO7 (trusted)", nbo7_correct), (0, "PySCF", os_correct)]:
        for xi, ok in zip(x, correct):
            color = COLOR_CORRECT if ok else COLOR_INCORRECT
            marker = "o" if ok else "x"
            if ok:
                ax.scatter(xi, row_y, marker=marker, s=140, color=color, zorder=3)
            else:
                ax.scatter(xi, row_y, marker=marker, s=110, color=color, linewidths=2.5, zorder=3)

    ax.set_yticks([0, 1])
    ax.set_yticklabels(["PySCF", "NBO7 (trusted)"])
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlim(-0.6, n - 0.4)

    xtick_labels = [f"{r['mol'].replace('mol_', '').replace('_E', '').replace('_Z', '')}\n({r['exp_outcome']})" for r in rows]
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels, fontsize=9)

    ax.set_title(
        "PySCF vs. NBO7 wCNmax-rule agreement with experiment, per molecule\n"
        "(both 9/17 overall -- tied, but wrong on different molecules)"
    )

    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_CORRECT,
                   markersize=10, label="Correctly predicted exp. outcome"),
        plt.Line2D([0], [0], marker="x", color=COLOR_INCORRECT, markersize=10,
                   markeredgewidth=2.5, linestyle="None", label="Incorrectly predicted"),
    ]
    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.32),
              ncol=2, frameon=False, fontsize=9)

    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(left=False)
    fig.tight_layout()

    out_path = analysis_dir / "plots" / "wcnmax_opensource_vs_nbo7_agreement.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"-- wrote {out_path}")


if __name__ == "__main__":
    main()
