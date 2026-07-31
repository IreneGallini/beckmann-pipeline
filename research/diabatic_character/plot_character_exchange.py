"""
Plot the diabatic character-exchange pattern (beckmann.dft.diabatic_character)
on Tetiana's 4-point reference scan (example_scans/5_s1_Me.log .. 5_s4_Me.log).
Not part of the regular prediction pipeline -- a supervisor-facing figure
showing (a) the character-exchange headline result and (b) the concrete
evidence for why branch_tracking.py's E>0 candidate filter was rejected: the
true max-w_CN acceptor MO goes negative-energy at 3 of the 4 scan points.

Produces:
  data/output/analysis/character_exchange_reference.csv -- one row per scan
      point: R(N-O), the max-w_CN and max-w_CC carrier MOs, their canonical
      energies, and their f_CN/f_CC fractions.
  data/output/analysis/plots/character_exchange_reference.png -- two panels
      sharing the R(N-O) x-axis:
        top:    f_CN of the max-w_CC-carrier MO -- the character-exchange
                curve itself (N-side/mixed -> C-C-routed)
        bottom: canonical MO energies of both carriers, with a zero-energy
                reference line, showing the max-w_CN carrier crossing
                negative at s2-s4
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

from beckmann_nbo.config import DATA_OUTPUT

from diabatic_character.diabatic_character import track_diabatic_character

EXAMPLE_DIR = Path(__file__).resolve().parents[2] / "example_scans"
ANALYSIS_DIR = DATA_OUTPUT / "analysis"
PLOTS_DIR = ANALYSIS_DIR / "plots"

# R(N-O) per reference scan point -- Detailed_Orbital_Character_Exchange_
# Handout.docx Section 5 (same values validate_branch_tracking.py's
# GATE2_TARGETS uses).
POINTS = [("s1", 1.55), ("s2", 1.70), ("s3", 1.75), ("s4", 1.80)]

COLOR_CN = "#2a78d6"   # categorical slot 1 (blue) -- max-w_CN carrier
COLOR_CC = "#eb6834"   # categorical slot 2 (orange) -- max-w_CC carrier
COLOR_MUTED = "#898781"
COLOR_INK = "#0b0b0b"


def build_rows() -> list[dict]:
    logs = [EXAMPLE_DIR / f"5_{point}_Me.log" for point, _ in POINTS]
    results = track_diabatic_character(logs)
    rows = []
    for (point, r_no), r in zip(POINTS, results):
        rows.append({
            "point": point, "R_NO": r_no,
            "mo_CN": r["mo_CN"], "E_CN_hartree": r["E_CN"], "f_CN_CN": r["f_CN_CN"],
            "mo_CC": r["mo_CC"], "E_CC_hartree": r["E_CC"],
            "f_CC_CC": r["f_CC_CC"], "f_CN_CC": r["f_CN_CC"],
        })
    return rows


def write_csv(rows: list[dict]) -> Path:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = ANALYSIS_DIR / "character_exchange_reference.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def plot(rows: list[dict]) -> Path:
    r_values = [row["R_NO"] for row in rows]

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(7, 6.5), sharex=True, height_ratios=[1, 1.2],
    )

    # -- top panel: the character-exchange curve --
    f_cn_cc = [row["f_CN_CC"] for row in rows]
    ax_top.plot(r_values, f_cn_cc, marker="o", markersize=7, linewidth=2, color=COLOR_CN)
    for r, row in zip(r_values, rows):
        ax_top.annotate(
            f"MO{row['mo_CC']}", (r, row["f_CN_CC"]), textcoords="offset points",
            xytext=(0, 10), ha="center", fontsize=9, color=COLOR_INK,
        )
    ax_top.set_ylabel("f$_{CN}$ of max-w$_{CC}$ MO", fontsize=10)
    ax_top.set_ylim(-0.05, 1.1)
    ax_top.set_title("Diabatic character exchange: N-side/mixed → C-C-routed", fontsize=11, color=COLOR_INK)
    ax_top.grid(True, color="#e1e0d9")
    ax_top.set_axisbelow(True)

    # -- bottom panel: carrier energies, showing the E>0 filter would fail --
    e_cn = [row["E_CN_hartree"] for row in rows]
    e_cc = [row["E_CC_hartree"] for row in rows]
    ax_bot.axhline(0, color=COLOR_MUTED, linewidth=1, linestyle="--")
    ax_bot.plot(r_values, e_cn, marker="o", markersize=7, linewidth=2, color=COLOR_CN, label="max-w$_{CN}$ carrier")
    ax_bot.plot(r_values, e_cc, marker="o", markersize=7, linewidth=2, color=COLOR_CC, label="max-w$_{CC}$ carrier")
    for r, row in zip(r_values, rows):
        ax_bot.annotate(
            f"MO{row['mo_CN']}", (r, row["E_CN_hartree"]), textcoords="offset points",
            xytext=(0, -14), ha="center", fontsize=9, color=COLOR_INK,
        )
        ax_bot.annotate(
            f"MO{row['mo_CC']}", (r, row["E_CC_hartree"]), textcoords="offset points",
            xytext=(0, 8), ha="center", fontsize=9, color=COLOR_INK,
        )
    ax_bot.set_xlabel("R(N-O)  (Å)", fontsize=10)
    ax_bot.set_ylabel("Canonical MO energy (a.u.)", fontsize=10)
    ax_bot.grid(True, color="#e1e0d9")
    ax_bot.set_axisbelow(True)
    ax_bot.legend(loc="center right", fontsize=9, frameon=False)

    fig.tight_layout()
    out_path = PLOTS_DIR / "character_exchange_reference.png"
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    rows = build_rows()
    csv_path = write_csv(rows)
    print(f"-- wrote {csv_path}")
    png_path = plot(rows)
    print(f"-- wrote {png_path}")


if __name__ == "__main__":
    main()
