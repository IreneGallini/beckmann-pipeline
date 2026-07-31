"""
Plot mol_006_E's wCNmax(R) trend at 0.05 Å scan resolution -- the finer scan
that revealed a real interior minimum the standard 0.1 Å 5-point grid steps
directly over (see Notes.md, "mol_006_E's missing wCNmax minimum: resolved").

Standalone, not part of the regular pipeline. Reuses parse_cmo.parse_log()
directly against the completed fine-scan log rather than duplicating its
table-finding logic.

Output: data/output/analysis/plots/mol006_finescan_wcnmax.png
"""
import matplotlib.pyplot as plt

from beckmann_nbo import parse_cmo
from beckmann_nbo.config import DATA_OUTPUT
from beckmann_nbo.descriptors import get_substituent_map
from beckmann_nbo.scan import oxime_atom_map_from_gjf

MOL = "mol_006_E"
FINESCAN_LOG = DATA_OUTPUT / "dft_opt_finescan" / f"{MOL}_finescan" / f"{MOL}_finescan_scan.log"
OLD_MOL_DIR  = DATA_OUTPUT / "dft_opt" / MOL

# R(N-O) values that the old 0.1 A / 5-point grid actually sampled -- used to
# mark, on this finer plot, exactly which points the standard scan would have
# seen (everything else here is invisible at standard resolution).
OLD_GRID_STEP = 0.1


def load_series() -> tuple[list[float], list[float]]:
    ci, ni, oi, _ = oxime_atom_map_from_gjf(OLD_MOL_DIR / f"{MOL}_opt.gjf")
    subst = get_substituent_map(MOL, OLD_MOL_DIR)
    c_aryl, c_alkyl = subst["c_aryl"], subst["c_alkyl"]

    by_r: dict[float, dict] = {}
    for row in parse_cmo.parse_log(FINESCAN_LOG, ci, ni, oi, c_aryl, c_alkyl):
        by_r[round(row["r_no"], 4)] = row  # last entry per R wins (post-opt, not the Stable=Opt seed)

    r0_row = parse_cmo.parse_log(OLD_MOL_DIR / f"{MOL}_nbo.log", ci, ni, oi, c_aryl, c_alkyl)[0]
    by_r[round(r0_row["r_no"], 4)] = r0_row

    xs = sorted(by_r)
    ys = [by_r[x]["wcnmax"] for x in xs]
    return xs, ys


def main() -> None:
    xs, ys = load_series()
    r0 = xs[0]
    old_grid_xs = [x for x in xs if abs((x - r0) / OLD_GRID_STEP - round((x - r0) / OLD_GRID_STEP)) < 0.01]
    old_grid_ys = [ys[xs.index(x)] for x in old_grid_xs]

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(xs, ys, marker="o", markersize=7, linewidth=2, color="tab:green")
    ax.scatter(
        old_grid_xs, old_grid_ys, marker="o", s=170, facecolors="none",
        edgecolors="dimgray", linewidths=1.5, zorder=3,
    )

    ax.set_xlabel("R(N-O)  (Å)")
    ax.set_ylabel("wCNmax  (nitrilium/CN channel)")
    ax.set_title(f"{MOL}: wCNmax vs. N-O distance at 0.05 Å resolution")
    fig.tight_layout()

    out_path = DATA_OUTPUT / "analysis" / "plots" / "mol006_finescan_wcnmax.png"
    fig.savefig(out_path, dpi=150)
    print(f"-- wrote {out_path}")


if __name__ == "__main__":
    main()
