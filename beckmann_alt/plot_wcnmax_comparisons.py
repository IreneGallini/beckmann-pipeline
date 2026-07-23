"""
Plot each covered molecule's wCNmax(R) trend from both methods: the trusted NBO7
series (data/output/analysis/cmo_channel_extraction.csv) and the PySCF
per-atom-pair series (data/output/analysis/wcnmax_channel_extraction_pyscf.csv,
beckmann_alt.pair_nbo). One PNG per molecule present in the PySCF extraction
CSV (32 of 34 -- mol_005_E/mol_032_E permanently skipped, bromine, see
Notes_pyscf_alt.md). Originally written for mol_006_E only (the one molecule in
the 6-molecule test-set comparison where the two methods disagreed on whether a
genuine interior wCNmax minimum exists -- NBO7 finds one, depth 0.1010, predicting
rearrangement, matching experiment; the PySCF series didn't), generalized here
to the full covered set.

Output: data/output/analysis/plots/mol{NNN}_pyscf_vs_nbo7_wcnmax.png
"""
import csv

import matplotlib.pyplot as plt

from beckmann.config import DATA_OUTPUT

# Okabe-Ito colorblind-safe pair.
COLOR_TRUSTED = "#0072B2"      # blue
COLOR_PYSCF = "#E69F00"  # orange


def load_cn_series(path, mol: str) -> tuple[list[float], list[float]]:
    by_r = {}
    for row in csv.DictReader(open(path)):
        if row["mol"] == mol and row["channel"] == "cn" and row["weight"] not in (None, "", "None"):
            by_r[round(float(row["R_NO"]), 4)] = float(row["weight"])
    xs = sorted(by_r)
    ys = [by_r[x] for x in xs]
    return xs, ys


def plot_comparison(mol: str, analysis_dir, plots_dir) -> None:
    trusted_xs, trusted_ys = load_cn_series(analysis_dir / "cmo_channel_extraction.csv", mol)
    pyscf_xs, pyscf_ys = load_cn_series(analysis_dir / "wcnmax_channel_extraction_pyscf.csv", mol)

    fig, ax = plt.subplots(figsize=(7, 5))

    if trusted_xs:
        ax.plot(
            trusted_xs, trusted_ys, marker="o", markersize=7, linewidth=2,
            color=COLOR_TRUSTED, label="NBO7 (trusted)",
        )
    else:
        print(f"-- {mol}: no trusted NBO7 'cn' data, plotting PySCF only")
    ax.plot(
        pyscf_xs, pyscf_ys, marker="o", markersize=7, linewidth=2,
        color=COLOR_PYSCF, label="PySCF (per-atom-pair)",
    )

    ax.set_xlabel("R(N-O)  (Å)")
    ax.set_ylabel("wCNmax  (nitrilium/CN channel)")
    ax.set_title(f"{mol}: wCNmax vs. N-O distance -- NBO7 vs. PySCF method")
    ax.legend(frameon=False)
    fig.tight_layout()

    mol_id = mol.split("_")[1]
    out_path = plots_dir / f"mol{mol_id}_pyscf_vs_nbo7_wcnmax.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"-- wrote {out_path}")


def main() -> None:
    analysis_dir = DATA_OUTPUT / "analysis"
    plots_dir = analysis_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    mols = sorted({
        row["mol"] for row in csv.DictReader(open(analysis_dir / "wcnmax_channel_extraction_pyscf.csv"))
    })
    for mol in mols:
        plot_comparison(mol, analysis_dir, plots_dir)
    print(f"\n{len(mols)} molecules plotted")


if __name__ == "__main__":
    main()
