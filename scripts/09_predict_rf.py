"""The payoff: R (rearrangement) vs F (fragmentation) prediction from the
wCNmax-minimum rule, plus the classical anti-periplanar dihedral baseline
for comparison (this project's core finding is that they often disagree
see ../README.md). Also writes a wCNmax-vs-R(N-O) plot. Needs
08_parse_descriptors.py to have already run.

Edit MOL_NAME below, then:
    python 09_predict_rf.py
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit.Chem import rdMolTransforms

from _common import QUERY_PREFIX, load_query_mol, local_substituent_map, sanitize_id, workdir_for

from beckmann_core.classical import get_oxime_atoms, predict as classical_predict
from beckmann_core.wcnmax_rule import find_wcnmax_minimum, predict_from_wcnmax, resolve_series
from beckmann_nbo.parse_cmo import collect_molecule as collect_cmo

MOL_NAME = "test1"


def plot_wcnmax(mol_name, channel_rows, minimum, out_path):
    by_stage = {
        r["stage"]: r for r in channel_rows
        if r["channel"] == "cn" and r["weight"] not in (None, "", "None")
    }
    rows = resolve_series(by_stage)
    pts = sorted((float(r["R_NO"]), float(r["weight"])) for r in rows)
    if not pts:
        return

    fig, ax = plt.subplots(figsize=(6, 4.2))
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color="#2a78d6", linewidth=2, marker="o", markersize=6)
    if minimum is not None:
        ax.plot(minimum["R_star"], minimum["w_star"], marker="o", markersize=10,
                 markerfacecolor="none", markeredgecolor="#e34948", markeredgewidth=2,
                 label="interior minimum")
        ax.legend(frameon=False)
    ax.set_xlabel("R(N-O) / Å")
    ax.set_ylabel("w$_{CN}^{max}$")
    ax.set_title(f"{mol_name}: wCNmax vs R(N-O)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="#d8d7d0", linewidth=0.6, alpha=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    workdir = workdir_for(mol_id)
    dft_opt_dir = workdir / "dft_opt"

    matches = sorted(dft_opt_dir.glob(f"{QUERY_PREFIX}_{mol_id}_*"))
    if not matches:
        print(f"ERROR: no directory matching {QUERY_PREFIX}_{mol_id}_* under {dft_opt_dir}", file=sys.stderr)
        sys.exit(1)
    mol_dir = matches[0]
    mol_name = mol_dir.name

    subst = local_substituent_map(mol_name, mol_dir, workdir)
    _, channel_rows = collect_cmo(mol_name, mol_dir, subst["c_aryl"], subst["c_alkyl"])

    minimum = find_wcnmax_minimum(mol_name, channel_rows)
    wcnmax_pred = predict_from_wcnmax(minimum)
    if minimum is not None:
        print(
            f"wCNmax prediction: {wcnmax_pred}  (interior minimum at "
            f"R(N-O)={minimum['R_star']:.3f} Å, w={minimum['w_star']:.4f}, "
            f"MO{minimum['MO_index']})"
        )
    else:
        print(f"wCNmax prediction: {wcnmax_pred}  (no interior wCNmax minimum found)")

    rdkit_mol = load_query_mol(mol_name, workdir)
    atom_ids = get_oxime_atoms(rdkit_mol)
    if atom_ids is None:
        classical_pred = "inspect"
        print("classical prediction: inspect (oxime atoms not identified)")
    else:
        cox_idx, nox_idx, oox_idx, c_aryl_idx, c_allyl_idx = atom_ids
        conf = rdkit_mol.GetConformer()
        d_aryl = abs(rdMolTransforms.GetDihedralDeg(conf, oox_idx, nox_idx, cox_idx, c_aryl_idx))
        d_allyl = abs(rdMolTransforms.GetDihedralDeg(conf, oox_idx, nox_idx, cox_idx, c_allyl_idx))
        classical_pred = classical_predict(d_aryl, d_allyl)
        print(
            f"classical prediction: {classical_pred}  "
            f"(dihedral_aryl={d_aryl:.1f}, dihedral_allyl={d_allyl:.1f})"
        )

    agreement = "unclear" if classical_pred == "inspect" else ("yes" if classical_pred == wcnmax_pred else "no")
    print(f"agreement: {agreement}")

    analysis_dir = workdir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    plot_path = analysis_dir / "wcnmax_vs_rno.png"
    plot_wcnmax(mol_name, channel_rows, minimum, plot_path)
    print(f"\nWrote {plot_path}")


if __name__ == "__main__":
    main()
