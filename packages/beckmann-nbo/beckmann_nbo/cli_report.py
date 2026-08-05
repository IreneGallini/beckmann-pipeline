"""`beckmann-nbo report` -- per-molecule plots + a classical-vs-wCNmax
comparison, written into --out. Every number comes from the existing
parse_wiberg/parse_cmo/parse_nbo/descriptors/wcnmax_rule/classical
functions, called in-process (no CSV round-trip); only the plotting code
and the classical-vs-wCNmax comparison text are new.

Colors: slot 1 (blue #2a78d6) / slot 2 (orange #eb6834) from the project's
validated categorical palette, assigned by series identity (never by rank),
one y-axis per plot, legend present whenever a plot has >1 series.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit.Chem import rdMolTransforms

from beckmann_core.classical import get_oxime_atoms, predict as classical_predict
from beckmann_core.wcnmax_rule import find_wcnmax_minimum, predict_from_wcnmax, resolve_series
from beckmann_nbo.descriptors import _load_mols, compute_psi_row, get_substituent_map
from beckmann_nbo.hpc import DEFAULT_LOCAL_DFT_DIR, mol_dirs
from beckmann_nbo.inputs import STEP_SCAN_SOURCES
from beckmann_nbo.parse_cmo import collect_molecule as collect_cmo, collect_molecule_stepscan as collect_cmo_stepscan
from beckmann_nbo.parse_nbo import collect_molecule as collect_e2pert, collect_molecule_stepscan as collect_e2pert_stepscan
from beckmann_nbo.parse_wiberg import collect_molecule as collect_wiberg, collect_molecule_stepscan as collect_wiberg_stepscan

SERIES_1 = "#2a78d6"  # blue
SERIES_2 = "#eb6834"  # orange


def _abs_dihedral(conf, i: int, j: int, k: int, l: int) -> float:
    return abs(rdMolTransforms.GetDihedralDeg(conf, i, j, k, l))


def _style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="#d8d7d0", linewidth=0.6, alpha=0.6)
    ax.set_axisbelow(True)


def plot_wcnmax(mol: str, channel_rows: list[dict], minimum: dict | None, out_path: Path) -> None:
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
    ax.plot(xs, ys, color=SERIES_1, linewidth=2, marker="o", markersize=6)
    if minimum is not None:
        ax.plot(minimum["R_star"], minimum["w_star"], marker="o", markersize=10,
                 markerfacecolor="none", markeredgecolor="#e34948", markeredgewidth=2,
                 label="interior minimum")
        ax.legend(frameon=False)
    ax.set_xlabel("R(N-O) / Å")
    ax.set_ylabel("w$_{CN}^{max}$")
    ax.set_title(f"{mol}: wCNmax vs R(N-O)")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_wiberg(mol: str, wiberg_rows: list[dict], out_path: Path) -> None:
    by_stage = {r["point"]: r for r in wiberg_rows}
    rows = resolve_series(by_stage)
    pts = sorted((float(r["r_no"]), float(r["bond_order_aryl"]), float(r["bond_order_alkyl"])) for r in rows)
    if not pts:
        return

    fig, ax = plt.subplots(figsize=(6, 4.2))
    xs, aryl, alkyl = zip(*pts)
    ax.plot(xs, aryl, color=SERIES_1, linewidth=2, marker="o", markersize=6, label="C-C(aryl)")
    ax.plot(xs, alkyl, color=SERIES_2, linewidth=2, marker="o", markersize=6, label="C-C(alkyl)")
    ax.legend(frameon=False)
    ax.set_xlabel("R(N-O) / Å")
    ax.set_ylabel("Wiberg bond order")
    ax.set_title(f"{mol}: C-C bond order vs R(N-O)")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_e2pert(mol: str, e2pert_rows: list[dict], summary_rows: list[dict],
                 ci: int, ni: int, oi: int, c_aryl: int, c_alkyl: int, out_path: Path) -> None:
    r_no_by_stage = {r["stage"]: r["r_no"] for r in summary_rows}
    by_stage_e2: dict[str, list[dict]] = {}
    for row in e2pert_rows:
        by_stage_e2.setdefault(row["stage"], []).append(row)

    psi_by_stage = {}
    for stage, rows in by_stage_e2.items():
        if stage not in r_no_by_stage:
            continue
        psi = compute_psi_row(rows, ci, ni, oi, c_aryl, c_alkyl)
        psi_by_stage[stage] = {"r_no": r_no_by_stage[stage], **psi}

    resolved = resolve_series(psi_by_stage)
    pts = sorted((float(r["r_no"]), r["k_anti"], r["k_frag"]) for r in resolved)
    if not pts:
        return

    fig, ax = plt.subplots(figsize=(6, 4.2))
    xs, k_anti, k_frag = zip(*pts)
    ax.plot(xs, k_anti, color=SERIES_1, linewidth=2, marker="o", markersize=6, label="K_anti (rearrangement channel)")
    ax.plot(xs, k_frag, color=SERIES_2, linewidth=2, marker="o", markersize=6, label="K_frag (fragmentation channel)")
    ax.legend(frameon=False)
    ax.set_xlabel("R(N-O) / Å")
    ax.set_ylabel("E(2) sum / kcal mol$^{-1}$")
    ax.set_title(f"{mol}: hyperconjugative channels vs R(N-O)")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def classical_vs_wcnmax_line(mol: str, c_map: dict, wcnmax_pred: str) -> str:
    mols = _load_mols()
    rdkit_mol = mols.get(mol)
    if rdkit_mol is None:
        return f"{mol}: classical=unavailable (not in best_per_substrate.sdf)  wcnmax={wcnmax_pred}  agreement=unknown"

    atom_ids = get_oxime_atoms(rdkit_mol)
    if atom_ids is None:
        return f"{mol}: classical=inspect (oxime atoms not identified)  wcnmax={wcnmax_pred}  agreement=unclear"

    cox_idx, nox_idx, oox_idx, c_aryl_idx, c_allyl_idx = atom_ids
    conf = rdkit_mol.GetConformer()
    d_aryl = _abs_dihedral(conf, oox_idx, nox_idx, cox_idx, c_aryl_idx)
    d_allyl = _abs_dihedral(conf, oox_idx, nox_idx, cox_idx, c_allyl_idx)
    classical_pred = classical_predict(d_aryl, d_allyl)

    if classical_pred == "inspect":
        agreement = "unclear"
    else:
        agreement = "yes" if classical_pred == wcnmax_pred else "no"

    return (
        f"{mol}: classical={classical_pred} (dihedral_aryl={d_aryl:.1f}, dihedral_allyl={d_allyl:.1f})  "
        f"wcnmax={wcnmax_pred}  agreement={agreement}"
    )


def cmd_report(args) -> None:
    local_dir = Path(args.dir) if args.dir else DEFAULT_LOCAL_DFT_DIR
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    dirs = mol_dirs(local_dir, args.mol)
    if not dirs:
        print(f"ERROR: no molecule directories found under {local_dir}", file=sys.stderr)
        sys.exit(1)

    comparison_lines = []
    for mol_dir in dirs:
        mol = mol_dir.name
        try:
            c_map = get_substituent_map(mol, mol_dir)
        except ValueError as e:
            print(f"SKIP {mol}: {e}", file=sys.stderr)
            continue
        c_aryl, c_alkyl = c_map["c_aryl"], c_map["c_alkyl"]

        mol_out = out_dir / mol
        mol_out.mkdir(parents=True, exist_ok=True)

        if mol in STEP_SCAN_SOURCES:
            summary_rows, channel_rows = collect_cmo_stepscan(mol, mol_dir, c_aryl, c_alkyl)
            wiberg_rows = collect_wiberg_stepscan(mol, mol_dir, c_aryl, c_alkyl)
        else:
            summary_rows, channel_rows = collect_cmo(mol, mol_dir, c_aryl, c_alkyl)
            wiberg_rows = collect_wiberg(mol, mol_dir, c_aryl, c_alkyl)

        minimum = find_wcnmax_minimum(mol, channel_rows)
        prediction = predict_from_wcnmax(minimum)

        plot_wcnmax(mol, channel_rows, minimum, mol_out / "wcnmax_vs_rno.png")
        plot_wiberg(mol, wiberg_rows, mol_out / "wiberg_bond_order_vs_rno.png")

        if args.advanced:
            e2pert_rows = collect_e2pert_stepscan(mol, mol_dir) if mol in STEP_SCAN_SOURCES else collect_e2pert(mol, mol_dir)
            plot_e2pert(
                mol, e2pert_rows, summary_rows,
                c_map["ci"], c_map["ni"], c_map["oi"], c_aryl, c_alkyl,
                mol_out / "e2pert_vs_rno.png",
            )

        comparison_lines.append(classical_vs_wcnmax_line(mol, c_map, prediction))
        print(f"{mol}: wrote plots to {mol_out}")

    (out_dir / "classical_vs_wcnmax.txt").write_text("\n".join(comparison_lines) + "\n")
    print(f"\nWrote {out_dir / 'classical_vs_wcnmax.txt'}")
