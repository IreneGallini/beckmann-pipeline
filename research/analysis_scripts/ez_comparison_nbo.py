"""
Step 3d of the priority E/Z re-run: parse each of the 9 "other isomer"
Stage 3 scan logs (data/output/dft_opt_ez_other/{name}/{name}_scan.log) into
a live wCNmax R/F prediction, combine with the already-run isomer's NBO
prediction (from wcnmax_rule_results_pyscf.csv's NBO column) into one
comparison table. Any of the 9 not yet finished within the 24h window gets
'pending' rather than blocking the rest of the table.

Reuses beckmann_nbo's own in-process prediction path (the same one
`beckmann-nbo status` uses, see beckmann_nbo/cli_status.py:_predict) rather
than reimplementing scan-log parsing:
  beckmann_nbo.parse_cmo.collect_molecule() -> channel_rows
  beckmann_core.wcnmax_rule.find_wcnmax_minimum/predict_from_wcnmax()

get_substituent_map() (cli_status.py's usual source for c_aryl/c_alkyl) reads
best_per_substrate.sdf, which only has the already-run isomer per molecule --
the "other" isomer isn't in there. So c_aryl/c_alkyl are derived the same way
ez_comparison_pyscf.py/prepare_ez_other_isomer_dft.py already do: straight
from get_oxime_atoms() on the isomer's own AIMNet2-optimized mol in
best_aimnet_optimized.sdf (same underlying geometry/atom-ordering that was
carried into the Stage 1 .gjf, so atom indices line up), cross-checked
against the independently-derived (ci, ni, oi) parsed from the molecule's own
.gjf title line -- mirrors get_substituent_map()'s own cross-check.
"""
import csv
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem

from beckmann_core.classical import get_oxime_atoms
from beckmann_core.wcnmax_rule import find_wcnmax_minimum, predict_from_wcnmax
from beckmann_nbo.config import DATA_INPUT, DATA_OUTPUT
from beckmann_nbo.log_diagnostics import FailureCategory, classify_scan
from beckmann_nbo.parse_cmo import collect_molecule
from beckmann_nbo.scan import oxime_atom_map_from_gjf

from ez_energy_check import energy_table

BEST_AIMNET_SDF = DATA_OUTPUT / "aimnet_optimized" / "best_aimnet_optimized.sdf"
DFT_OPT_DIR = DATA_OUTPUT / "dft_opt_ez_other"
EXISTING_RESULTS_CSV = DATA_OUTPUT / "analysis" / "wcnmax_rule_results_pyscf.csv"
BENCHMARK_CSV = DATA_INPUT / "benchmark.csv"
OUT_CSV = DATA_OUTPUT / "analysis" / "ez_comparison_nbo.csv"

FIELDS = ["mol_id", "lowest_energy_isomer", "other_isomer",
          "prediction_lowest_energy", "prediction_other_isomer",
          "pct_product_A", "pct_product_B"]


def load_existing_predictions(csv_path=EXISTING_RESULTS_CSV) -> dict[str, str]:
    out = {}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            out[row["mol"]] = row["NBO"]
    return out


def load_product_pcts(csv_path=BENCHMARK_CSV) -> dict[str, tuple[float, float]]:
    """{'mol_001': (33.0, 67.0), ...} -- (% product A, % product B) per molecule id."""
    out = {}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            mol_id = f"mol_{int(row['id']):03d}"
            out[mol_id] = (float(row["% product A"]), float(row["% product B"]))
    return out


def substituent_indices(name: str, mol_dir: Path, mols_by_name: dict) -> tuple[int, int]:
    """(c_aryl, c_alkyl), 1-based -- from the isomer's own AIMNet2-optimized
    mol, cross-checked against the .gjf title line's (ci, ni, oi)."""
    mol = mols_by_name[name]
    atom_ids = get_oxime_atoms(mol)
    if atom_ids is None:
        raise ValueError(f"{name}: oxime substructure not found in AIMNet2 geometry")
    cox, nox, oox, c_aryl, c_allyl = (idx + 1 for idx in atom_ids)

    gjf_ci, gjf_ni, gjf_oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{name}_opt.gjf")
    if (cox, nox, oox) != (gjf_ci, gjf_ni, gjf_oi):
        raise ValueError(
            f"{name}: oxime atom map mismatch -- AIMNet2 sdf ({cox},{nox},{oox}) "
            f"vs .gjf title ({gjf_ci},{gjf_ni},{gjf_oi})"
        )
    return c_aryl, c_allyl


def nbo_prediction_for(name: str, mols_by_name: dict) -> str:
    """'R'/'F'/'pending' -- 'pending' if the Stage 3 scan log isn't a clean,
    Normal-terminated run yet."""
    mol_dir = DFT_OPT_DIR / name
    scan_log = mol_dir / f"{name}_scan.log"
    if not scan_log.exists():
        return "pending"

    diagnoses = classify_scan(mol_dir, name)
    scan_diag = next((d for d in diagnoses if d.stage == "scan"), None)
    if scan_diag is None or scan_diag.category != FailureCategory.NORMAL:
        return "pending"

    c_aryl, c_alkyl = substituent_indices(name, mol_dir, mols_by_name)
    _summary_rows, channel_rows = collect_molecule(name, mol_dir, c_aryl, c_alkyl)
    cn_points = {
        r["stage"] for r in channel_rows
        if r["channel"] == "cn" and r["weight"] not in (None, "", "None")
    }
    if len(cn_points) < 3:
        return "pending"  # not enough resolved points for a reliable extremum search yet

    minimum = find_wcnmax_minimum(name, channel_rows)
    return predict_from_wcnmax(minimum)


def main(mol_ids=None) -> None:
    from ez_energy_check import ALREADY_RUN
    rows = energy_table(mol_ids=mol_ids if mol_ids is not None else tuple(ALREADY_RUN))
    existing = load_existing_predictions()
    product_pcts = load_product_pcts()

    suppl = Chem.SDMolSupplier(str(BEST_AIMNET_SDF), removeHs=False)
    mols_by_name = {m.GetProp("_Name"): m for m in suppl if m is not None}

    out_rows = []
    pending = []
    for r in rows:
        mol_id = r["mol_id"]
        lowest, other = r["lowest_energy_isomer"], r["other_isomer"]
        already_run_isomer = r["already_run_isomer"]

        other_prediction = nbo_prediction_for(f"{mol_id}_{other}", mols_by_name)
        if other_prediction == "pending":
            pending.append(mol_id)
        already_run_prediction = existing.get(f"{mol_id}_{already_run_isomer}", "")

        preds_by_isomer = {already_run_isomer: already_run_prediction, other: other_prediction}
        pct_a, pct_b = product_pcts.get(mol_id, ("", ""))

        out_row = {
            "mol_id": mol_id,
            "lowest_energy_isomer": lowest,
            "other_isomer": other,
            "prediction_lowest_energy": preds_by_isomer.get(lowest, ""),
            "prediction_other_isomer": preds_by_isomer.get(other, ""),
            "pct_product_A": pct_a,
            "pct_product_B": pct_b,
        }
        out_rows.append(out_row)
        print(f"  {out_row['mol_id']}: lowest={out_row['lowest_energy_isomer']}({out_row['prediction_lowest_energy']}) "
              f"other={out_row['other_isomer']}({out_row['prediction_other_isomer']}) "
              f"pct_A={out_row['pct_product_A']} pct_B={out_row['pct_product_B']}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"\nWrote {OUT_CSV}")
    if pending:
        print(f"PENDING (Stage 3 scan not yet complete): {', '.join(pending)}")
    else:
        print("All 9 molecules' Stage 3 scans complete -- no pending entries.")


if __name__ == "__main__":
    main()
