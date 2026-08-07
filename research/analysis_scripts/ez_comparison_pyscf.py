"""
Step 2 of the priority E/Z re-run: run the PySCF-native wCNmax scan on the
*other* isomer's already-optimized AIMNet2 geometry (from
best_aimnet_optimized.sdf) for each of the 9 priority molecules, and combine
with the already-computed prediction for the already-run isomer (read from
wcnmax_rule_results_pyscf.csv) into one comparison table.

Deliberately bypasses beckmann_pyscf.pipeline.predict() -- that re-runs
SMILES -> Auto3D -> AIMNet2 from scratch and auto-picks the lower-energy
isomer, which is exactly the redundant work/auto-selection this task is
checking for a conformer-selection artifact around. Instead it calls
run_scan_series()/predict_outcome() directly on the isomer we choose,
mirroring beckmann_pyscf.pipeline.__init__.predict()'s own call chain
(packages/beckmann-pyscf/backend/beckmann_pyscf/pipeline/__init__.py:70-86).

run_scan_series() uses a fixed absolute R(N-O) window (1.50-1.80 A), not an
R0-relative one -- see its module docstring. Left unmodified here so these 9
results stay comparable to the existing 34-molecule
wcnmax_rule_results_pyscf.csv; the "R0 to 1.85 A" window from the task spec
is NBO-side-native and applied only in ez_comparison_nbo.py.
"""
import csv
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from rdkit import Chem

from beckmann_core.classical import get_oxime_atoms
from beckmann_nbo.config import DATA_INPUT, DATA_OUTPUT

from ez_energy_check import energy_table

BEST_AIMNET_SDF = DATA_OUTPUT / "aimnet_optimized" / "best_aimnet_optimized.sdf"
EXISTING_RESULTS_CSV = DATA_OUTPUT / "analysis" / "wcnmax_rule_results_pyscf.csv"
BENCHMARK_CSV = DATA_INPUT / "benchmark.csv"
OUT_CSV = DATA_OUTPUT / "analysis" / "ez_comparison_pyscf.csv"

FIELDS = ["mol_id", "lowest_energy_isomer", "other_isomer",
          "prediction_lowest_energy", "prediction_other_isomer",
          "pct_product_A", "pct_product_B"]


def load_existing_predictions(csv_path=EXISTING_RESULTS_CSV) -> dict[str, str]:
    """{'mol_001_E': 'R', ...} from the PySCF column of the existing 34-molecule results."""
    out = {}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            out[row["mol"]] = row["PySCF"]
    return out


def load_product_pcts(csv_path=BENCHMARK_CSV) -> dict[str, tuple[float, float]]:
    """{'mol_001': (33.0, 67.0), ...} -- (% product A, % product B) per molecule id."""
    out = {}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            mol_id = f"mol_{int(row['id']):03d}"
            out[mol_id] = (float(row["% product A"]), float(row["% product B"]))
    return out


def run_other_isomer_pyscf(mol_id: str, other_isomer: str, mols_by_name: dict) -> str:
    from beckmann_pyscf.pipeline.wcnmax_pyscf import run_scan_series
    from beckmann_pyscf.pipeline.predict import predict_outcome

    name = f"{mol_id}_{other_isomer}"
    mol = mols_by_name[name]
    atom_ids = get_oxime_atoms(mol)
    if atom_ids is None:
        return "inspect"
    cox, nox, oox, c_aryl, c_allyl = atom_ids
    print(f"  [{name}] running PySCF wCNmax scan...", flush=True)
    rows = run_scan_series(mol, ci=cox + 1, ni=nox + 1, oi=oox + 1,
                            c_aryl=c_aryl + 1, c_alkyl=c_allyl + 1, name=name)
    prediction, _minimum = predict_outcome(name, rows)
    print(f"  [{name}] prediction = {prediction}", flush=True)
    return prediction


def load_done_mol_ids(csv_path=OUT_CSV) -> set[str]:
    """Resume support: mol_ids already written to OUT_CSV (e.g. after an
    interrupted prior run) are skipped rather than recomputed."""
    if not csv_path.exists():
        return set()
    with open(csv_path, newline="") as f:
        return {row["mol_id"] for row in csv.DictReader(f)}


def main(mol_ids=None) -> None:
    from ez_energy_check import ALREADY_RUN
    rows = energy_table(mol_ids=mol_ids if mol_ids is not None else tuple(ALREADY_RUN))
    existing = load_existing_predictions()
    product_pcts = load_product_pcts()
    done = load_done_mol_ids()
    if done:
        print(f"Resuming: already have results for {sorted(done)}")

    suppl = Chem.SDMolSupplier(str(BEST_AIMNET_SDF), removeHs=False)
    mols_by_name = {m.GetProp("_Name"): m for m in suppl if m is not None}

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    write_header = not OUT_CSV.exists()
    with open(OUT_CSV, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if write_header:
            writer.writeheader()

        for r in rows:
            mol_id = r["mol_id"]
            if mol_id in done:
                continue
            lowest, other = r["lowest_energy_isomer"], r["other_isomer"]
            already_run_isomer = r["already_run_isomer"]

            other_prediction = run_other_isomer_pyscf(mol_id, other, mols_by_name)
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
            writer.writerow(out_row)
            f.flush()
            print(f"  {out_row['mol_id']}: lowest={out_row['lowest_energy_isomer']}({out_row['prediction_lowest_energy']}) "
                  f"other={out_row['other_isomer']}({out_row['prediction_other_isomer']}) "
                  f"pct_A={out_row['pct_product_A']} pct_B={out_row['pct_product_B']}")

    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
