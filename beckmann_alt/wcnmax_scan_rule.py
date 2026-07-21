"""
Reproduce beckmann.dft.wcnmax_rule's "genuine interior wCNmax minimum ->
predict rearrangement" benchmark using the open-source per-atom-pair wCNmax
(beckmann_alt.pair_nbo) instead of NBO7, for the 6 main-pipeline test-set
molecules (mol_002/006/014/020/021/029).

For each molecule: compute wCNmax at every point of its R(N-O) scan series
(beckmann_alt.pair_nbo.run_test_set_scan_series), feed the resulting rows
into beckmann.dft.descriptors.find_wcnmax_minimum() -- the SAME minimum-
detection code the trusted NBO7 benchmark uses (beckmann/dft/wcnmax_rule.py)
-- and predict R/F via beckmann.dft.wcnmax_rule.predict_from_wcnmax(). This
checks whether the open-source method's own R(N-O) trend reproduces the
rule's classification, not just whether its equilibrium wCNmax value is
numerically close to NBO7's (see beckmann_alt/compare_wcnmax.py for that
narrower single-point check).

Trusted comparison numbers (NBO7-based minimum_found/predicted/agreement)
are read directly from the already-generated
data/output/analysis/wcnmax_rule_results.csv, never recomputed differently.

Each molecule is a real DFT scan (6-13 PySCF single-points), so this is slow
(order of an hour total on a laptop, most of it mol_020_E's 13-point series)
-- see beckmann_alt/_compute_scan_cache.py for running molecules in parallel
as separate processes instead of serially through this script.

Output:
  data/output/analysis/wcnmax_channel_extraction_opensource.csv (per-point
    'cn'-channel rows, same shape as beckmann.dft.parse_cmo's trusted
    cmo_channel_extraction.csv)
  data/output/analysis/wcnmax_rule_results_opensource.csv (condensed
    per-molecule result, same shape as the trusted wcnmax_rule_results.csv
    plus trusted_*/matches_trusted_prediction comparison columns)
"""
import csv
import json

from beckmann.config import DATA_INPUT, DATA_OUTPUT
from beckmann.dft.descriptors import find_wcnmax_minimum
from beckmann.dft.inputs import TEST_IDS
from beckmann.dft.wcnmax_rule import predict_from_wcnmax

from beckmann_alt.pair_nbo import run_test_set_scan_series

EXTRACTION_FIELDS = [
    "mol", "stage", "channel", "R_NO", "MO_index", "epsilon_i_star", "coefficient",
    "weight", "delta_lumo", "in_window",
]
RESULTS_FIELDS = [
    "mol", "n_points", "minimum_found", "R_star", "depth", "predicted", "exp_outcome",
    "agreement", "trusted_minimum_found", "trusted_predicted", "trusted_agreement",
    "matches_trusted_prediction",
]


def trusted_row(mol_name: str) -> dict:
    path = DATA_OUTPUT / "analysis" / "wcnmax_rule_results.csv"
    for row in csv.DictReader(open(path)):
        if row["mol"] == mol_name:
            return row
    raise ValueError(f"{mol_name}: no row in {path}")


def main() -> None:
    outcomes = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())
    analysis_dir = DATA_OUTPUT / "analysis"

    header = (
        f"{'mol':<12} {'n_pts':>5} {'min_found':>10} {'predicted':>9} {'exp':>4} {'agree':>6}"
        f"   || {'trusted min':>11} {'trusted pred':>12} {'trusted agree':>13}"
    )
    print(header)
    print("-" * len(header))

    all_extraction_rows = []
    result_rows = []
    open_source_agree = 0
    for mol_id in sorted(TEST_IDS):
        print(f"-- running {mol_id}...", flush=True)
        rows = run_test_set_scan_series(mol_id)
        mol_name = rows[0]["mol"]
        all_extraction_rows.extend(rows)

        minimum = find_wcnmax_minimum(mol_name, rows)
        predicted = predict_from_wcnmax(minimum)
        mol_num = mol_name.split("_")[1]
        exp = outcomes.get(f"mol_{mol_num}", {}).get("exp_outcome", "")
        agreement = "yes" if predicted == exp else "no"
        if agreement == "yes":
            open_source_agree += 1

        trusted = trusted_row(mol_name)
        result_rows.append({
            "mol": mol_name, "n_points": len(rows),
            "minimum_found": minimum is not None,
            "R_star": f"{minimum['R_star']:.4f}" if minimum else "",
            "depth": f"{minimum['depth']:.4f}" if minimum else "",
            "predicted": predicted, "exp_outcome": exp, "agreement": agreement,
            "trusted_minimum_found": trusted["minimum_found"],
            "trusted_predicted": trusted["predicted"],
            "trusted_agreement": trusted["agreement"],
            "matches_trusted_prediction": "yes" if predicted == trusted["predicted"] else "no",
        })

        print(
            f"{mol_name:<12} {len(rows):>5} {str(minimum is not None):>10} {predicted:>9} "
            f"{exp:>4} {agreement:>6}   || "
            f"{trusted['minimum_found']:>11} {trusted['predicted']:>12} {trusted['agreement']:>13}"
        )

    print(f"\nopen-source wCNmax rule agreement with experiment: {open_source_agree}/{len(result_rows)}")

    extraction_path = analysis_dir / "wcnmax_channel_extraction_opensource.csv"
    with open(extraction_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=EXTRACTION_FIELDS)
        writer.writeheader()
        writer.writerows(all_extraction_rows)
    print(f"-> {extraction_path}")

    results_path = analysis_dir / "wcnmax_rule_results_opensource.csv"
    with open(results_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=RESULTS_FIELDS)
        writer.writeheader()
        writer.writerows(result_rows)
    print(f"-> {results_path}")


if __name__ == "__main__":
    main()
