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
"""
import csv
import json

from beckmann.config import DATA_INPUT, DATA_OUTPUT
from beckmann.dft.descriptors import find_wcnmax_minimum
from beckmann.dft.inputs import TEST_IDS
from beckmann.dft.wcnmax_rule import predict_from_wcnmax

from beckmann_alt.pair_nbo import run_test_set_scan_series


def trusted_row(mol_name: str) -> dict:
    path = DATA_OUTPUT / "analysis" / "wcnmax_rule_results.csv"
    for row in csv.DictReader(open(path)):
        if row["mol"] == mol_name:
            return row
    raise ValueError(f"{mol_name}: no row in {path}")


def main() -> None:
    outcomes = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    header = (
        f"{'mol':<12} {'n_pts':>5} {'min_found':>10} {'predicted':>9} {'exp':>4} {'agree':>6}"
        f"   || {'trusted min':>11} {'trusted pred':>12} {'trusted agree':>13}"
    )
    print(header)
    print("-" * len(header))

    open_source_agree = 0
    total = 0
    for mol_id in sorted(TEST_IDS):
        print(f"-- running {mol_id}...", flush=True)
        rows = run_test_set_scan_series(mol_id)
        mol_name = rows[0]["mol"]

        minimum = find_wcnmax_minimum(mol_name, rows)
        predicted = predict_from_wcnmax(minimum)
        mol_num = mol_name.split("_")[1]
        exp = outcomes.get(f"mol_{mol_num}", {}).get("exp_outcome", "")
        agreement = "yes" if predicted == exp else "no"
        if agreement == "yes":
            open_source_agree += 1
        total += 1

        trusted = trusted_row(mol_name)

        print(
            f"{mol_name:<12} {len(rows):>5} {str(minimum is not None):>10} {predicted:>9} "
            f"{exp:>4} {agreement:>6}   || "
            f"{trusted['minimum_found']:>11} {trusted['predicted']:>12} {trusted['agreement']:>13}"
        )

    print(f"\nopen-source wCNmax rule agreement with experiment: {open_source_agree}/{total}")


if __name__ == "__main__":
    main()
