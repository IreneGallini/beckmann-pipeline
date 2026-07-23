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
    per-molecule result: exp/NBO/PySCF predicted R/F labels side by side,
    correctness/cross-method-match columns, then n_points/R_star/R_depth
    -- PySCF's own scan diagnostics -- as the last three columns)
"""
import csv
import json
import sys
from pathlib import Path

from beckmann.config import DATA_INPUT, DATA_OUTPUT
from beckmann.dft.descriptors import find_wcnmax_minimum
from beckmann.dft.inputs import ALL_IDS, TEST_IDS
from beckmann.dft.wcnmax_rule import predict_from_wcnmax

from beckmann_alt.pair_nbo import run_test_set_scan_series

EXTRACTION_FIELDS = [
    "mol", "stage", "channel", "R_NO", "MO_index", "epsilon_i_star", "coefficient",
    "weight", "delta_lumo", "in_window",
]
RESULTS_FIELDS = [
    "mol", "exp", "NBO", "PySCF",
    "NBO_correct", "PySCF_correct", "NBO_PySCF_match",
    "NBO_minimum_found", "PySCF_minimum_found",
    "n_points", "R_star", "R_depth",
]

# mol_034_E's STEP_SCAN_SOURCES-merged series has 12 points (2x every other
# molecule's 6), which kept getting caught mid-run by interruptions on this
# machine. Every-other-point subset, shifted to land on scan_6 (R=1.669), the
# exact R of NBO7's own trusted interior minimum for this molecule (R_star=
# 1.6690 in wcnmax_rule_results.csv) -- keeps the minimum visible with points
# on both sides at roughly 6-point-series resolution.
STAGE_OVERRIDE = {"034": ["scan_2", "scan_4", "scan_6", "scan_8", "scan_10", "scan_12"]}


def trusted_row(mol_name: str) -> dict:
    path = DATA_OUTPUT / "analysis" / "wcnmax_rule_results.csv"
    for row in csv.DictReader(open(path)):
        if row["mol"] == mol_name:
            return row
    raise ValueError(f"{mol_name}: no row in {path}")


def merge_write_csv(path: Path, fieldnames: list[str], new_rows: list[dict]) -> None:
    """Read-merge-write instead of blind overwrite: keeps every existing row whose
    'mol' isn't in new_rows (in its existing order), then appends new_rows. Called
    once per molecule (not once at the end of the whole run) so a later molecule's
    crash can't discard an earlier molecule's already-computed result."""
    existing = list(csv.DictReader(open(path))) if path.exists() else []
    new_mols = {r["mol"] for r in new_rows}
    merged = [r for r in existing if r["mol"] not in new_mols] + new_rows
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(merged)


def main() -> None:
    mol_ids = sys.argv[1:] or sorted(TEST_IDS)
    bad = [m for m in mol_ids if m not in ALL_IDS]
    if bad:
        raise ValueError(f"not in ALL_IDS: {bad}")

    outcomes = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())
    analysis_dir = DATA_OUTPUT / "analysis"
    extraction_path = analysis_dir / "wcnmax_channel_extraction_opensource.csv"
    results_path = analysis_dir / "wcnmax_rule_results_opensource.csv"

    header = (
        f"{'mol':<12} {'n_pts':>5} {'min_found':>10} {'pyscf':>6} {'exp':>4} {'agree':>6}"
        f"   || {'nbo min':>8} {'nbo':>6} {'nbo agree':>10}"
    )
    print(header)
    print("-" * len(header))

    open_source_agree = 0
    n_done = 0
    for mol_id in mol_ids:
        print(f"-- running {mol_id}...", flush=True)
        try:
            rows = run_test_set_scan_series(mol_id, stages=STAGE_OVERRIDE.get(mol_id))
            mol_name = rows[0]["mol"]

            minimum = find_wcnmax_minimum(mol_name, rows)
            predicted = predict_from_wcnmax(minimum)
            mol_num = mol_name.split("_")[1]
            exp = outcomes.get(f"mol_{mol_num}", {}).get("exp_outcome", "")
            agreement = "yes" if predicted == exp else "no"

            trusted = trusted_row(mol_name)
            result_row = {
                "mol": mol_name, "exp": exp,
                "NBO": trusted["predicted"], "PySCF": predicted,
                "NBO_correct": trusted["agreement"], "PySCF_correct": agreement,
                "NBO_PySCF_match": "yes" if predicted == trusted["predicted"] else "no",
                "NBO_minimum_found": trusted["minimum_found"], "PySCF_minimum_found": minimum is not None,
                "n_points": len(rows),
                "R_star": f"{minimum['R_star']:.4f}" if minimum else "",
                "R_depth": f"{minimum['depth']:.4f}" if minimum else "",
            }
        except Exception as exc:
            print(f"-- {mol_id}: FAILED ({exc}), skipping")
            continue

        merge_write_csv(extraction_path, EXTRACTION_FIELDS, rows)
        merge_write_csv(results_path, RESULTS_FIELDS, [result_row])
        n_done += 1
        if agreement == "yes":
            open_source_agree += 1

        print(
            f"{mol_name:<12} {len(rows):>5} {str(minimum is not None):>10} {predicted:>6} "
            f"{exp:>4} {agreement:>6}   || "
            f"{trusted['minimum_found']:>8} {trusted['predicted']:>6} {trusted['agreement']:>10}"
        )

    print(f"\nPySCF wCNmax rule agreement with experiment (this run): {open_source_agree}/{n_done}")
    print(f"-> {extraction_path}")
    print(f"-> {results_path}")


if __name__ == "__main__":
    main()
