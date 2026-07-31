"""
wCNmax minima rule benchmark: does an interior wCNmax(R) minimum predict
Beckmann rearrangement, across the full 34-molecule benchmark set?
find_wcnmax_minimum()/predict_from_wcnmax() themselves moved to
beckmann_core.wcnmax_rule (geometry-source-agnostic, used live by both
products) -- this file has only the benchmark-CSV-driving harness around
them, moved verbatim from beckmann/dft/wcnmax_rule.py's
run_wcnmax_benchmark()/main().

Output: data/output/analysis/wcnmax_rule_results.csv
"""
import csv
import json

from beckmann_core.wcnmax_rule import find_wcnmax_minimum, predict_from_wcnmax
from beckmann_nbo.config import DATA_INPUT, DATA_OUTPUT
from beckmann_nbo.descriptors import load_series

FIELDS = ["mol", "n_points", "minimum_found", "R_star", "depth", "predicted", "exp_outcome", "agreement"]


def run_wcnmax_benchmark(channel_rows: list[dict], extraction_rows: list[dict],
                          outcomes: dict) -> tuple[list[dict], float]:
    """Run the wCNmax minima rule over every molecule with a usable series
    in channel_rows. Returns (result rows, accuracy fraction)."""
    mols = sorted({row["mol"] for row in channel_rows})
    rows = []
    for mol in mols:
        r_values, _ = load_series(mol, channel_rows)
        if not r_values:
            continue
        minimum = find_wcnmax_minimum(mol, extraction_rows)
        pred = predict_from_wcnmax(minimum)

        mol_id = mol.split("_")[1]
        exp = outcomes.get(f"mol_{mol_id}", {}).get("exp_outcome", "")
        agreement = "yes" if pred == exp else "no"

        rows.append({
            "mol": mol, "n_points": len(r_values),
            "minimum_found": minimum is not None,
            "R_star": f"{minimum['R_star']:.4f}" if minimum else "",
            "depth": f"{minimum['depth']:.4f}" if minimum else "",
            "predicted": pred, "exp_outcome": exp, "agreement": agreement,
        })
        print(f"  {mol}  minimum={'yes' if minimum else 'no ':<3}  pred={pred}  exp={exp}  -> {agreement}")

    total = len(rows)
    agree = sum(1 for r in rows if r["agreement"] == "yes")
    accuracy = agree / total if total else 0.0
    return rows, accuracy


def main() -> None:
    analysis_dir = DATA_OUTPUT / "analysis"
    with open(analysis_dir / "channel_descriptors.csv", newline="") as f:
        channel_rows = list(csv.DictReader(f))
    with open(analysis_dir / "cmo_channel_extraction.csv", newline="") as f:
        extraction_rows = list(csv.DictReader(f))
    outcomes = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    rows, accuracy = run_wcnmax_benchmark(channel_rows, extraction_rows, outcomes)

    out_path = analysis_dir / "wcnmax_rule_results.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    agree = sum(1 for r in rows if r["agreement"] == "yes")
    print(f"\nResults -> {out_path}")
    print(f"wCNmax rule accuracy: {agree}/{len(rows)} ({100 * accuracy:.0f}%)")


if __name__ == "__main__":
    main()
