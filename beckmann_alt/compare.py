"""
Compare both open-source prototypes (LIVVO, crude AO-projection) against the trusted
NBO7 numbers for the two reference cases. Trusted numbers are PULLED from existing
artifacts, never recomputed by a different method here:
  - mol_002: read directly from the main branch's
    data/output/analysis/cmo_channel_extraction.csv (nbo/R0 stage).
  - 5_s0_Me: no pre-computed CSV row exists for this external, non-benchmark reference
    file, so its trusted numbers come from re-running the SAME already-validated
    beckmann.dft.parse_cmo extraction functions that scripts/analysis/
    validate_reference_descriptors.py itself uses (find_cmo_sections/parse_cmo_table/
    virtual_window/max_weight_for_target) directly on 5_s0_Me.log -- this is calling
    the existing trusted parser, not deriving the numbers a different way.

See Notes_open_source_alt.md for the full writeup of what these numbers mean.
"""
import csv
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft.parse_cmo import (
    find_cmo_sections, max_weight_for_target, parse_cmo_table, virtual_window,
)

from beckmann_alt import ao_projection, pyscf_livvo
from beckmann_alt.geometry import REFERENCE_CASES, REFERENCE_LOG

CHANNEL_ATOMS = {  # channel -> (atom_a, atom_b) key names in a REFERENCE_CASES entry
    "cn": ("ci", "ni"), "w17": ("ci", "c_aryl"), "w78": ("ci", "c_alkyl"),
}


def trusted_mol_002() -> dict:
    """nbo/R0-stage cn/17/78 rows for mol_002_E from the already-generated CSV."""
    path = DATA_OUTPUT / "analysis" / "cmo_channel_extraction.csv"
    rows = {r["channel"]: r for r in csv.DictReader(open(path)) if r["mol"] == "mol_002_E" and r["stage"] == "nbo"}
    return {
        "cn":  {"wmax": float(rows["cn"]["weight"]),  "mo_index": rows["cn"]["MO_index"]},
        "w17": {"wmax": float(rows["17"]["weight"]),  "mo_index": rows["17"]["MO_index"]},
        "w78": {"wmax": float(rows["78"]["weight"]),  "mo_index": rows["78"]["MO_index"]},
    }


def trusted_5_s0_me() -> dict:
    """Re-run the existing, already-validated parse_cmo extraction on 5_s0_Me.log
    (Tier 1 of validate_reference_descriptors.py) -- same functions, not a new method."""
    case = REFERENCE_CASES["5_s0_Me"]
    lines = REFERENCE_LOG.read_text().splitlines()
    cmo_starts = find_cmo_sections(lines)
    table = parse_cmo_table(lines, cmo_starts[-1])  # last section, see Stable=Opt note
    window = virtual_window(table)

    out = {}
    for key, (a_name, b_name) in CHANNEL_ATOMS.items():
        a, b = case[a_name], case[b_name]
        wmax, mo, eps, coeff, _, _ = max_weight_for_target(window, a, b)
        out[key] = {"wmax": wmax, "mo_index": mo}
    return out


def rank(d: dict) -> list[str]:
    """Channel names ordered by descending wmax (None sorts last)."""
    return sorted(d, key=lambda k: (d[k]["wmax"] is None, -(d[k]["wmax"] or 0)))


def report_case(name: str, trusted: dict) -> None:
    print(f"\n{'='*70}\n{name}\n{'='*70}")
    livvo = pyscf_livvo.run_case(name)
    livvo = {"cn": livvo["cn"], "w17": livvo["w17"], "w78": livvo["w78"]}
    crude = ao_projection.run_case(name)

    print(f"{'channel':<8} {'NBO7 (trusted)':>16} {'LIVVO':>16} {'crude AO-proj':>16}")
    for ch in ["cn", "w17", "w78"]:
        t = trusted[ch]["wmax"]
        l = livvo[ch]["wmax"] if livvo[ch].get("wmax") is not None else livvo[ch].get("livvo", {}).get("combined_pop")
        l = livvo[ch]["wmax"]
        c = crude[ch]["wmax"]
        fmt = lambda v: f"{v:.4f}" if v is not None else "None"
        print(f"{ch:<8} {fmt(t):>16} {fmt(l):>16} {fmt(c):>16}")

    print(f"\nranking (highest wmax first):")
    print(f"  NBO7 (trusted): {rank(trusted)}")
    print(f"  LIVVO:          {rank(livvo)}")
    print(f"  crude AO-proj:  {rank(crude)}")

    livvo_agrees = rank(trusted)[0] == rank(livvo)[0]
    crude_agrees = rank(trusted)[0] == rank(crude)[0]
    print(f"\ntop channel agrees with NBO7?  LIVVO: {livvo_agrees}   crude AO-proj: {crude_agrees}")

    for label, r in [("LIVVO", livvo), ("crude AO-proj", crude)]:
        for ch in ["cn", "w17", "w78"]:
            if r[ch]["wmax"] is None:
                print(f"  NOTE: {label} could not identify a channel orbital for '{ch}' at all")


def main() -> None:
    print("Pulling trusted NBO7 numbers (not recomputed)...")
    t_mol002 = trusted_mol_002()
    t_ref    = trusted_5_s0_me()

    report_case("mol_002", t_mol002)
    report_case("5_s0_Me", t_ref)


if __name__ == "__main__":
    main()
