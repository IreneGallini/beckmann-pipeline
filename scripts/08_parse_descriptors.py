"""Parse the downloaded Stage 1-3 .log files into per-molecule descriptor
CSVs: E2PERT (nbo_e2pert.csv) and the CMO/wCNmax channel data
(cmo_descriptors.csv, cmo_channel_extraction.csv the latter is what
09_predict_rf.py reads the R/F prediction off of). Needs Stage 3 downloaded
(05_download_results.py after 07_upload_submit_stage3.py).

Edit MOL_NAME below, then:
    python 08_parse_descriptors.py
"""
import csv
import sys

from _common import QUERY_PREFIX, local_substituent_map, sanitize_id, workdir_for

from beckmann_nbo.parse_cmo import EXTRACTION_FIELDS, FIELDS as CMO_FIELDS, collect_molecule as collect_cmo
from beckmann_nbo.parse_nbo import FIELDS as NBO_FIELDS, collect_molecule as collect_e2pert

MOL_NAME = "test1"


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    dft_opt_dir = workdir_for(mol_id) / "dft_opt"
    analysis_dir = workdir_for(mol_id) / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    matches = sorted(dft_opt_dir.glob(f"{QUERY_PREFIX}_{mol_id}_*"))
    if not matches:
        print(f"ERROR: no directory matching {QUERY_PREFIX}_{mol_id}_* under {dft_opt_dir}", file=sys.stderr)
        sys.exit(1)
    mol_dir = matches[0]
    mol_name = mol_dir.name

    subst = local_substituent_map(mol_name, mol_dir, workdir_for(mol_id))
    print(f"{mol_name} (aryl=C{subst['c_aryl']}, alkyl=C{subst['c_alkyl']})")

    e2pert_rows = collect_e2pert(mol_name, mol_dir)
    e2pert_path = analysis_dir / "nbo_e2pert.csv"
    with open(e2pert_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=NBO_FIELDS)
        writer.writeheader()
        writer.writerows(e2pert_rows)
    print(f"{len(e2pert_rows)} E2PERT rows -> {e2pert_path}")

    cmo_rows, channel_rows = collect_cmo(mol_name, mol_dir, subst["c_aryl"], subst["c_alkyl"])
    cmo_path = analysis_dir / "cmo_descriptors.csv"
    with open(cmo_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CMO_FIELDS)
        writer.writeheader()
        writer.writerows(cmo_rows)
    print(f"{len(cmo_rows)} CMO stage rows -> {cmo_path}")

    extraction_path = analysis_dir / "cmo_channel_extraction.csv"
    with open(extraction_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=EXTRACTION_FIELDS)
        writer.writeheader()
        writer.writerows(channel_rows)
    print(f"{len(channel_rows)} channel-extraction rows -> {extraction_path}")

    print(f"\nNext: edit MOL_NAME in 09_predict_rf.py to '{MOL_NAME}' and run it.")


if __name__ == "__main__":
    main()
