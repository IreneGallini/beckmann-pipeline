"""`beckmann-nbo predict` -- SMILES (or a benchmark.csv-shaped CSV) to a
submitted Stage 1+2 DFT job, and later (via --continue) Stage 3.

This is a thin orchestration layer: every step below calls straight into
an existing beckmann_core/beckmann_nbo function. Non-blocking by design --
submitting returns immediately with a name to poll via `beckmann-nbo status`.

Stage 3 (the N-O scan `status`/`report` need for a wCNmax prediction) can
only be generated after Stage 1's converged geometry has been downloaded
from the cluster, so it isn't part of the initial `predict` call -- run
`predict --continue <query_id>` once `status` shows Stage 1 complete.
"""
import csv
import sys
from pathlib import Path

from rdkit import Chem

from beckmann_core.oximes import KETONE_PAT, enumerate_ez, ketone_to_protonated_oximes
from beckmann_core.conformers import generate_conformers
from beckmann_core.optimize import select_and_optimize
from beckmann_nbo.config import DATA_OUTPUT
from beckmann_nbo.hpc import cmd_submit_opt, cmd_submit_scan, cmd_upload, load_config, require_config
from beckmann_nbo.inputs import prepare_opt, prepare_scan_rigid

QUERY_PREFIX = "mol"  # must be "mol" -- beckmann_nbo.hpc.mol_dirs()/_remote_dir_names()
                       # hardcode a "mol_*" glob (not configurable, and out of scope to
                       # change here), so query directories have to use that same prefix
                       # or upload/submit-opt/submit-scan would silently find nothing.


def _sanitize_id(query_id: str) -> str:
    """prepare_opt()'s test_ids filter does name.split('_')[1] on the
    3-token '{prefix}_{id}_{E|Z}' convention (mirroring 'mol_002_E') -- a
    query_id containing '_' would shift that split, so collapse underscores
    before use. Prefixed with 'q' so a query never collides with a real
    3-digit benchmark id (001-034), then zfill(3)'d to match exactly what
    beckmann_nbo.hpc.mol_dirs() does internally to whatever id it's passed
    (mol.zfill(3)) -- doing it ourselves up front keeps the directory name
    we create and the name hpc.py later globs for byte-identical."""
    return ("q" + query_id.replace("_", "-")).zfill(3)


def smiles_to_oximes(smiles: str) -> list:
    """SMILES -> list of protonated-oxime RDKit Mols. Raises ValueError on
    an unparseable SMILES or one with no ketone/oxime-convertible group --
    same validation beckmann-pyscf's parse_and_check_ketone() does, kept
    inline here rather than importing across the product-package boundary
    (beckmann-nbo's own dependency footprint stays beckmann-core + rdkit)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"could not parse SMILES: {smiles!r}")
    if not mol.HasSubstructMatch(KETONE_PAT):
        raise ValueError(f"no ketone group found in: {smiles!r}")
    oximes = ketone_to_protonated_oximes(mol)
    if not oximes:
        raise ValueError(f"ketone-to-oxime conversion produced no products for: {smiles!r}")
    return oximes


def write_smi(smiles: str, query_id: str, workdir: Path) -> Path:
    """Writes '{QUERY_PREFIX}_{query_id}_{E|Z}' isomer names -- the same
    3-token '{prefix}_{id}_{isomer}' shape the benchmark pipeline uses
    (e.g. 'mol_002_E'), required by prepare_opt()'s test_ids filter and by
    select_and_optimize()'s per-substrate grouping (name.rsplit('_', 1)[0])."""
    oximes = smiles_to_oximes(smiles)
    smi_path = workdir / f"{QUERY_PREFIX}_{query_id}.smi"
    workdir.mkdir(parents=True, exist_ok=True)
    with open(smi_path, "w") as f:
        for ox in oximes:
            for iso, ez_label in enumerate_ez(ox):
                f.write(f"{Chem.MolToSmiles(iso)} {QUERY_PREFIX}_{query_id}_{ez_label}\n")
    return smi_path


def submit_stage_1_2(query_id: str, smi_path: Path, workdir: Path, dry_run: bool) -> None:
    print(f"\n[1/4] Generating conformers for {query_id}...")
    sdf_path = generate_conformers(smi_path, workdir / "conformers")

    print(f"[2/4] AIMNet2 geometry optimization for {query_id}...")
    _, sub_sdf = select_and_optimize(sdf_path, workdir / "aimnet_optimized")

    print(f"[3/4] Writing Stage 1+2 Gaussian input files for {query_id}...")
    dft_opt_dir = workdir / "dft_opt"
    prepare_opt(sub_sdf, dft_opt_dir, test_ids={query_id})

    print(f"[4/4] Uploading and submitting Stage 1 (geometry opt) on the cluster...")
    config = load_config()
    require_config(config)
    cmd_upload(config, dry_run, query_id, dft_opt_dir)
    cmd_submit_opt(config, dry_run, query_id, dft_opt_dir)

    print(
        f"\nSubmitted. Poll with:\n"
        f"  beckmann-nbo status --mol {query_id} --dir {dft_opt_dir}\n"
        f"Once Stage 1 shows Normal termination, continue to Stage 3 with:\n"
        f"  beckmann-nbo predict --continue {query_id} --dir {dft_opt_dir}"
    )


def continue_stage_3(query_id: str, dft_opt_dir: Path, dry_run: bool) -> None:
    matches = sorted(dft_opt_dir.glob(f"{QUERY_PREFIX}_{query_id}_*"))
    if not matches:
        print(f"ERROR: no directory matching {QUERY_PREFIX}_{query_id}_* under {dft_opt_dir}", file=sys.stderr)
        sys.exit(1)
    mol_dir = matches[0]
    mol_name = mol_dir.name

    opt_log = mol_dir / f"{mol_name}_opt.log"
    if not opt_log.exists() or "Normal termination" not in opt_log.read_text()[-2000:]:
        print(
            f"ERROR: {opt_log} not found or Stage 1 hasn't reached Normal termination yet.\n"
            f"Check with: beckmann-nbo status --mol {query_id} --dir {dft_opt_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Generating Stage 3 (N-O scan) for {mol_name}...")
    prepare_scan_rigid(mol_dir, mol_name)

    config = load_config()
    require_config(config)
    cmd_upload(config, dry_run, query_id, dft_opt_dir)
    cmd_submit_scan(config, dry_run, query_id, dft_opt_dir)
    print(
        f"\nSubmitted. Poll with:\n"
        f"  beckmann-nbo status --mol {query_id} --dir {dft_opt_dir}"
    )


def cmd_predict(args) -> None:
    dry_run = args.dry_run
    workdir = Path(args.workdir) if args.workdir else DATA_OUTPUT / "query_predictions"

    if args.continue_name:
        query_id = _sanitize_id(args.continue_name)
        dft_opt_dir = Path(args.dir) if args.dir else workdir / query_id / "dft_opt"
        continue_stage_3(query_id, dft_opt_dir, dry_run)
        return

    if args.csv:
        csv_path = Path(args.csv)
        with open(csv_path, newline="") as f:
            rows = list(csv.DictReader(f))
        if not rows or "id" not in rows[0] or "SMILES" not in rows[0]:
            print(
                "ERROR: --csv file must have 'id' and 'SMILES' columns "
                "(same shape as data/input/benchmark.csv)",
                file=sys.stderr,
            )
            sys.exit(1)
        for row in rows:
            query_id = _sanitize_id(str(row["id"]))
            mol_workdir = workdir / query_id
            try:
                smi_path = write_smi(row["SMILES"], query_id, mol_workdir)
            except ValueError as e:
                print(f"SKIP {query_id}: {e}", file=sys.stderr)
                continue
            submit_stage_1_2(query_id, smi_path, mol_workdir, dry_run)
        return

    if not args.smiles:
        print("ERROR: pass --smiles (with --name) or --csv", file=sys.stderr)
        sys.exit(1)

    query_id = _sanitize_id(args.name if args.name else "query")
    mol_workdir = workdir / query_id
    try:
        smi_path = write_smi(args.smiles, query_id, mol_workdir)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    submit_stage_1_2(query_id, smi_path, mol_workdir, dry_run)
