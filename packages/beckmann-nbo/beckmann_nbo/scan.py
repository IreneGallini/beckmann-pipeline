"""
LEGACY (old internal-walk scan architecture only — see RIGID_SCAN_MIGRATION.md
for the current architecture and why this was superseded). Extract converged
intermediate geometries from a Gaussian scan log and write single-point NBO
input files for the missing scan points.

Under the old architecture, the scan (Stage 3) only ran NBO at R0 and R0+0.4 Å.
This module extracted the converged geometries at R0+0.1, R0+0.2, R0+0.3 and
created three single-point .gjf files for upload to Citadel. The current
rigid-scan architecture (`beckmann/dft/inputs.py::prepare_scan_rigid()`) gets
full NBO natively at every point, so no test molecule needs this anymore —
kept for reference, not part of the active pipeline.

Output: data/output/dft_opt/{mol}/{mol}_sp{N}.gjf  (N = 2, 3, 4)
"""
import argparse
import re
from pathlib import Path

from beckmann_core.constants import CHARGE, MULTIPLICITY
from beckmann_core.geometry import (
    ATOMIC_SYMBOLS, displace_leaving_group, find_leaving_group, no_distance,
)
from beckmann_nbo.config import (
    DATA_OUTPUT,
    FUNCTIONAL, BASIS, NPROC, MEM_GB,
    NBO_KEYWORDS, SOLVENT,
)
from beckmann_nbo.geometry import parse_standard_orientations
from beckmann_nbo.inputs import TEST_IDS, resolve_mol_name

OXIME_LABEL_RE = re.compile(r"\[oxime:\s*C(\d+)=N(\d+)-O(\d+)\]")


def gjf_sp(job_name: str, atoms: list, oxime_label: str, r_no: float) -> str:
    """Generate a single-point NBO7 .gjf for one intermediate scan geometry."""
    coord_block = "\n".join(
        f"{sym:<3}  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
        for sym, x, y, z in atoms
    )
    # No blank line between charge/mult and cartesian coordinates —
    # Gaussian treats a blank line as end of the molecule specification.
    return (
        f"%chk={job_name}.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} sp pop=nbo7read {SOLVENT}\n"
        f"\n"
        f"{job_name}  R(N-O)={r_no:.4f}A  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        f"{coord_block}\n"
        f"\n"
        f"$NBO {NBO_KEYWORDS} $END\n"
        f"\n\n"
    )


def extract_scan_sp(
    log_path: Path,
    out_dir: Path,
    ni: int,
    oi: int,
    oxime_label: str,
    mol_name: str,
    step: float = 0.1,
) -> list[Path]:
    """Extract converged geometries at R0+step, R0+2*step, R0+3*step from a scan log.

    Returns list of .gjf paths written.
    """
    print(f"Reading {log_path.name}...")
    lines = log_path.read_text().splitlines()

    so_blocks = parse_standard_orientations(lines)
    print(f"  {len(so_blocks)} Standard orientation blocks found")

    r0 = no_distance(so_blocks[0][1], ni, oi)
    print(f"  R0(N-O) = {r0:.4f} Å, step = {step} Å")

    written: list[Path] = []
    for step_n in (1, 2, 3):
        target = round(r0 + step_n * step, 4)
        candidates = [
            (idx, atoms) for idx, atoms in so_blocks
            if abs(no_distance(atoms, ni, oi) - target) < 1e-3
        ]
        if not candidates:
            print(f"  WARNING: no geometry found with R(N-O) ≈ {target:.4f} Å")
            continue

        # Take the LAST SO with this R value — Gaussian prints one extra SO at
        # the already-incremented R before starting the next optimization cycle.
        _, atoms = max(candidates, key=lambda x: x[0])
        r_actual = no_distance(atoms, ni, oi)
        job_name = f"{mol_name}_sp{step_n + 1}"

        text     = gjf_sp(job_name, atoms, oxime_label, r_actual)
        out_file = out_dir / f"{job_name}.gjf"
        out_file.write_text(text)
        written.append(out_file)
        print(f"  R0+{step_n * step:.1f} Å → R(N-O) = {r_actual:.4f} Å  →  {out_file.name}")

    return written


def oxime_indices_from_gjf(gjf_path: Path) -> tuple[int, int, str]:
    """Parse '[oxime: C{ci}=N{ni}-O{oi}]' out of a .gjf title line."""
    ci, ni, oi, _ = oxime_atom_map_from_gjf(gjf_path)
    return ni, oi, f"[oxime: C{ci}=N{ni}-O{oi}]"


def oxime_atom_map_from_gjf(gjf_path: Path) -> tuple[int, int, int, str]:
    """Parse '[oxime: C{ci}=N{ni}-O{oi}]' out of a .gjf title line, including ci."""
    match = OXIME_LABEL_RE.search(gjf_path.read_text())
    if not match:
        raise ValueError(f"{gjf_path}: no '[oxime: C#=N#-O#]' label found")
    ci, ni, oi = (int(g) for g in match.groups())
    return ci, ni, oi, f"[oxime: C{ci}=N{ni}-O{oi}]"


def process_molecule(mol: str) -> list[Path]:
    """Run extract_scan_sp() for one test-set molecule, e.g. 'mol_002_E'."""
    mol_dir  = DATA_OUTPUT / "dft_opt" / mol
    log_path = mol_dir / f"{mol}_scan.log"
    ni, oi, oxime_label = oxime_indices_from_gjf(mol_dir / f"{mol}_opt.gjf")

    written = extract_scan_sp(
        log_path    = log_path,
        out_dir     = mol_dir,
        ni          = ni,
        oi          = oi,
        oxime_label = oxime_label,
        mol_name    = mol,
    )

    print(f"\n{len(written)} files written to {mol_dir}")
    if written:
        mol_id = mol.split("_")[1]
        print("\nUpload and submit on Citadel:")
        print(f"  python scripts/dft/hpc_sync.py --mol {mol_id} upload")
        sp_names = " ".join(f.name for f in written)
        print(f"  # then submit: {sp_names}")
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mol", metavar="ID",
        help="Test-set molecule ID (e.g. 002). Default: all molecules in "
             "TEST_IDS with a completed _scan.log but no _sp2.gjf yet.",
    )
    args = parser.parse_args()

    if args.mol:
        mol_ids = [args.mol]
    else:
        mol_ids = sorted(TEST_IDS)

    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    for mol_id in mol_ids:
        mol = resolve_mol_name(mol_id, dft_opt_dir)
        if mol is None:
            print(f"-- mol_{mol_id.zfill(3)}: no directory, skipping")
            continue
        mol_dir = dft_opt_dir / mol
        scan_log = mol_dir / f"{mol}_scan.log"
        sp2_gjf  = mol_dir / f"{mol}_sp2.gjf"
        if not scan_log.exists():
            print(f"-- {mol}: no _scan.log yet, skipping")
            continue
        if sp2_gjf.exists() and not args.mol:
            print(f"-- {mol}: _sp2.gjf already exists, skipping (pass --mol to force)")
            continue

        print(f"\n== {mol} ==")
        process_molecule(mol)


if __name__ == "__main__":
    main()
