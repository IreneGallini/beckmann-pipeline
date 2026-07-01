"""
Extract converged intermediate geometries from a Gaussian scan log and write
single-point NBO input files for the missing scan points.

For mol_002_E, the scan ran NBO only at R0 and R0+0.4 Å.  This module
extracts the converged geometries at R0+0.1, R0+0.2, R0+0.3 and creates
three single-point .gjf files for upload to Citadel.

Output: data/output/dft_opt/mol_002_E/{name}_sp{N}.gjf  (N = 2, 3, 4)
"""
import math
import re
from pathlib import Path

from beckmann.config import (
    DATA_OUTPUT,
    FUNCTIONAL, BASIS, NPROC, MEM_GB, CHARGE, MULTIPLICITY,
    NBO_KEYWORDS_EQ,
)

ATOMIC_SYMBOLS = {
    1: "H",  6: "C",  7: "N",  8: "O",  9: "F",
    16: "S", 17: "Cl", 35: "Br",
}


def parse_standard_orientations(lines: list[str]) -> list[tuple[int, list]]:
    """Return [(header_line_idx, atoms)] for every Standard orientation block."""
    blocks = []
    i = 0
    while i < len(lines):
        if "Standard orientation:" in lines[i]:
            j = i + 5
            atoms = []
            while j < len(lines) and "---" not in lines[j]:
                parts = lines[j].split()
                if len(parts) == 6:
                    atomic_num = int(parts[1])
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    sym = ATOMIC_SYMBOLS.get(atomic_num, f"X{atomic_num}")
                    atoms.append((sym, x, y, z))
                j += 1
            if atoms:
                blocks.append((i, atoms))
            i = j
        else:
            i += 1
    return blocks


def no_distance(atoms: list, i: int, j: int) -> float:
    """Euclidean distance between 1-based atom indices i and j."""
    a, b = atoms[i - 1], atoms[j - 1]
    return math.sqrt((a[1]-b[1])**2 + (a[2]-b[2])**2 + (a[3]-b[3])**2)


def gjf_sp(job_name: str, atoms: list, oxime_label: str, r_no: float) -> str:
    """Generate a single-point NBO .gjf for one intermediate scan geometry."""
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
        f"#p {FUNCTIONAL}/{BASIS} sp pop=nboread\n"
        f"\n"
        f"{job_name}  R(N-O)={r_no:.4f}A  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        f"{coord_block}\n"
        f"\n"
        f"$NBO {NBO_KEYWORDS_EQ} $END\n"
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


def main() -> None:
    mol      = "mol_002_E"
    mol_dir  = DATA_OUTPUT / "dft_opt" / mol
    log_path = mol_dir / f"{mol}_scan.log"

    NI, OI      = 12, 13
    oxime_label = "[oxime: C11=N12-O13]"

    written = extract_scan_sp(
        log_path    = log_path,
        out_dir     = mol_dir,
        ni          = NI,
        oi          = OI,
        oxime_label = oxime_label,
        mol_name    = mol,
    )

    print(f"\n{len(written)} files written to {mol_dir}")
    print("\nUpload and submit on Citadel:")
    print("  python scripts/dft/hpc_sync.py --mol 002 upload")
    sp_names = " ".join(f.name for f in written)
    print(f"  # then submit: {sp_names}")


if __name__ == "__main__":
    main()
