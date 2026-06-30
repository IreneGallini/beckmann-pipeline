"""
Extract converged intermediate geometries from a Gaussian scan log and write
single-point NBO input files for the missing scan points.

Usage:
    python scripts/dft/extract_scan_sp.py

For mol_002_E, the scan ran NBO only at R0 and R0+0.4 A.  This script
extracts the converged geometries at R0+0.1, R0+0.2, R0+0.3 and creates
three single-point .gjf files for upload to Citadel.

Output: data/output/dft_opt/mol_002_E/{name}_sp{N}.gjf  (N = 2, 3, 4)
"""

import re
import math
from pathlib import Path

FUNCTIONAL   = "wB97XD"
BASIS        = "6-311+G(d,p)"
NPROC        = 8
MEM_GB       = 16
CHARGE       = 1
MULTIPLICITY = 1
NBO_KEYWORDS = "E2PERT BNDIDX NBOSUM"

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
            j = i + 5  # skip: dashes, column header, dashes
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


def find_scan_point_starts(lines: list[str]) -> dict[int, int]:
    """Return {scan_point_number: line_of_first_step} for each scan point."""
    starts: dict[int, int] = {}
    for i, line in enumerate(lines):
        m = re.search(r"Step number\s+1 out of.*scan point\s+(\d+) out of", line)
        if m:
            sp = int(m.group(1))
            if sp not in starts:
                starts[sp] = i
    return starts


def no_distance(atoms: list, i: int, j: int) -> float:
    a, b = atoms[i - 1], atoms[j - 1]
    return math.sqrt((a[1]-b[1])**2 + (a[2]-b[2])**2 + (a[3]-b[3])**2)


def gjf_sp(job_name: str, atoms: list, oxime_label: str, r_no: float) -> str:
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
        f"$NBO {NBO_KEYWORDS} $END\n"
        f"\n\n"
    )


def main() -> None:
    root     = Path(__file__).parent.parent.parent
    mol      = "mol_002_E"
    log_path = root / "data/output/dft_opt" / mol / f"{mol}_scan.log"
    out_dir  = root / "data/output/dft_opt" / mol

    NI, OI      = 12, 13           # 1-based Gaussian atom indices (N and O)
    oxime_label = "[oxime: C11=N12-O13]"

    print(f"Reading {log_path.name}...")
    lines = log_path.read_text().splitlines()

    so_blocks = parse_standard_orientations(lines)
    print(f"  {len(so_blocks)} Standard orientation blocks found")

    # Determine R0 from the first SO block (the DFT equilibrium geometry).
    r0 = no_distance(so_blocks[0][1], NI, OI)
    step = 0.1
    print(f"  R0(N-O) = {r0:.4f} Å, step = {step} Å")

    # For each intermediate point, pick the LAST SO block whose N-O distance
    # matches the target within 1e-3 Å.  Gaussian prints one extra SO at each
    # constraint increment (already showing the NEW distance) before starting
    # the next optimisation; taking the LAST SO with the CORRECT distance
    # avoids picking that incremented-but-not-yet-relaxed geometry.
    written = []
    for step_n in (1, 2, 3):          # R0+0.1, R0+0.2, R0+0.3
        target = round(r0 + step_n * step, 4)
        candidates = [
            (idx, atoms) for idx, atoms in so_blocks
            if abs(no_distance(atoms, NI, OI) - target) < 1e-3
        ]
        if not candidates:
            print(f"  WARNING: no geometry found with R(N-O) ≈ {target:.4f} Å")
            continue

        _, atoms = max(candidates, key=lambda x: x[0])
        r_actual = no_distance(atoms, NI, OI)
        job_name = f"{mol}_sp{step_n + 1}"   # sp2, sp3, sp4 to match scan point labels

        text     = gjf_sp(job_name, atoms, oxime_label, r_actual)
        out_file = out_dir / f"{job_name}.gjf"
        out_file.write_text(text)
        written.append((step_n, r_actual, out_file.name))
        print(f"  R0+{step_n * step:.1f} Å → R(N-O) = {r_actual:.4f} Å  →  {out_file.name}")

    print(f"\n{len(written)} files written to {out_dir}")
    print("\nUpload and submit on Citadel:")
    print("  python scripts/dft/hpc_sync.py --mol 002 upload")
    sp_names = " ".join(f"{mol}_sp{n+1}.gjf" for n, _, _ in written)
    print(f"  ssh citadel 'cd <HPC_REMOTE_DIR>/{mol} &&")
    print(f"    for f in {sp_names}; do")
    print(f"      nohup /opt/g16/g16 < $f > ${{f%.gjf}}.log 2>&1 &")
    print(f"    done'")


if __name__ == "__main__":
    main()
