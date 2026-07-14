"""
Geometry primitives shared by beckmann.dft.scan and beckmann.dft.inputs.

Split out into its own module (rather than living in scan.py, where they
originated) so beckmann.dft.inputs can use displace_leaving_group() for the
rigid-scan architecture without a circular import (scan.py imports TEST_IDS
etc. from inputs.py).
"""
import math

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


def find_leaving_group(atoms: list, oi: int, bond_cutoff: float = 1.3) -> set[int]:
    """1-based indices of atoms bonded directly to O (other than through the
    N-O bond itself) -- e.g. the two H's of a protonated oxime's -OH2+ -- so
    they ride along with O when it's rigidly displaced. Determined purely by
    distance in the given geometry (O-H ~0.96-1.0 A comfortably clears the
    default cutoff; O-N, the bond being stretched, is excluded since it's
    always >=1.4 A at any scan point relevant here) -- no external bonding
    annotation (e.g. a $CHOOSE block) needed."""
    o_atom = atoms[oi - 1]
    group = {oi}
    for i, (_, x, y, z) in enumerate(atoms, start=1):
        if i == oi:
            continue
        dist = math.sqrt((x - o_atom[1]) ** 2 + (y - o_atom[2]) ** 2 + (z - o_atom[3]) ** 2)
        if dist < bond_cutoff:
            group.add(i)
    return group


def displace_leaving_group(atoms: list, ni: int, oi: int, delta: float) -> list:
    """Rigidly translate O (+ whatever's bonded to it, e.g. -OH2+'s two H's)
    by `delta` Angstroms along the N->O unit vector; every other atom is
    untouched. Used to build each rigid-scan point's starting geometry
    independently from a single fixed base structure (Stage 1's converged
    geometry), rather than chaining from the previous point or extracting
    from an internal Gaussian scan walk -- see the rigid-scan architecture
    notes in JOB_ISSUES.md / Notes.md for why.
    """
    leaving_group = find_leaving_group(atoms, oi)
    n_atom, o_atom = atoms[ni - 1], atoms[oi - 1]
    vec = (o_atom[1] - n_atom[1], o_atom[2] - n_atom[2], o_atom[3] - n_atom[3])
    vlen = math.sqrt(sum(c * c for c in vec))
    unit = tuple(c / vlen for c in vec)
    shift = tuple(c * delta for c in unit)
    return [
        (sym, x + shift[0], y + shift[1], z + shift[2]) if i in leaving_group
        else (sym, x, y, z)
        for i, (sym, x, y, z) in enumerate(atoms, start=1)
    ]
