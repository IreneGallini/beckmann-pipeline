"""
Gaussian-log-format-specific geometry parsing. The method-agnostic geometry
math (displace_leaving_group, no_distance, find_leaving_group, ATOMIC_SYMBOLS)
lives in beckmann_core.geometry -- this file only has the piece that's
inherently tied to Gaussian's own log output format.
"""
from beckmann_core.geometry import ATOMIC_SYMBOLS


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
