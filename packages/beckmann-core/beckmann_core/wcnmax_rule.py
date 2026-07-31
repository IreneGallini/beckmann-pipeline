"""
The wCNmax-minimum R/F rule: does an interior wCNmax(R) minimum predict
Beckmann rearrangement? Operates on a generic row shape (mol/stage/channel/
R_NO/weight/MO_index/epsilon_i_star) that doesn't care whether the rows came
from NBO7 or PySCF -- this is genuinely shared logic, not NBO-specific
(confirmed: beckmann-pyscf's own scan-row adapter feeds PySCF-computed rows
through find_wcnmax_minimum() unmodified). Moved here verbatim from
beckmann/dft/descriptors.py + beckmann/dft/wcnmax_rule.py; the NBO-specific
half of the old descriptors.py (get_substituent_map, compute_slopes,
load_series -- Psi/Lambda are NBO-only quantities) stays in beckmann-nbo.

Rule (as specified, not derived): if a molecule's wCNmax(R) series has a
genuine interior local minimum, predict 'R' (rearrangement); otherwise
predict 'F' (fragmentation).
"""


def resolve_series(by_stage: dict) -> list[dict]:
    """Pick the R(N-O) series from a {stage: row} map: 'nbo' (R0) followed by
    every 'scan_N' stage present, sorted numerically by N. Dynamic rather
    than a fixed-length list so it naturally covers however many stretched
    points a molecule's series actually contains."""
    series = []
    if "nbo" in by_stage:
        series.append(by_stage["nbo"])
    scan_stages = sorted(
        (s for s in by_stage if s.startswith("scan_")),
        key=lambda s: int(s.split("_")[1]),
    )
    series.extend(by_stage[s] for s in scan_stages)
    return series


def find_wcnmax_extremum(mol: str, extraction_rows: list[dict]) -> dict | None:
    """R_star/w_star/MO_index/epsilon_i_star at the MOST PROMINENT interior
    wCNmax extremum -- min OR max, largest |depth| -- if any. See
    find_wcnmax_minimum() below for the minimum-only filter the wCNmax rule
    actually needs.

    Scans every interior point and keeps the one with the largest |depth|
    (not just the first found), since a small local wobble can sit right
    before the real, much deeper extremum in a finer-resolution series.
    """
    by_stage = {
        r["stage"]: r for r in extraction_rows
        if r["mol"] == mol and r["channel"] == "cn" and r["weight"] not in (None, "", "None")
    }
    rows = resolve_series(by_stage)
    pts = [(float(r["R_NO"]), float(r["weight"]), r["MO_index"], r["epsilon_i_star"]) for r in rows]
    if len(pts) < 3:
        return None
    pts.sort(key=lambda p: p[0])
    best = None
    for i in range(1, len(pts) - 1):
        _, w_prev, _, _ = pts[i - 1]
        r_cur, w_cur, mo_cur, eps_cur = pts[i]
        _, w_next, _, _ = pts[i + 1]
        if (w_cur < w_prev and w_cur < w_next) or (w_cur > w_prev and w_cur > w_next):
            # depth = how far w_cur sits below (positive) or above (negative) the
            # midpoint of its two neighbors -- the yes/no extremum flag alone can't
            # tell a barely-there wobble from a deep collapse.
            depth = (w_prev + w_next) / 2 - w_cur
            if best is None or abs(depth) > abs(best["depth"]):
                best = {
                    "R_star": r_cur, "w_star": w_cur, "MO_index": mo_cur,
                    "epsilon_i_star": float(eps_cur) if eps_cur not in (None, "", "None") else None,
                    "depth": depth,
                }
    return best


def find_wcnmax_minimum(mol: str, extraction_rows: list[dict]) -> dict | None:
    """Same as find_wcnmax_extremum(), but only a genuine interior MINIMUM
    counts (depth > 0 -- w_cur sits below both neighbors) -- a local maximum
    (depth < 0) returns None here. This is the specific signature
    predict_from_wcnmax() uses: an interior minimum predicts rearrangement
    ('R'), its absence predicts fragmentation ('F')."""
    extremum = find_wcnmax_extremum(mol, extraction_rows)
    if extremum is None or extremum["depth"] <= 0:
        return None
    return extremum


def predict_from_wcnmax(minimum: dict | None) -> str:
    """'R' if a genuine interior wCNmax minimum was found, else 'F'."""
    return "R" if minimum is not None else "F"
