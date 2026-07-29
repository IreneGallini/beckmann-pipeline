"""
Weighted-family branch tracking for the CN/C1 avoided-crossing diagram --
Tetiana's reference method (SI_avoided_crossing_diagram.docx,
Weighted_Orbital_Character_Crossing_Handout.docx), validated against
example_scans/5_s1_Me.log .. 5_s4_Me.log (scripts/dft/validate_branch_tracking.py).

Deliberately separate from beckmann.dft.parse_cmo's compute_cn_extras()/
classify_crossing() -- that machinery independently re-picks a 'runner-up' MO
by single largest coefficient at each geometry, which is exactly the failure
mode this method replaces (handout Section 2): a canonical MO carrying real
target character can be outranked by an unrelated MO purely because of how
character happens to be split at that one geometry, and the same underlying
orbital can trade energy order with another entirely (a real avoided crossing)
between adjacent geometries -- tracking by MO index or single largest
coefficient hides that. This module only reuses the low-level CMO log parsing
(find_cmo_sections/parse_cmo_table/virtual_window/is_bond_between) from
parse_cmo.py, not any of its per-channel selection logic.

Family definitions per Detailed_Orbital_Character_Exchange_Handout.docx (which
supersedes SI_avoided_crossing_diagram.docx / Weighted_Orbital_Character_
Crossing_Handout.docx's Section 4 table -- Tetiana found her own script used
unsquared/unnormalized coefficients past the first scan point):

  C-C-side family (was "C1-side"): LP(n) C{c1}, RY*(n) C{c1}, BD*(n) C{c1}-X
                                    (any X, either side)
  C=N-side family (BROADENED from the old handout's narrow "just BD*(n)
    C{cn_c}-N{cn_n}"): LP(n) N{cn_n}, RY*(n) N{cn_n}, BD*(n) involving N{cn_n}
    with ANY other atom (not just C{cn_c}) -- e.g. BD*(1) N{cn_n}-O{oi} (the
    breaking N-O bond's own antibond) now counts. Same shape as the C-C-side
    family, just centered on N{cn_n} instead of C{c1} -- see
    _is_family_member().
  w_ref: a separate third quantity, NOT part of the two-branch candidate
    selection -- the fixed C{cn_c}-C{ref} antibond weight specifically (still
    matched via is_bond_between(), a fixed pair, unlike the two broadened
    families above). Maps directly onto this project's existing three NBO
    channels: C-C-side -> w17max, C=N-side -> wCNmax, w_ref -> w78max.

w_CC(i)/w_CN(i) = sum of squared coefficients of matching contribs in MO i.
w_target = w_CC + w_CN (NOT + w_ref -- every worked example in the handout's
own Section 6 computes f_CC/f_CN with only w_CC+w_CN in the denominator,
despite its prose formula section stating otherwise; implemented per the
worked math, not the prose).
f_CC = w_CC/w_target, f_CN = w_CN/w_target.
Coefficient signs are discarded for these (squared away) -- deliberately, per
the handout. track_branches() below does NOT do this for its own similarity
metric -- see its docstring for why the sign has to survive there.

Atom numbers (c1_atom/cn_c_atom/cn_n_atom/ref_atom) are parameterized, not
hardcoded to 1/7/17/8, but this module is only ever called with those values
(Tetiana's own reference-molecule numbering) in the current validation pass --
the benchmark set uses a different per-molecule numbering
(oxime_atom_map_from_gjf()/get_substituent_map()) and extending this method to
it is explicitly out of scope here.
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from beckmann.dft.parse_cmo import (
    find_cmo_sections, is_bond_between, parse_cmo_table, virtual_window,
)

HARTREE_TO_EV = 27.2114


def _is_family_member(label: str, element: str, atom_num: int) -> bool:
    """LP(n) {element}{atom_num} (not LP*), RY*(n) {element}{atom_num}, or BD*(n)
    involving {element}{atom_num} on either side -- the symmetric family-matching
    rule shared by both C-C-side (element='C', atom_num=c1_atom) and the
    broadened C=N-side (element='N', atom_num=cn_n_atom) per
    Detailed_Orbital_Character_Exchange_Handout.docx. Atom-mention regex is
    unanchored to position, so it matches both 'C1-X' and 'X-C1' orderings
    (checked directly against the reference logs: C-C-side matches always show
    up as 'C1-X' there, but the C=N-side match is always 'X-N17', i.e. both
    orderings genuinely occur and need to be handled)."""
    if not re.search(rf"{element}\s*{atom_num}(?!\d)", label):
        return False
    return bool(re.match(r"LP\s*\(", label) or re.match(r"RY\*\(", label) or re.match(r"BD\*\(", label))


def extract_family_weights(
    log_path: Path, c1_atom: int = 1, cn_c_atom: int = 7, cn_n_atom: int = 17,
    ref_atom: int = 8,
) -> list[dict]:
    """Per virtual MO in the 0.4 a.u.-above-LUMO window (virtual_window()), the
    C-C-side/C=N-side family weights and fractions, the separate w_ref
    (C{cn_c}-C{ref_atom} antibond) quantity, plus the MO's full raw contribs
    list (needed by track_branches() for continuity matching). Uses the LAST
    CMO section in the log if there's more than one (Stable=Opt seed vs. final
    pass -- same dedup convention as parse_cmo.parse_log()).

    Returns a list of {mo, energy, w_CC, w_CN, w_target, f_CC, f_CN, w_ref,
    contribs} dicts, one per virtual MO in the window, in the log's own
    (ascending energy) order. f_CC/f_CN are None when w_target == 0 (no family
    character at all in that MO). w_target = w_CC + w_CN only -- w_ref is
    tracked but deliberately excluded from it, see module docstring."""
    lines = Path(log_path).read_text().splitlines()
    starts = find_cmo_sections(lines)
    if not starts:
        return []
    table = parse_cmo_table(lines, starts[-1])
    vir = virtual_window(table)

    rows = []
    for mo in vir:
        w_cc = sum(c ** 2 for c, label in mo["contribs"] if _is_family_member(label, "C", c1_atom))
        w_cn = sum(c ** 2 for c, label in mo["contribs"] if _is_family_member(label, "N", cn_n_atom))
        w_ref = sum(
            c ** 2 for c, label in mo["contribs"]
            if re.match(r"BD\*\(", label) and is_bond_between(label, cn_c_atom, ref_atom)
        )
        w_target = w_cc + w_cn
        rows.append({
            "mo": mo["mo"], "energy": mo["energy"],
            "w_CC": w_cc, "w_CN": w_cn, "w_target": w_target,
            "f_CC": w_cc / w_target if w_target > 0 else None,
            "f_CN": w_cn / w_target if w_target > 0 else None,
            "w_ref": w_ref,
            "contribs": mo["contribs"],
        })
    return rows


def _signed_dot(mo_a: dict, mo_b: dict) -> float:
    """Similarity between two MOs' full (not just C1/CN-family) NBO contributions,
    for continuity tracking across geometries -- signed coefficients, matched by
    exact label text, summed only over labels present in BOTH MOs. Unnormalized
    (plain dot product, not cosine).

    Signed, not squared: family weights (above) deliberately discard sign because
    an antibond's phase doesn't matter for "how much C1/CN character is here."
    But continuity identity is a different question -- whether two MOs computed at
    adjacent geometries are 'the same' orbital -- and the relative sign PATTERN
    across a MO's several printed labels (which spectator bonds are in-phase vs.
    out-of-phase with the target antibond) is exactly the fingerprint that
    survives a small geometry step even when energy order doesn't. Verified by
    hand against the s1->s2 transition (the actual crossing step) before writing
    this: squared-coefficient dot product gives the WRONG branch assignment
    there, signed dot product gives the right one -- see
    beckmann/dft/branch_tracking.py's module history / the plan this was built
    from for the worked arithmetic."""
    b_by_label = dict(zip((label for _, label in mo_b["contribs"]), (c for c, _ in mo_b["contribs"])))
    return sum(
        coeff * b_by_label[label]
        for coeff, label in mo_a["contribs"]
        if label in b_by_label
    )


def track_branches(
    scan_logs: list[Path], c1_atom: int = 1, cn_c_atom: int = 7, cn_n_atom: int = 17,
) -> list[dict]:
    """Continuity-tracked branch identity (A/B) across an ordered list of scan-point
    logs, following handout Section 3's similarity-based reassignment instead of
    MO index or single-largest-coefficient.

    At each geometry, the candidate pair is the top-2 MOs by w_target AMONG THOSE
    WITH POSITIVE CANONICAL ENERGY, not the full virtual manifold. This positive-
    energy filter is an empirical finding, not something either handout states
    explicitly -- the handout's own virtual-orbital window is relative to the LUMO
    (0 < eps_k - eps_LUMO <= 0.4 a.u.), which virtual_window() already enforces on
    every MO passed in here; requiring the absolute eps_k > 0 on top of that is a
    stricter, separate condition. Without it, top-2-by-w_target over the full
    manifold lets a persistent, structurally unrelated MO (the N-O 'activation
    coordinate' antibond, always negative-energy in the s2-s4 range of the
    reference case) hijack a candidate slot and crowd out the true branch member --
    see the plan this was built from for the full characterization. WITH the
    filter, this reproduces the reference case's Section 5 table exactly at all 4
    points (branch identity AND A/B order) -- but that's one 4-point validation,
    not independent confirmation this is Tetiana's actual selection rule. Treat as
    provisional until confirmed with her, especially before extending to molecules
    where a candidate MO's energy sign relative to this cutoff isn't yet checked.

    At the first geometry, candidates become branch A (higher w_target) and B
    directly. At every later geometry, S_direct = S(A_old,cand1)+S(B_old,cand2) is
    compared against S_swapped = S(A_old,cand2)+S(B_old,cand1) (_signed_dot() for
    S), and whichever total is larger determines which candidate continues which
    branch -- this comparison logic is unchanged from before the positive-energy
    filter was added.

    Returns one row per scan point: {mo_A, E_A, f_CC_A, f_CN_A, mo_B, E_B,
    f_CC_B, f_CN_B}. E_* is in Hartree (convert with HARTREE_TO_EV for plotting)."""
    if not scan_logs:
        return []

    def _row(label, mo):
        return {
            f"mo_{label}": mo["mo"], f"E_{label}": mo["energy"],
            f"f_CC_{label}": mo["f_CC"], f"f_CN_{label}": mo["f_CN"],
        }

    results = []
    a_mo = b_mo = None
    for log_path in scan_logs:
        weights = extract_family_weights(log_path, c1_atom, cn_c_atom, cn_n_atom)
        candidates = [m for m in weights if m["energy"] > 0]
        by_target = sorted(candidates, key=lambda m: m["w_target"], reverse=True)
        cand1, cand2 = by_target[0], by_target[1]

        if a_mo is None:
            a_mo, b_mo = cand1, cand2
        else:
            s_direct = _signed_dot(a_mo, cand1) + _signed_dot(b_mo, cand2)
            s_swapped = _signed_dot(a_mo, cand2) + _signed_dot(b_mo, cand1)
            a_mo, b_mo = (cand1, cand2) if s_direct >= s_swapped else (cand2, cand1)

        results.append({**_row("A", a_mo), **_row("B", b_mo)})
    return results


def plot_branch_diagram(
    scan_labels: list[str] | list[float],
    branch_a: list[tuple[int, float]],
    branch_b: list[tuple[int, float]],
    branch_ref: list[tuple[int, float]] | None = None,
    title: str = "Crossing revealed by weighted orbital-character tracking",
    x_label: str = "Scan geometry",
    y_label: str = "Canonical MO eigenvalue (eV)",
    numeric_x: bool = False,
) -> Figure:
    """Branch energy vs. scan geometry, connected by branch label (A/B) -- NOT
    by MO number, since the whole point of this method is that the MO number
    carrying a branch's character can change while the branch itself is
    continuous. Each point is annotated with its carrier MO number directly on
    the plot.

    branch_a/branch_b: one (mo_index, energy) pair per scan point, same
    order/length as scan_labels -- e.g. track_branches()'s output, converted to
    eV via HARTREE_TO_EV or left in Hartree, caller's choice (see y_label).
    branch_ref: optional third series (e.g. w_ref's carrier MO) -- plotted
    dashed/gray, MO-annotated like the other two, but not one of the two
    tracked branches.
    numeric_x: when True, scan_labels are used as real numeric x-coordinates
    (e.g. R(N-O) in Angstroms, not necessarily evenly spaced) instead of
    evenly-spaced categorical positions with scan_labels as tick text."""
    fig, ax = plt.subplots(figsize=(10, 6))
    x = list(scan_labels) if numeric_x else list(range(len(scan_labels)))

    series_specs = [
        ("Weighted-character branch A", branch_a, "tab:blue", "-", "o"),
        ("Weighted-character branch B", branch_b, "tab:orange", "-", "o"),
    ]
    if branch_ref is not None:
        series_specs.append(("C7-C8 reference", branch_ref, "gray", "--", "s"))

    for label, series, color, linestyle, marker in series_specs:
        ys = [e for _, e in series]
        ax.plot(x, ys, marker=marker, linestyle=linestyle, color=color, label=label)
        for xi, (mo, e) in zip(x, series):
            ax.annotate(
                f"MO {mo}", (xi, e), textcoords="offset points", xytext=(0, 8),
                ha="center", fontsize=9,
            )

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if not numeric_x:
        ax.set_xticks(x)
        ax.set_xticklabels(scan_labels)
    ax.grid(True, color="lightgray")
    ax.set_axisbelow(True)

    # Headroom above the highest point/annotation so the upper-right legend has
    # empty space to sit in instead of covering the last points' "MO n" labels.
    ymin, ymax = ax.get_ylim()
    pad = (ymax - ymin) * 0.18 or 0.1
    ax.set_ylim(ymin, ymax + pad)

    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig
