"""
Task 1 follow-up: a per-atom-pair local construction, closer to NBO's actual Step 4
(see beckmann_alt/pyscf_livvo.py's module docstring for the reference-material
breakdown) than LIVVO's single global SVD.

Motivation (see Notes_open_source_alt.md's "Follow-up: does a better channel-
identification metric fix it?"): trying a better SCORING rule over LIVVO's fixed
~29-30-orbital set did not resolve the channel collisions found in the first pass --
both the mol_002 and 5_s0_Me cases kept mapping multiple different target antibonds
onto the *same* orbital, because all candidates were being drawn from one fixed,
global, non-iterative set. The diagnosis there: no scoring rule over a fixed set can
manufacture a better-localized orbital if the set doesn't contain one.

This module tests the actual next step implied by that diagnosis: instead of picking
from a pre-built global set, build a FRESH, LOCAL subspace for each atom pair
individually, direct from that pair's own 2x2-block-of-atoms slice of the density
matrix -- structurally the same operation as NBO's own per-atom-pair block
diagonalization (Step 4), just for a single requested pair rather than NBO's full
iterative multi-pair deflation across the whole molecule (that fuller construction,
which subtracts each accepted bond's occupancy before moving to the next atom pair, is
a meaningfully bigger undertaking and not attempted here -- see caveat below).

Construction:
  1. Build IAOs on the occupied space (pyscf.lo.iao.iao) and Löwdin-orthogonalize them
     (pyscf.lo.orth.vec_lowdin) -- same as the first step of pyscf.lo.vvo.vvo, but we
     keep every IAO column (the full minimal-basis-sized set, not just the
     valence-virtual-sized rotation LIVVO reduces to) and their per-atom labels
     (pyscf.lo.iao.reference_mol's ao_labels, which are in 1:1 order correspondence
     with the IAO columns).
  2. Transform the full AO-basis density matrix into this orthonormal IAO basis:
     D_IAO = C_iao^T S D_AO S C_iao (standard density-matrix basis change under an
     S-orthonormal transformation).
  3. For the two target atoms, slice out just their IAOs' rows/columns from D_IAO --
     a small, LOCAL block (typically ~9-18 dimensional for two heavy atoms in
     6-311+G(d,p)'s minao-sized IAO space), fresh for every atom pair, not shared
     with any other channel's search.
  4. Diagonalize that local block directly (a legitimate, standard NBO-style
     approximation -- NBO itself diagonalizes each 2-center block as if it were an
     isolated 2-atom system). Every eigenvector below ANTIBOND_OCC_THRESHOLD is kept as
     an antibond CANDIDATE, not just the single lowest-occupation one -- see that
     constant's comment: for a double-bonded pair (C=N), the local block has TWO
     low-occupation eigenvectors (sigma*-like and pi*-like), and it is NOT always the
     lowest-occupation one that actually carries the physically dominant character.
     The highest-occupation eigenvector is the bond-like (BD) combination, kept only as
     a diagnostic (should read ~2.0 for a real localized bond -- and does, in both test
     cases, in every channel).
  5. Rotate each candidate back into the full AO basis, project every canonical virtual
     MO onto it (beckmann_alt.pyscf_livvo.project_virtuals_onto_livvo -- same projection
     math regardless of where the target vector came from), and keep whichever candidate
     gives the largest projected weight -- mirroring the real wX^max definition's own
     "both BD*(1) and BD*(2) are eligible, whichever gives the larger squared
     coefficient wins" rule (Notes.md) exactly, just applied to locally-built candidates
     instead of NBO's own pre-labeled BD*(1)/BD*(2) pair.

Caveat -- what run_case()/compute_channel_descriptor() (the ORIGINAL, non-deflated
functions above) do NOT reproduce: NBO's real algorithm is iterative across the WHOLE
molecule -- once a pair's bond/antibond is accepted, its occupancy is subtracted from
the density matrix before searching the next pair, so two different target pairs
sharing an atom (e.g. our cn=(ci,ni) and w17=(ci,c_aryl), which both include ci) don't
end up drawing on the same undiminished density at ci. Those functions search each pair
completely independently with no deflation between them, and this caused a real,
measured problem: w78 landed on the exact same canonical MO as cn in both reference
cases (see Notes_open_source_alt.md's "Second follow-up"). The "Third follow-up:
iterative deflation" section below (run_case_deflated/compute_channel_descriptor_
deflating/deflate_bond) implements deflation between these three specific channel
searches to test whether that fixes it -- see Notes_open_source_alt.md's "Third
follow-up" section for whether it did. The ORIGINAL non-deflated functions are kept
as-is, unmodified in behavior, as the direct before/after baseline.
"""
import numpy as np
from pyscf import gto, lo

from beckmann_alt.geometry import load_case
from beckmann_alt.pyscf_livvo import build_mol, run_scf, project_virtuals_onto_livvo


def build_local_iaos(mf) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Full (not valence-virtual-reduced) Löwdin-orthogonalized IAO set, plus the
    0-based atom index owning each IAO column. Returns (iaos_orth, s, atom_of_iao)."""
    mol = mf.mol
    s = mol.intor_symmetric("int1e_ovlp")
    nocc = mol.nelectron // 2
    orbocc = mf.mo_coeff[:, :nocc]

    iaos = lo.iao.iao(mol, orbocc)
    iaos_orth = lo.orth.vec_lowdin(iaos, s)

    ref_mol = lo.iao.reference_mol(mol)
    atom_of_iao = np.array([label[0] for label in ref_mol.ao_labels(fmt=False)])
    return iaos_orth, s, atom_of_iao


def pair_density_matrix(mf) -> np.ndarray:
    """AO-basis density matrix from cached/converged mo_coeff (closed-shell: 2 * occ @ occ.T).
    Avoids requiring a live mf.make_rdm1() so cached (mo_coeff-only) results work too."""
    nocc = mf.mol.nelectron // 2
    orbocc = mf.mo_coeff[:, :nocc]
    return 2.0 * orbocc @ orbocc.T


ANTIBOND_OCC_THRESHOLD = 1.0  # local eigenvectors below this occupation are treated as
# antibond CANDIDATES, not just the single lowest one. Found necessary empirically: for
# a double-bonded pair (e.g. C=N), the local block has TWO low-occupation eigenvectors
# (sigma* and pi*-like), not one -- mol_002's cn=(C11,N12) block gives occupations
# [0.016, 0.219, 0.568, 1.02, 1.16, 1.98, ...], and it's the SECOND-lowest (0.219), not
# the lowest (0.016), that actually projects onto the virtual manifold close to NBO7's
# trusted wCNmax (0.460 vs trusted 0.436, MO47 vs trusted MO48 -- taking only the single
# lowest eigenvector gave 0.087, roughly 5x too small). This mirrors the real wX^max
# definition directly: "both the BD*(1) (sigma) and BD*(2) (pi) components ... whichever
# gives the larger squared coefficient wins" (Notes.md) -- we don't know in advance
# which local eigenvector plays that BD*(1)/BD*(2) role, so every sub-threshold
# candidate is projected and the best one is kept, same as NBO7's own convention.

BOND_OCC_THRESHOLD = 1.9  # local eigenvectors AT/ABOVE this occupation are treated as
# genuine, fully-formed bonds -- safe to fully deflate (subtract) from the working
# density matrix. Deliberately NOT the same as ANTIBOND_OCC_THRESHOLD (1.0), despite
# first trying that symmetric choice: it caused severe PSD violations (local
# occupations down to -0.56/-0.60, far past PSD_NEGATIVE_THRESHOLD) because every local
# block has an "ambiguous middle" cluster of eigenvectors (occupation roughly 0.5-1.4)
# that are NOT clean, fully-localized 2-electron bonds -- fully subtracting them as if
# they were over-deflates. Both molecules, all three channels, show the exact same
# pattern: the ambiguous middle always tops out below ~1.4, and genuine bonds always
# cluster at ~1.96-2.00 -- e.g. mol_002's w17 block: [0.03, 0.53, 0.84, 0.94, 1.00,
# 1.17, 1.37, 1.96, 2.00, 2.00]. 1.9 sits cleanly in that gap in every case observed.
# Eigenvectors strictly between ANTIBOND_OCC_THRESHOLD and BOND_OCC_THRESHOLD are
# neither an antibond candidate nor deflated -- left untouched, since neither
# interpretation is well-supported for them.


def local_pair_antibonds(
    mf, s: np.ndarray, iaos_orth: np.ndarray, atom_of_iao: np.ndarray, atom_a: int, atom_b: int,
    dm_ao: np.ndarray | None = None,
) -> dict:
    """Build every antibond-like CANDIDATE local combination for one atom pair (1-based
    atom_a/atom_b) -- every local eigenvector below ANTIBOND_OCC_THRESHOLD, not just the
    single lowest (see that constant's comment for why). Returns the AO-basis candidate
    vectors plus diagnostics (occupations of every local eigenvector, so a caller can
    see how cleanly separated the "antibond" cluster is from the rest of the local
    spectrum).

    dm_ao defaults to pair_density_matrix(mf) (the original, undeflated density) --
    passing a deflated density (see deflate_bond) is what lets a later channel's search
    see an earlier channel's accepted bond(s) as already "used up." min(all_occupations)
    is worth checking by callers after passing a deflated dm_ao: a deflated bond
    vector's overlap onto a later pair's IAO columns is a partial (non-unit) vector, so
    the resulting rank-1 downdate can in principle push a local eigenvalue negative --
    see Notes_open_source_alt.md's "Third follow-up" for whether/when this was observed.
    """
    a0, b0 = atom_a - 1, atom_b - 1
    local_idx = np.where((atom_of_iao == a0) | (atom_of_iao == b0))[0]
    if len(local_idx) < 2:
        raise ValueError(f"atoms {atom_a}/{atom_b}: fewer than 2 IAOs found -- can't form a pair block")

    if dm_ao is None:
        dm_ao = pair_density_matrix(mf)
    local_iaos = iaos_orth[:, local_idx]
    dm_local = local_iaos.T @ s @ dm_ao @ s @ local_iaos  # (n_local, n_local), IAO basis is S-orthonormal

    occupations, evecs = np.linalg.eigh(dm_local)  # ascending occupation order
    n_candidates = max(1, int(np.searchsorted(occupations, ANTIBOND_OCC_THRESHOLD)))
    candidates_ao = [local_iaos @ evecs[:, k] for k in range(n_candidates)]
    # Bond side: every eigenvector AT/ABOVE BOND_OCC_THRESHOLD (NOT the same threshold
    # as the antibond side -- see that constant's comment: reusing ANTIBOND_OCC_THRESHOLD
    # here caused severe over-deflation by treating ambiguous mid-spectrum eigenvectors
    # as if they were clean 2-electron bonds). Still deflates every genuine bond found,
    # not just the single highest-occupation one -- a double-bonded pair like C=N has
    # TWO near-2.0-occupation eigenvectors (sigma- and pi-like), and deflating only the
    # top one would leave the other undeflated.
    bond_start = int(np.searchsorted(occupations, BOND_OCC_THRESHOLD))
    bond_candidates_ao = [local_iaos @ evecs[:, k] for k in range(bond_start, len(occupations))]

    return {
        "candidates_ao": candidates_ao,
        "candidate_occupations": occupations[:n_candidates],
        "bond_candidates_ao": bond_candidates_ao,
        "bond_candidate_occupations": occupations[bond_start:],
        "bond_occupation": occupations[-1],
        "n_local_iaos": len(local_idx),
        "all_occupations": occupations,
    }


def compute_channel_descriptor(mf, s: np.ndarray, iaos_orth: np.ndarray, atom_of_iao: np.ndarray, atom_a: int, atom_b: int) -> dict:
    pair = local_pair_antibonds(mf, s, iaos_orth, atom_of_iao, atom_a, atom_b)
    best = None
    for occ, vec in zip(pair["candidate_occupations"], pair["candidates_ao"]):
        proj = project_virtuals_onto_livvo(mf, s, vec)
        if best is None or proj["weight"] > best["weight"]:
            best = {**proj, "antibond_occupation": occ}
    return {
        "pair": pair,
        "wmax": best["weight"], "mo_index": best["mo_index"],
        "epsilon": best["epsilon"], "coefficient": best["coefficient"],
        "antibond_occupation": best["antibond_occupation"],
    }


def run_case(name: str) -> dict:
    case = load_case(name)
    mol = build_mol(case)
    mf = run_scf(mol)

    iaos_orth, s, atom_of_iao = build_local_iaos(mf)

    cn  = compute_channel_descriptor(mf, s, iaos_orth, atom_of_iao, case["ci"], case["ni"])
    c17 = compute_channel_descriptor(mf, s, iaos_orth, atom_of_iao, case["ci"], case["c_aryl"])
    c78 = compute_channel_descriptor(mf, s, iaos_orth, atom_of_iao, case["ci"], case["c_alkyl"])

    return {"case": name, "basis_note": case["basis_note"], "cn": cn, "w17": c17, "w78": c78}


# ── Third follow-up: iterative deflation (see Notes_open_source_alt.md) ─────────────
#
# Motivation: run_case()/compute_channel_descriptor() above search all three channels
# against the SAME undeflated density matrix. Since cn=(ci,ni), w17=(ci,c_aryl), and
# w78=(ci,c_alkyl) all share atom ci, this lets w78's search rediscover density already
# "claimed" by cn's accepted bond -- diagnosed as the cause of w78 landing on the exact
# same canonical MO as cn in both reference cases. NBO's real algorithm avoids this by
# deflating (subtracting) each accepted bond's density before searching the next atom
# pair. This section implements that for these three channels specifically -- not a
# full molecule-wide iterative NBO search, just deflation between these three searches.

PSD_NEGATIVE_THRESHOLD = -0.01  # a post-deflation local occupation below this is
# flagged as a real PSD violation, not float noise -- see local_pair_antibonds'
# docstring: a deflated bond vector's overlap onto a later pair's IAO columns is a
# partial (non-unit) vector, so the resulting rank-1 downdate can in principle drive a
# later channel's local eigenvalue negative. Checked and reported explicitly rather
# than silently letting a negative-occupation "antibond candidate" compete in the
# projection step.

CHANNEL_ATOM_KEYS = {"cn": ("ci", "ni"), "w17": ("ci", "c_aryl"), "w78": ("ci", "c_alkyl")}


def deflate_bond(dm_ao: np.ndarray, s: np.ndarray, bond_vec_ao: np.ndarray, occupation: float) -> np.ndarray:
    """Subtract one accepted bond's density contribution: dm_ao - occupation * |v><v|.

    Uses the MEASURED occupation (not a hardcoded 2.0) -- exact for the trace/electron-
    count identity trace((dm-occ*vv^T) @ S) == trace(dm @ S) - occ (since bond_vec_ao is
    S-normalized, v^T @ S @ v == 1), asserted below as a cheap correctness check on the
    outer-product bookkeeping itself. Does not, and cannot, catch the separate local-PSD
    risk (see PSD_NEGATIVE_THRESHOLD) -- that's a property of the *local* block for a
    different atom pair, not of this global trace identity.
    """
    new_dm_ao = dm_ao - occupation * np.outer(bond_vec_ao, bond_vec_ao)
    expected = np.trace(dm_ao @ s) - occupation
    actual = np.trace(new_dm_ao @ s)
    assert abs(actual - expected) < 1e-8, f"deflate_bond trace mismatch: {actual} vs {expected}"
    return new_dm_ao


def compute_channel_descriptor_deflating(
    mf, s: np.ndarray, iaos_orth: np.ndarray, atom_of_iao: np.ndarray, atom_a: int, atom_b: int, dm_ao: np.ndarray,
) -> tuple[dict, np.ndarray]:
    """Like compute_channel_descriptor, but takes/returns the working (possibly already
    deflated) density matrix explicitly, and deflates every above-threshold bond
    candidate from THIS pair's own local block before returning -- so the next call in
    a sequence sees this channel's accepted bond(s) as already used up."""
    pair = local_pair_antibonds(mf, s, iaos_orth, atom_of_iao, atom_a, atom_b, dm_ao=dm_ao)

    min_occ = pair["all_occupations"].min()
    psd_violation = min_occ < PSD_NEGATIVE_THRESHOLD

    best = None
    for occ, vec in zip(pair["candidate_occupations"], pair["candidates_ao"]):
        proj = project_virtuals_onto_livvo(mf, s, vec)
        if best is None or proj["weight"] > best["weight"]:
            best = {**proj, "antibond_occupation": occ}

    for bond_occ, bond_vec in zip(pair["bond_candidate_occupations"], pair["bond_candidates_ao"]):
        dm_ao = deflate_bond(dm_ao, s, bond_vec, bond_occ)

    result = {
        "pair": pair,
        "wmax": best["weight"], "mo_index": best["mo_index"],
        "epsilon": best["epsilon"], "coefficient": best["coefficient"],
        "antibond_occupation": best["antibond_occupation"],
        "min_local_occupation": min_occ,
        "psd_violation": psd_violation,
    }
    return result, dm_ao


def _channel_order_by_occupation(mf, s, iaos_orth, atom_of_iao, case: dict) -> tuple[str, ...]:
    """First pass (undeflated): each channel's own max bond-candidate occupation,
    sorted descending. Printed by callers so a near-tie (expected here -- see
    Notes_open_source_alt.md, all three channels' bond occupations were ~2.0000 in the
    non-deflated data) is visible rather than presented as a meaningful signal."""
    dm_ao = pair_density_matrix(mf)
    occs = {}
    for label, (a_key, b_key) in CHANNEL_ATOM_KEYS.items():
        pair = local_pair_antibonds(mf, s, iaos_orth, atom_of_iao, case[a_key], case[b_key], dm_ao=dm_ao)
        occs[label] = pair["bond_candidate_occupations"].max()
    order = tuple(sorted(occs, key=lambda k: -occs[k]))
    return order, occs


def run_case_deflated(name: str, order: tuple[str, ...] | None = None) -> dict:
    """order=None -> determine order by decreasing bond occupation (see
    _channel_order_by_occupation); pass an explicit tuple (e.g. ("cn","w17","w78")) to
    compare against a fixed order instead."""
    case = load_case(name)
    mol = build_mol(case)
    mf = run_scf(mol)
    iaos_orth, s, atom_of_iao = build_local_iaos(mf)

    occ_order, occ_values = _channel_order_by_occupation(mf, s, iaos_orth, atom_of_iao, case)
    if order is None:
        order = occ_order

    dm_ao = pair_density_matrix(mf)
    results = {}
    for label in order:
        a_key, b_key = CHANNEL_ATOM_KEYS[label]
        result, dm_ao = compute_channel_descriptor_deflating(
            mf, s, iaos_orth, atom_of_iao, case[a_key], case[b_key], dm_ao,
        )
        results[label] = result

    return {
        "case": name, "basis_note": case["basis_note"], "order": order,
        "occupation_order": occ_order, "occupation_values": occ_values,
        **results,
    }


def _print_result(result: dict, title: str) -> None:
    print(f"\n=== {result['case']} ({result['basis_note']}) -- {title} ===")
    for label, key in [("wCNmax", "cn"), ("w17max", "w17"), ("w78max", "w78")]:
        r = result[key]
        p = r["pair"]
        psd_note = ""
        if "psd_violation" in r:
            flag = "*** PSD VIOLATION ***" if r["psd_violation"] else "ok"
            psd_note = f", min_local_occ={r['min_local_occupation']:.4f} [{flag}]"
        print(
            f"  {label}: wmax={r['wmax']:.4f}  MO={r['mo_index']}  "
            f"[local block: {p['n_local_iaos']} IAOs, {len(p['candidates_ao'])} antibond candidate(s), "
            f"winning_occ={r['antibond_occupation']:.4f}, bond_occ={p['bond_occupation']:.4f}{psd_note}]"
        )


if __name__ == "__main__":
    import sys
    name = sys.argv[1] if len(sys.argv) > 1 else "mol_002"

    baseline = run_case(name)
    _print_result(baseline, "baseline, no deflation")

    deflated_occ_order = run_case_deflated(name, order=None)
    print(f"\n  (occupation-order channels: {deflated_occ_order['occupation_order']}, "
          f"values: { {k: round(v, 4) for k, v in deflated_occ_order['occupation_values'].items()} })")
    _print_result(deflated_occ_order, f"deflated, decreasing-occupation order {deflated_occ_order['order']}")

    deflated_cn_first = run_case_deflated(name, order=("cn", "w17", "w78"))
    _print_result(deflated_cn_first, f"deflated, cn-first order {deflated_cn_first['order']}")
