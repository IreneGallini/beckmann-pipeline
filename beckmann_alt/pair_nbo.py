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

Caveat -- what this does NOT reproduce: NBO's real algorithm is iterative across the
WHOLE molecule -- once a pair's bond/antibond is accepted, its occupancy is subtracted
from the density matrix before searching the next pair, so two different target pairs
sharing an atom (e.g. our cn=(ci,ni) and w17=(ci,c_aryl), which both include ci) don't
end up drawing on the same undiminished density at ci. This module searches each pair
completely independently with no deflation between them -- closer to NBO's Step 4 than
LIVVO, but still a simplification of it, and channel collisions from double-counting
ci's density across pairs are still possible in principle. Reported plainly in
Notes_open_source_alt.md's results, not assumed away.
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


def local_pair_antibonds(
    mf, s: np.ndarray, iaos_orth: np.ndarray, atom_of_iao: np.ndarray, atom_a: int, atom_b: int,
) -> dict:
    """Build every antibond-like CANDIDATE local combination for one atom pair (1-based
    atom_a/atom_b) -- every local eigenvector below ANTIBOND_OCC_THRESHOLD, not just the
    single lowest (see that constant's comment for why). Returns the AO-basis candidate
    vectors plus diagnostics (occupations of every local eigenvector, so a caller can
    see how cleanly separated the "antibond" cluster is from the rest of the local
    spectrum)."""
    a0, b0 = atom_a - 1, atom_b - 1
    local_idx = np.where((atom_of_iao == a0) | (atom_of_iao == b0))[0]
    if len(local_idx) < 2:
        raise ValueError(f"atoms {atom_a}/{atom_b}: fewer than 2 IAOs found -- can't form a pair block")

    dm_ao = pair_density_matrix(mf)
    local_iaos = iaos_orth[:, local_idx]
    dm_local = local_iaos.T @ s @ dm_ao @ s @ local_iaos  # (n_local, n_local), IAO basis is S-orthonormal

    occupations, evecs = np.linalg.eigh(dm_local)  # ascending occupation order
    n_candidates = max(1, int(np.searchsorted(occupations, ANTIBOND_OCC_THRESHOLD)))
    candidates_ao = [local_iaos @ evecs[:, k] for k in range(n_candidates)]

    return {
        "candidates_ao": candidates_ao,
        "candidate_occupations": occupations[:n_candidates],
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


if __name__ == "__main__":
    import sys
    name = sys.argv[1] if len(sys.argv) > 1 else "mol_002"
    result = run_case(name)
    print(f"\n=== {result['case']} ({result['basis_note']}) -- local per-atom-pair antibond ===")
    for label, key in [("wCNmax", "cn"), ("w17max", "w17"), ("w78max", "w78")]:
        r = result[key]
        p = r["pair"]
        print(
            f"  {label}: wmax={r['wmax']:.4f}  MO={r['mo_index']}  "
            f"[local block: {p['n_local_iaos']} IAOs, {len(p['candidates_ao'])} antibond candidate(s), "
            f"winning_occ={r['antibond_occupation']:.4f}, bond_occ={p['bond_occupation']:.4f}]"
        )
