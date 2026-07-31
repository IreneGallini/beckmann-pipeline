"""
Local per-atom-pair construction for an PySCF wCNmax analog: instead of drawing
from a pre-built global set of localized virtual orbitals, build a fresh, local
subspace for each atom pair directly from that pair's own block of the density matrix
(in a Löwdin-orthogonalized IAO basis) -- structurally the same operation as NBO's own
per-atom-pair antibond search, just applied to one requested pair at a time.

Scoped to wCNmax only -- the actual predictive descriptor this project needs. Earlier
exploration (LIVVO-based and crude-AO-projection channel identification, w17max/w78max,
an iterative-deflation follow-up) is preserved in `Notes_pyscf_alt.md` and in this
branch's git history (`git log -- beckmann_alt/pair_nbo.py`), not repeated here.

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
  3. Slice out just the oxime C/N atoms' IAOs' rows/columns from D_IAO -- a small,
     LOCAL block (10-dimensional for two heavy atoms in 6-311+G(d,p)'s minao-sized
     IAO space).
  4. Diagonalize that local block directly (a legitimate, standard NBO-style
     approximation -- NBO itself diagonalizes each 2-center block as if it were an
     isolated 2-atom system). Every eigenvector below ANTIBOND_OCC_THRESHOLD is kept as
     an antibond CANDIDATE, not just the single lowest-occupation one -- see that
     constant's comment: C=N is a double bond, so the local block has TWO
     low-occupation eigenvectors (sigma*-like and pi*-like), and it is NOT always the
     lowest-occupation one that actually carries the physically dominant character.
  5. Rotate each candidate back into the full AO basis, project every canonical virtual
     MO onto it (beckmann_pyscf.engine.pyscf_livvo.project_virtuals_onto_livvo -- same projection
     math regardless of where the target vector came from), and keep whichever candidate
     gives the largest projected weight -- mirroring the real wX^max definition's own
     "both BD*(1) and BD*(2) are eligible, whichever gives the larger squared
     coefficient wins" rule (Notes.md) exactly, just applied to locally-built candidates
     instead of NBO's own pre-labeled BD*(1)/BD*(2) pair.

Result (Notes_pyscf_alt.md, "Second follow-up"): wCNmax within 1-6% of the
trusted NBO7 value on both reference cases (mol_002_E, 5_s0_Me) -- the best result of
any method tried in this exploration, and the only one that reproduces wCNmax's central
diagnostic signature (clearly the dominant channel) in both cases.
"""
import numpy as np
from pyscf import gto, lo

from beckmann_pyscf.engine.pyscf_livvo import build_mol, run_scf, project_virtuals_onto_livvo


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
# the C=N oxime bond, the local block has TWO low-occupation eigenvectors (sigma* and
# pi*-like), not one -- mol_002's cn=(C11,N12) block gives occupations [0.016, 0.219,
# 0.568, 1.02, 1.16, 1.98, ...], and it's the SECOND-lowest (0.219), not the lowest
# (0.016), that actually projects onto the virtual manifold close to NBO7's trusted
# wCNmax (0.460 vs trusted 0.436, MO47 vs trusted MO48 -- taking only the single lowest
# eigenvector gave 0.087, roughly 5x too small). This mirrors the real wX^max definition
# directly: "both the BD*(1) (sigma) and BD*(2) (pi) components ... whichever gives the
# larger squared coefficient wins" (Notes.md) -- we don't know in advance which local
# eigenvector plays that BD*(1)/BD*(2) role, so every sub-threshold candidate is
# projected and the best one is kept, same as NBO7's own convention.


def local_pair_antibonds(
    mf, s: np.ndarray, iaos_orth: np.ndarray, atom_of_iao: np.ndarray, atom_a: int, atom_b: int,
) -> dict:
    """Build every antibond-like CANDIDATE local combination for one atom pair (1-based
    atom_a/atom_b) -- every local eigenvector below ANTIBOND_OCC_THRESHOLD, not just the
    single lowest (see that constant's comment for why). Returns the AO-basis candidate
    vectors plus diagnostics (occupations of every local eigenvector, so a caller can
    see how cleanly separated the "antibond" cluster is from the rest of the local
    spectrum, and bond_occupation as a sanity check that this pair forms a genuine
    localized bond -- should read ~2.0)."""
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


def compute_wcnmax(mf, s: np.ndarray, iaos_orth: np.ndarray, atom_of_iao: np.ndarray, ci: int, ni: int) -> dict:
    """Note: 'second' (the runner-up MO for the winning antibond candidate, from
    project_virtuals_onto_livvo) is threaded through into the returned dict
    unchanged -- the candidate-selection loop below (sigma*-like vs. pi*-like local
    antibond) is a separate choice from MO-tracking and is untouched."""
    pair = local_pair_antibonds(mf, s, iaos_orth, atom_of_iao, ci, ni)
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
        "second": best["second"],
    }


def compute_aryl_coeff_at_mos(
    mf, s: np.ndarray, iaos_orth: np.ndarray, atom_of_iao: np.ndarray,
    ci: int, c_aryl: int, mo_indices: list[int],
) -> dict[int, float]:
    """The aryl-migrating C-C antibond's (local per-atom-pair candidate for (ci,
    c_aryl), built with the exact same local_pair_antibonds() machinery
    compute_wcnmax() uses for (ci, ni) -- it's already generic over any atom
    pair) own signed coefficient projected onto specific, already-known MOs --
    e.g. the CN channel's winning and runner-up MO indices -- rather than
    searched for its own best match. Picks whichever local candidate
    (sigma*-/pi*-like) gives the largest projected weight across the whole
    virtual manifold once (mirroring compute_wcnmax's own candidate-selection
    convention), then reports that one candidate's coefficient at each
    requested MO index.

    Returns {mo_index: coefficient} -- entries for None in mo_indices are
    skipped (callers pass None when there's no runner-up MO, see
    compute_wcnmax()'s 'second')."""
    pair = local_pair_antibonds(mf, s, iaos_orth, atom_of_iao, ci, c_aryl)
    best_weight = None
    best_candidate = None
    for vec in pair["candidates_ao"]:
        proj = project_virtuals_onto_livvo(mf, s, vec)
        if best_weight is None or proj["weight"] > best_weight:
            best_weight, best_candidate = proj["weight"], vec

    mo_coeff = mf.mo_coeff
    return {
        mo_index: mo_coeff[:, mo_index] @ s @ best_candidate
        for mo_index in mo_indices
        if mo_index is not None
    }


def run_from_case(case: dict) -> dict:
    """Everything from a loaded case dict (atom_spec/charge/spin/ci/ni/c_aryl/
    name/basis_note -- see beckmann-pyscf's own backend/pipeline/wcnmax_pyscf.py
    for how a case is built from an AIMNet2 geometry, or
    research/pyscf_validation/geometry.py for the benchmark/reference-log
    loaders used in validation) through to a wCNmax result. Nothing in this
    function's body references any fixed benchmark set -- any case dict with
    these keys works.

    Also computes the aryl-migrating C-C antibond's own coefficient at the CN
    channel's winning MO (and runner-up MO, if one exists) via
    compute_aryl_coeff_at_mos(). Returned under "aryl_coeffs" as
    {mo_index: coefficient}."""
    mol = build_mol(case)
    mf = run_scf(mol)
    iaos_orth, s, atom_of_iao = build_local_iaos(mf)
    cn = compute_wcnmax(mf, s, iaos_orth, atom_of_iao, case["ci"], case["ni"])

    mo_indices = [cn["mo_index"]]
    if cn["second"] is not None:
        mo_indices.append(cn["second"]["mo_index"])
    aryl_coeffs = compute_aryl_coeff_at_mos(
        mf, s, iaos_orth, atom_of_iao, case["ci"], case["c_aryl"], mo_indices,
    )
    return {
        "case": case["name"], "basis_note": case["basis_note"], "cn": cn,
        "aryl_coeffs": aryl_coeffs,
    }
