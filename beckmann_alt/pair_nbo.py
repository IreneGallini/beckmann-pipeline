"""
Local per-atom-pair construction for an open-source wCNmax analog: instead of drawing
from a pre-built global set of localized virtual orbitals, build a fresh, local
subspace for each atom pair directly from that pair's own block of the density matrix
(in a Löwdin-orthogonalized IAO basis) -- structurally the same operation as NBO's own
per-atom-pair antibond search, just applied to one requested pair at a time.

Scoped to wCNmax only -- the actual predictive descriptor this project needs. Earlier
exploration (LIVVO-based and crude-AO-projection channel identification, w17max/w78max,
an iterative-deflation follow-up) is preserved in `Notes_open_source_alt.md` and in this
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
     MO onto it (beckmann_alt.pyscf_livvo.project_virtuals_onto_livvo -- same projection
     math regardless of where the target vector came from), and keep whichever candidate
     gives the largest projected weight -- mirroring the real wX^max definition's own
     "both BD*(1) and BD*(2) are eligible, whichever gives the larger squared
     coefficient wins" rule (Notes.md) exactly, just applied to locally-built candidates
     instead of NBO's own pre-labeled BD*(1)/BD*(2) pair.

Result (Notes_open_source_alt.md, "Second follow-up"): wCNmax within 1-6% of the
trusted NBO7 value on both reference cases (mol_002_E, 5_s0_Me) -- the best result of
any method tried in this exploration, and the only one that reproduces wCNmax's central
diagnostic signature (clearly the dominant channel) in both cases.
"""
import numpy as np
from pyscf import gto, lo

from beckmann_alt.geometry import TEST_IDS, load_case, load_test_set_case, load_test_set_scan_series
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
    }


def run_from_case(case: dict) -> dict:
    """Shared by run_case()/run_test_set_case() below -- everything from a loaded case
    dict (either shape, see beckmann_alt.geometry.load_case/load_test_set_case) through
    to a wCNmax result."""
    mol = build_mol(case)
    mf = run_scf(mol)
    iaos_orth, s, atom_of_iao = build_local_iaos(mf)
    cn = compute_wcnmax(mf, s, iaos_orth, atom_of_iao, case["ci"], case["ni"])
    return {"case": case["name"], "basis_note": case["basis_note"], "cn": cn}


def run_case(name: str) -> dict:
    """One of the two hand-picked reference cases (mol_002, 5_s0_Me) -- see
    beckmann_alt.geometry.REFERENCE_CASES."""
    return run_from_case(load_case(name))


def run_test_set_case(mol_id: str) -> dict:
    """Any of the six main-pipeline test-set molecules (mol_002/006/014/020/021/029) --
    see beckmann_alt.geometry.load_test_set_case."""
    return run_from_case(load_test_set_case(mol_id))


def run_test_set_scan_series(mol_id: str, stages: list[str] | None = None) -> list[dict]:
    """wCNmax at every R(N-O) point of a test-set molecule's scan series
    (beckmann_alt.geometry.load_test_set_scan_series) -- one PySCF SCF +
    local per-atom-pair projection per point, not just the single
    equilibrium geometry run_test_set_case() computes. Returns one row per
    point shaped like a beckmann.dft.parse_cmo 'cn'-channel extraction row
    (mol/stage/channel/R_NO/MO_index/weight/...) so
    beckmann.dft.descriptors.find_wcnmax_minimum() can be called on the
    result directly, reusing the main pipeline's own interior-minimum
    criterion rather than reimplementing it here.

    stages, if given, restricts the run to just those stage labels (e.g.
    mol_034_E's STEP_SCAN_SOURCES-merged series has 12 points, ~2x every
    other molecule's 6 -- pass a 6-point subset to keep runtime comparable).
    """
    cases = load_test_set_scan_series(mol_id)
    if stages is not None:
        cases = [c for c in cases if c["stage"] in stages]
    rows = []
    for case in cases:
        result = run_from_case(case)
        cn = result["cn"]
        rows.append({
            "mol": case["name"], "stage": case["stage"], "channel": "cn",
            "R_NO": case["r_no"], "MO_index": cn["mo_index"],
            "epsilon_i_star": cn["epsilon"], "coefficient": cn["coefficient"],
            "weight": cn["wmax"], "delta_lumo": None, "in_window": None,
        })
    return rows


def _print_wcnmax(result: dict) -> None:
    r, p = result["cn"], result["cn"]["pair"]
    print(f"\n=== {result['case']} ({result['basis_note']}) -- wCNmax, local per-atom-pair antibond ===")
    print(
        f"  wCNmax: wmax={r['wmax']:.4f}  MO={r['mo_index']}  "
        f"[local block: {p['n_local_iaos']} IAOs, {len(p['candidates_ao'])} antibond candidate(s), "
        f"winning_occ={r['antibond_occupation']:.4f}, bond_occ={p['bond_occupation']:.4f}]"
    )


if __name__ == "__main__":
    import sys
    name = sys.argv[1] if len(sys.argv) > 1 else "mol_002"
    result = run_test_set_case(name) if name in TEST_IDS else run_case(name)
    _print_wcnmax(result)
