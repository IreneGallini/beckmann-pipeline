"""
PySCF + LIVVO (Localized Intrinsic Valence Virtual Orbitals) prototype for an
open-source wCNmax/w17max/w78max analog -- see Notes_open_source_alt.md for the full
writeup, caveats, and comparison against the trusted NBO7 numbers.

Pipeline, mirroring the NBO7 CMO step conceptually (not numerically):
  1. Single-point RKS wB97X-D/6-311+G(d,p) + SMD(water) on the frozen DFT geometry
     Gaussian already converged to (see beckmann_alt.geometry).
  2. Build IAOs on the occupied space, then LIVVOs (pyscf.lo.vvo) -- a small, valence-
     sized set of localized virtual orbitals near the frontier.
  3. Identify which LIVVO corresponds to each target channel (CN / aryl / alkyl) by
     Mulliken population weight on the two atoms defining that antibond.
  4. Project every canonical virtual MO in the LUMO..LUMO+0.4 a.u. window onto the
     identified LIVVO, square the overlap, take the max -- the wX^max analog.

Known, deliberate deviations from Gaussian (do not silently treat as equivalent):
  - Solvent: attempted pyscf.solvent.smd first -- its Python API and solvent parameter
    table (water entry present) are real, but its non-electrostatic CDS term requires a
    compiled extension (libsolvent, built with -DENABLE_SMD=ON) that this pip wheel
    does not ship (RuntimeError: "SMD module is not available"). Building pyscf from
    source to get it is a bigger ask than the plain `pip install` that was approved, so
    this prototype falls back to ddCOSMO instead. ddCOSMO is NOT equivalent to SMD --
    different cavity construction (fixed atomic-radii spheres vs. SMD's own tessellated
    solvent-accessible surface), no explicit non-electrostatic (cavitation/dispersion/
    solvent-structure) term at all. Treat the solvent term as a rough implicit-solvent
    approximation, not a reproduction of Gaussian's SMD/water calculation.
  - Functional: 'wb97x-d' is blacklisted by PySCF's high-level RKS dispersion dispatch
    (pyscf.scf.dispersion._black_list) -- that blacklist exists to stop PySCF silently
    stacking an extra D3/D4 correction on top of a functional that already has its own
    dispersion baked in (which is exactly what Gaussian's wB97XD is). We remove it from
    the blacklist locally (see _enable_wb97xd()) rather than substitute a different,
    D3(BJ)-corrected variant (wb97x-d3bj) that PySCF does support out of the box --
    that variant uses a different dispersion correction than Gaussian's own, which
    would be a bigger, less-transparent deviation than bypassing an overcautious
    string check for a functional libxc already implements correctly.
  - Density fitting: enabled for tractable runtime (a 358-basis-function RSH-hybrid
    SCF without it did not converge in a reasonable time on this machine). This is a
    standard, well-controlled approximation to the two-electron integrals, not
    expected to be a significant source of divergence -- but it is a deviation from
    Gaussian's conventional (non-DF) SCF and is flagged here for completeness.
"""
from pathlib import Path

import numpy as np
from pyscf import dft, gto, lo, solvent
from pyscf.scf import dispersion

from beckmann_alt.geometry import load_case

WINDOW_AU = 0.4  # LUMO .. LUMO + WINDOW_AU, same window convention as beckmann/dft/parse_cmo.py
CHANNEL_POP_THRESHOLD = 0.5  # minimum combined Mulliken population on the target atom
                             # pair for a LIVVO to be accepted as "the" channel orbital


def _enable_wb97xd() -> None:
    """See module docstring -- 'wb97x-d' is a real libxc functional, just blacklisted
    by PySCF's higher-level dispersion-dispatch safety check. Idempotent."""
    dispersion._black_list.discard("wb97x-d")


def build_mol(case: dict, basis: str = "6-311+g(d,p)") -> gto.Mole:
    return gto.M(
        atom=case["atom_spec"], basis=basis,
        charge=case["charge"], spin=case["spin"],
        unit="Angstrom", verbose=0,
    )


def run_scf(mol: gto.Mole, xc: str = "wb97x-d", solvent_name: str = "water"):
    """Single-point RKS with density fitting + ddCOSMO solvent (see module docstring
    for why ddCOSMO, not SMD, despite our own pipeline using SMD/water). Returns the
    converged mf."""
    _enable_wb97xd()
    mf = dft.RKS(mol).density_fit()
    mf.xc = xc
    mf = solvent.ddCOSMO(mf)
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("SCF did not converge")
    return mf


def build_livvos(mf) -> tuple[np.ndarray, np.ndarray]:
    """IAOs (occupied space) -> LIVVOs (pyscf.lo.vvo). Returns (livvo_coeffs, overlap S).

    livvo_coeffs has shape (nao, n_livvo); columns are S-orthonormal (built from an
    orthogonal rotation of the already S-orthonormal canonical virtual MOs, see
    pyscf.lo.vvo.vvo), so Mulliken population and MO-projection formulas below can
    treat each column as a normalized MO-like vector directly.
    """
    mol = mf.mol
    s = mol.intor_symmetric("int1e_ovlp")
    nocc = mol.nelectron // 2
    orbocc  = mf.mo_coeff[:, :nocc]
    orbvirt = mf.mo_coeff[:, nocc:]
    iaos = lo.iao.iao(mol, orbocc)
    livvos = lo.vvo.vvo(mol, orbocc, orbvirt, iaos=iaos, s=s)
    return livvos, s


def atom_populations(orb: np.ndarray, s: np.ndarray, mol: gto.Mole) -> np.ndarray:
    """Mulliken population of a single MO-like vector (AO basis, S-normalized) on each
    atom. pop[A] = sum over AOs mu on atom A of orb[mu] * (S @ orb)[mu]."""
    aoslices = mol.aoslice_by_atom()
    sc = s @ orb
    pop = np.array([
        (orb[p0:p1] * sc[p0:p1]).sum()
        for p0, p1 in ((sl[2], sl[3]) for sl in aoslices)
    ])
    return pop


def identify_channel_livvo(livvos: np.ndarray, s: np.ndarray, mol: gto.Mole, atom_a: int, atom_b: int) -> dict:
    """Which LIVVO has the most combined Mulliken population on atoms atom_a/atom_b
    (1-based, matching the .gjf/[oxime: ...] convention used throughout beckmann/dft/).

    Returns the winning LIVVO's index, its combined population score, and per-atom
    breakdown -- callers should treat combined_pop < CHANNEL_POP_THRESHOLD as "no LIVVO
    clearly represents this channel" rather than trusting a low-confidence argmax.
    """
    a, b = atom_a - 1, atom_b - 1
    best = {"index": None, "combined_pop": -1.0, "pop_a": None, "pop_b": None}
    for k in range(livvos.shape[1]):
        pop = atom_populations(livvos[:, k], s, mol)
        combined = pop[a] + pop[b]
        if combined > best["combined_pop"]:
            best = {"index": k, "combined_pop": combined, "pop_a": pop[a], "pop_b": pop[b]}
    # Not hard-rejected -- callers get the number either way, with "accepted" as an
    # explicit low-confidence flag to report alongside it. Returning None here for a
    # sub-threshold channel would silently prevent any ranking/ordering comparison for
    # that channel at all, which is worse than reporting a flagged, low-confidence
    # number (see compute_channel_descriptor and Notes_open_source_alt.md -- in
    # practice combined_pop landed around 0.41-0.57 for mol_002's three channels, i.e.
    # LIVVOs are less sharply two-center-localized than NBO's explicitly per-atom-pair-
    # searched antibonds; that's a real, reportable difference between the methods, not
    # a bug to threshold away).
    best["accepted"] = best["combined_pop"] >= CHANNEL_POP_THRESHOLD
    return best


def interatomic_overlap_population(orb: np.ndarray, s: np.ndarray, mol: gto.Mole, atom_a: int, atom_b: int) -> float:
    """Mulliken OFF-diagonal (interatomic) overlap population between atom_a/atom_b
    (1-based) for one MO-like vector: 2 * sum_{mu in A, nu in B} orb[mu]*S[mu,nu]*orb[nu].

    This is the standard bonding/antibonding diagnostic (large positive = bonding,
    large negative = antibonding character concentrated specifically between these two
    atoms) -- unlike atom_populations()/identify_channel_livvo()'s plain sum of
    on-atom (diagonal) populations, which cannot see the relative PHASE between the two
    atoms' contributions at all. Two atoms both carrying large density in-phase (a
    generic delocalized orbital) and two atoms carrying the same density out-of-phase
    (an actual A-B antibond) score identically under a diagonal-population sum but very
    differently here -- which is the suspected cause of identify_channel_livvo's channel
    collisions (see Notes_open_source_alt.md's follow-up investigation).
    """
    aoslices = mol.aoslice_by_atom()
    a, b = atom_a - 1, atom_b - 1
    pa, pb = aoslices[a][2:4], aoslices[b][2:4]
    return 2.0 * orb[pa[0]:pa[1]] @ s[pa[0]:pa[1], pb[0]:pb[1]] @ orb[pb[0]:pb[1]]


def identify_channel_livvo_v2(livvos: np.ndarray, s: np.ndarray, mol: gto.Mole, atom_a: int, atom_b: int) -> dict:
    """Same role as identify_channel_livvo, but scores each LIVVO by how NEGATIVE its
    interatomic_overlap_population(atom_a, atom_b) is, rather than by summed on-atom
    population -- see that function's docstring for why. Most-negative overlap
    population wins (the strongest A-B antibonding character); ties/near-zero values
    mean no LIVVO shows meaningful A-B antibonding character at all.
    """
    a0, b0 = atom_a - 1, atom_b - 1
    best = {"index": None, "overlap_pop": None}
    for k in range(livvos.shape[1]):
        op = interatomic_overlap_population(livvos[:, k], s, mol, atom_a, atom_b)
        if best["overlap_pop"] is None or op < best["overlap_pop"]:
            best = {"index": k, "overlap_pop": op}
    return best


def project_virtuals_onto_livvo(mf, s: np.ndarray, target_livvo: np.ndarray, window_au: float = WINDOW_AU) -> dict:
    """Analog of beckmann.dft.parse_cmo.max_weight_for_target: for every virtual MO,
    project onto target_livvo (both S-normalized MO-like vectors in the same AO basis,
    so the overlap <MO_i|S|target> is the direct analog of an NBO7 CMO expansion
    coefficient), square it, take the max.

    Searches the FULL virtual manifold, not just LUMO..LUMO+window_au, even though
    window_au is still accepted and recorded as delta_lumo/in_window per result. The
    main pipeline's own wX^max (beckmann/dft/parse_cmo.py) used to cap this search at
    the window and that was identified as a bug there -- for some substrates the real
    antibond character peaks beyond LUMO+0.4 a.u., and mol_002's own trusted w17max/
    w78max (cmo_channel_extraction.csv) are themselves from MOs at delta_lumo=0.50/0.44,
    both OUTSIDE this window. Capping the search here would silently reproduce a bug
    the main pipeline already fixed, and would make any comparison against those two
    channels meaningless (comparing against a value the windowed search could not
    possibly find). window_au is kept only as a diagnostic in the returned dict.
    """
    mol = mf.mol
    nocc = mol.nelectron // 2
    mo_energy = mf.mo_energy
    mo_coeff  = mf.mo_coeff
    lumo_e = mo_energy[nocc]

    best = {"weight": None, "mo_index": None, "epsilon": None, "coefficient": None}
    for i in range(nocc, len(mo_energy)):
        coeff = mo_coeff[:, i] @ s @ target_livvo
        weight = coeff ** 2
        if best["weight"] is None or weight > best["weight"]:
            best = {"weight": weight, "mo_index": i, "epsilon": mo_energy[i], "coefficient": coeff}
    if best["epsilon"] is not None:
        delta_lumo = best["epsilon"] - lumo_e
        best["delta_lumo"] = delta_lumo
        best["in_window"] = delta_lumo <= window_au + 1e-9
    else:
        best["delta_lumo"] = best["in_window"] = None
    return best


def compute_channel_descriptor(mf, livvos: np.ndarray, s: np.ndarray, atom_a: int, atom_b: int) -> dict:
    """One channel end-to-end: identify the LIVVO, then wX^max via projection onto it.

    Always computes and returns wmax, even when identify_channel_livvo's population
    threshold isn't met -- "accepted" (nested under "livvo") is the flag callers should
    check/report, not a gate that silently withholds the number (see
    identify_channel_livvo's docstring).
    """
    channel_livvo = identify_channel_livvo(livvos, s, mf.mol, atom_a, atom_b)
    target = livvos[:, channel_livvo["index"]]
    proj = project_virtuals_onto_livvo(mf, s, target)
    return {
        "livvo": channel_livvo,
        "wmax": proj["weight"], "mo_index": proj["mo_index"],
        "epsilon": proj["epsilon"], "coefficient": proj["coefficient"],
        "delta_lumo": proj["delta_lumo"], "in_window": proj["in_window"],
    }


def run_case(name: str) -> dict:
    case = load_case(name)
    mol = build_mol(case)
    mf = run_scf(mol)

    livvos, s = build_livvos(mf)

    cn    = compute_channel_descriptor(mf, livvos, s, case["ci"], case["ni"])
    c17   = compute_channel_descriptor(mf, livvos, s, case["ci"], case["c_aryl"])
    c78   = compute_channel_descriptor(mf, livvos, s, case["ci"], case["c_alkyl"])

    return {
        "case": name, "basis_note": case["basis_note"],
        "n_livvo": livvos.shape[1],
        "cn": cn, "w17": c17, "w78": c78,
    }


if __name__ == "__main__":
    import sys
    name = sys.argv[1] if len(sys.argv) > 1 else "mol_002"
    result = run_case(name)
    print(f"\n=== {result['case']} ({result['basis_note']}) ===")
    print(f"n_livvo = {result['n_livvo']}")
    for label, key in [("wCNmax", "cn"), ("w17max", "w17"), ("w78max", "w78")]:
        r = result[key]
        acc = "accepted" if r["livvo"]["accepted"] else "REJECTED (below pop threshold)"
        print(
            f"  {label}: livvo#{r['livvo']['index']} combined_pop={r['livvo']['combined_pop']:.3f} [{acc}]  "
            f"wmax={r['wmax']}  MO={r['mo_index']}  eps={r['epsilon']}"
        )
