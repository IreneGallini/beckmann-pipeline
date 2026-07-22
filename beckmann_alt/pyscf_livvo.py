"""
Shared PySCF SCF/projection helpers for the open-source wCNmax prototype
(`beckmann_alt/pair_nbo.py`) -- see Notes_open_source_alt.md for the full writeup,
caveats, and comparison against the trusted NBO7 numbers.

Pipeline, mirroring the NBO7 CMO step conceptually (not numerically):
  1. Single-point RKS wB97X-D/6-311+G(d,p) + ddCOSMO(water) on the frozen DFT geometry
     Gaussian already converged to (see beckmann_alt.geometry).
  2. Build a local per-atom-pair antibond candidate (beckmann_alt.pair_nbo).
  3. Project every canonical virtual MO in the full virtual manifold onto that
     candidate, square the overlap, take the max -- the wX^max analog
     (`project_virtuals_onto_livvo` below).

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
from pyscf import dft, gto, solvent
from pyscf.scf import dispersion

WINDOW_AU = 0.4  # LUMO .. LUMO + WINDOW_AU, same window convention as beckmann/dft/parse_cmo.py
                 # -- kept only as a diagnostic on the returned result, see
                 # project_virtuals_onto_livvo's docstring for why the search itself
                 # isn't capped to this window.


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


def project_virtuals_onto_livvo(mf, s, target_livvo, window_au: float = WINDOW_AU) -> dict:
    """Analog of beckmann.dft.parse_cmo.max_weight_for_target: for every virtual MO,
    project onto target_livvo (both S-normalized MO-like vectors in the same AO basis,
    so the overlap <MO_i|S|target> is the direct analog of an NBO7 CMO expansion
    coefficient), square it, take the max.

    Searches the FULL virtual manifold, not just LUMO..LUMO+window_au, even though
    window_au is still accepted and recorded as delta_lumo/in_window per result. The
    main pipeline's own wX^max (beckmann/dft/parse_cmo.py) used to cap this search at
    the window and that was identified as a bug there -- for some substrates the real
    antibond character peaks beyond LUMO+0.4 a.u. Capping the search here would
    silently reproduce a bug the main pipeline already fixed. window_au is kept only
    as a diagnostic in the returned dict.
    """
    mol = mf.mol
    nocc = mol.nelectron // 2
    mo_energy = mf.mo_energy
    mo_coeff = mf.mo_coeff
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
