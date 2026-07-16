"""
Crude AO-projection fallback for wCNmax/w17max/w78max -- a same-afternoon sanity-check
baseline, NOT a second rigorous method. See Notes_open_source_alt.md.

Construction: for the two atoms defining a target antibond (e.g. C{ci}/N{ni} for the CN
channel), build each atom's "p-orbital pointing along the bond axis" by summing ALL of
that atom's p-type basis functions (every p-shell in the basis -- 6-311+G(d,p) has
several: 2p/3p/4p/5p for a first-row atom) dotted with the bond unit vector, equal
weight per shell. This mixes different radial extents with no attempt at a proper
hybridization/weighting scheme -- deliberately crude, see module-level caveat below.

Then form the standard two-center LCAO sigma/sigma* combination:
    trial_bond     = p_A(toward B) + p_B(toward A)   (in-phase / bonding-like)
    trial_antibond = p_A(toward B) - p_B(toward A)   (out-of-phase / antibonding-like)
S-normalize trial_antibond, and reuse pyscf_livvo.project_virtuals_onto_livvo (the
projection math is identical regardless of where the target vector came from) to get
the wX^max analog.

Caveat -- this is expected to be the weaker of the two prototypes: summing every p-shell
equally ignores that a real antibond is dominated by one particular radial extent (the
valence p-shell), not an unweighted mix of core-like and diffuse p functions. No
orthogonalization against the occupied space is performed either (unlike the IAO/VVO
route, which orthogonalizes by construction) -- the trial vector can have some overlap
with occupied orbitals that a proper virtual-space projection would remove first. Report
results as a rough same/different-ballpark check, not as comparably trustworthy to the
LIVVO prototype.
"""
import numpy as np
from pyscf import gto

from beckmann_alt.geometry import load_case
from beckmann_alt.pyscf_livvo import (
    WINDOW_AU, build_mol, project_virtuals_onto_livvo, run_scf,
)


def _atom_p_orbital_toward(mol: gto.Mole, atom_idx0: int, direction: np.ndarray) -> np.ndarray:
    """AO-basis vector for atom `atom_idx0` (0-based)'s p-character summed over every
    p-shell present in the basis, projected onto `direction` (a unit vector). Crude by
    construction -- see module docstring."""
    labels = mol.ao_labels(fmt=False)
    vec = np.zeros(mol.nao)
    for ao_i, (atom_id, _elem, shell, comp) in enumerate(labels):
        if atom_id != atom_idx0 or not shell.endswith("p") or comp not in ("x", "y", "z"):
            continue
        vec[ao_i] = direction["xyz".index(comp)]
    return vec


def build_trial_antibond(mol: gto.Mole, atom_a: int, atom_b: int) -> np.ndarray:
    """atom_a/atom_b are 1-based (matching the .gjf/[oxime: ...] convention)."""
    a, b = atom_a - 1, atom_b - 1
    coords = mol.atom_coords()  # Bohr, atomic units -- direction only, units cancel
    u = coords[b] - coords[a]
    u = u / np.linalg.norm(u)

    p_a_toward_b = _atom_p_orbital_toward(mol, a, u)
    p_b_toward_a = _atom_p_orbital_toward(mol, b, -u)

    trial_antibond = p_a_toward_b - p_b_toward_a
    return trial_antibond


def s_normalize(vec: np.ndarray, s: np.ndarray) -> np.ndarray:
    norm = np.sqrt(vec @ s @ vec)
    if norm < 1e-8:
        raise ValueError("trial vector has ~zero norm -- no p-type AOs found on these atoms?")
    return vec / norm


def compute_channel_descriptor(mf, s: np.ndarray, atom_a: int, atom_b: int) -> dict:
    trial = build_trial_antibond(mf.mol, atom_a, atom_b)
    trial = s_normalize(trial, s)
    proj = project_virtuals_onto_livvo(mf, s, trial, window_au=WINDOW_AU)
    return {
        "wmax": proj["weight"], "mo_index": proj["mo_index"],
        "epsilon": proj["epsilon"], "coefficient": proj["coefficient"],
    }


def run_case(name: str) -> dict:
    case = load_case(name)
    mol = build_mol(case)
    mf = run_scf(mol)
    s = mol.intor_symmetric("int1e_ovlp")

    cn  = compute_channel_descriptor(mf, s, case["ci"], case["ni"])
    c17 = compute_channel_descriptor(mf, s, case["ci"], case["c_aryl"])
    c78 = compute_channel_descriptor(mf, s, case["ci"], case["c_alkyl"])

    return {"case": name, "basis_note": case["basis_note"], "cn": cn, "w17": c17, "w78": c78}


if __name__ == "__main__":
    import sys
    name = sys.argv[1] if len(sys.argv) > 1 else "mol_002"
    result = run_case(name)
    print(f"\n=== {result['case']} ({result['basis_note']}) -- crude AO-projection ===")
    for label, key in [("wCNmax", "cn"), ("w17max", "w17"), ("w78max", "w78")]:
        r = result[key]
        print(f"  {label}: wmax={r['wmax']}  MO={r['mo_index']}  eps={r['epsilon']}")
