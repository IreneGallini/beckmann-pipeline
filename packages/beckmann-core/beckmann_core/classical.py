"""
Classical anti-periplanar Beckmann rule: the substituent ANTI
(|dihedral| >= 150 deg) to the N-O leaving group migrates. Aryl anti ->
rearrangement (R); alkyl anti -> fragmentation (F). A quick geometric
baseline both products can show alongside their real (NBO7- or
PySCF-derived) prediction -- not a substitute for either.

get_oxime_atoms() is also the shared atom-mapping primitive used to resolve
ci/ni/oi/c_aryl/c_alkyl for a molecule from its RDKit Mol alone, no file I/O.
"""
from rdkit import Chem

# Uses neutral [O;!R] (not [O+]) since callers see a mix of protonation
# states depending on which stage of a product's pipeline built this Mol;
# the protonation state itself doesn't matter for dihedral geometry.
OXIME_PAT   = Chem.MolFromSmarts('[C:1]=[N:2]-[O;!R:3]')
ANTI_THRESH = 150.0


def get_oxime_atoms(mol: Chem.Mol):
    """(cox_idx, nox_idx, oox_idx, c_aryl_idx, c_allyl_idx), 0-based RDKit
    indices, or None if the oxime substructure or its aryl/alkyl neighbors
    aren't found."""
    match = mol.GetSubstructMatch(OXIME_PAT)
    if not match:
        return None
    cox_idx, nox_idx, oox_idx = match
    c_aryl_idx = c_allyl_idx = None
    for nbr in mol.GetAtomWithIdx(cox_idx).GetNeighbors():
        if nbr.GetIdx() == nox_idx:
            continue
        if nbr.GetIsAromatic():
            c_aryl_idx = nbr.GetIdx()
        else:
            c_allyl_idx = nbr.GetIdx()
    if c_aryl_idx is None or c_allyl_idx is None:
        return None
    return cox_idx, nox_idx, oox_idx, c_aryl_idx, c_allyl_idx


def predict(d_aryl: float, d_allyl: float) -> str:
    aryl_anti  = d_aryl  >= ANTI_THRESH
    allyl_anti = d_allyl >= ANTI_THRESH
    if aryl_anti and not allyl_anti:
        return 'R'
    if allyl_anti and not aryl_anti:
        return 'F'
    if aryl_anti and allyl_anti:
        return 'R' if d_aryl >= d_allyl else 'F'
    return 'inspect'
