"""
Ketone -> protonated activated oxime conversion.

Two-step conversion:
  1. C=O  ->  C=N-OH   (neutral oxime, via reaction SMARTS)
  2. N-OH  ->  N-[OH2+] (set formal charge +1 on O; sanitization gives 2 implicit H)
"""
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdchem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers, StereoEnumerationOptions,
)

KETONE_PAT  = Chem.MolFromSmarts('[C;$([C](=O)([c,C])[c,C])]')
NEUTRAL_RXN = rdChemReactions.ReactionFromSmarts(
    '[#6:2]-[C:1](=[O])-[#6:3]>>[#6:2]-[C:1](=NO)-[#6:3]'
)
OXIME_O_PAT = Chem.MolFromSmarts('[C]=[N]-[OH1]')


def get_cn_stereo(mol: Chem.Mol) -> str:
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for bond in mol.GetBonds():
        if bond.GetBondType() != rdchem.BondType.DOUBLE:
            continue
        if {bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()} == {'C', 'N'}:
            stereo = bond.GetStereo()
            if stereo == rdchem.BondStereo.STEREOE:
                return 'E'
            if stereo == rdchem.BondStereo.STEREOZ:
                return 'Z'
    return 'noEZ'


def ketone_to_protonated_oximes(mol: Chem.Mol) -> list[Chem.Mol]:
    """Two-step: C=O -> C=N-OH -> C=N-[OH2+]."""
    seen: set[str] = set()
    result: list[Chem.Mol] = []

    for prod_tuple in NEUTRAL_RXN.RunReactants((mol,)):
        p = prod_tuple[0]
        try:
            Chem.SanitizeMol(p)
        except Exception:
            continue

        match = p.GetSubstructMatch(OXIME_O_PAT)
        if not match:
            continue

        o_idx = match[2]
        rw = Chem.RWMol(p)
        rw.GetAtomWithIdx(o_idx).SetFormalCharge(1)
        try:
            Chem.SanitizeMol(rw)
            protonated = rw.GetMol()
            smi = Chem.MolToSmiles(protonated)
            if smi not in seen:
                seen.add(smi)
                result.append(protonated)
        except Exception:
            pass

    return result


def enumerate_ez(mol: Chem.Mol) -> list[tuple[Chem.Mol, str]]:
    opts = StereoEnumerationOptions(unique=True, onlyUnassigned=True)
    isomers = list(EnumerateStereoisomers(mol, options=opts))
    if not isomers:
        isomers = [mol]
    return [(iso, get_cn_stereo(iso)) for iso in isomers]
