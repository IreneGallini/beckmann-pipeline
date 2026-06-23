# Running script 0, 1, and 2 on chemdraw.txt
All 22 entries parse correctly with a valid oxime group. The ester (COC(=O)) in molecules 1004-271 and 274-830 was left untouched — only the ring ketone C=O was converted.

One thing worth noting about 234-667: the original SMILES [O]C(C=C12)=C(F)C=C2CCC1=O parsed fine in RDKit, which canonicalized the [O] into [O]c1cc2c(cc1F)CC/C2=N\O (a phenol-like oxygen on an aromatic ring). If Auto3D has trouble with the [O] notation, you can manually replace it with O (add an implicit H) in molecules.smi for that pair of entries. Everything else is ready to pass into script 01.

Summary of scripts/00_ketones_to_oximes.py:
- Parses chemdraw.txt (ID / SMILES / blank format)
- Detects ketone C=O via SMARTS (correctly ignores esters and phenols)
- Converts C=O → C=N-OH using a reaction SMARTS (ring connectivity preserved)
- Enumerates RDKit tautomers, keeping only true oxime forms
- Generates both E and Z isomers by reading BondStereo on the C=N bond
- Deduplicates by canonical SMILES and writes data/input/molecules.smi

It's running well. Key status:
- All 22 SMILES validated
- 83 total 3D conformers generated (before ranking to top 5)
- Now in AIMNet2 optimization — at step 46/5000, ~1.1 s/step

- 58 conformers across all 22 oxime structures (up to 5 per molecule, AIMNet2-ranked)
- All 22 E/Z oxime pairs made it through — every structure is accounted for
- Rigid/symmetric structures (374-658, 234-667, 0544-891) got 1–2 conformers; flexible ones with bulky substituents (314-235, 0924-630) got the full 5


# Week 3 

Week 1 deliverable: Submit: (i) optimized E/Z structures for all substrates, (ii) one CSV summary table, and (iii) a short note identifying which substrates agree or disagree with experiment and which require manual inspection. 

- Draw substrates in chemdraw and add to existing chemdraw.txt --> figure out consistent id numbering 
- Protonate oxime
- extract min energy isomer -> is it E or Z?
- dihedral extraction of output (of OH group use it to predict rearrangement vs fragmentation) --> compare with experimental results.  O_oxime–N_oxime–C_oxime–C_aryl  Chem.rdMolTransforms.GetDihedralDeg()
- consistent atom mapping 

- Translate a structure into an atom-index map: oxime O, oxime N, oxime C, aryl/allyl carbon, and ring atoms. Important for orbital descriptors 
    - Make atom_map_template.csv for two molecules and explain each index
- run NBO7 descriptors
    - Parse one log into CSV and confirm exactly 5 scan rows
    - Plot or tabulate each descriptor versus R_NO
- Run the workflow in a controlled folder structure and never overwrite raw logs --> Use AIMNet/Auto3D for initial conformer filtering, then DFT optimization, then constrained N-O elongation with NBO7 at each relaxed point.
- Document basis, charge, multiplicity, atom map, completion status, and failure mode for every molecule

# Week 4
GOAL 1: run all best conformers through Gaussian 
- Message:
    - pick one representative molecule mol 19 or 22 clean 100% R, 6-membered ring, no heteroatoms 
    - ask which level of theory and basis set
    - show input 
- add testing for script 03 
    - charge should be +1
    - atom map .gjf should match oxime atom indices 


### Atom mapping
Gaussian numbers atoms by their position in the coordinate block atom 1 is the first line, atom 2 is the second, and so on. NBO7's output (E2PERT donor/acceptor table, bond indices, etc.) refers back to these same numbers. But which atom is the oxime carbon? Which is nitrogen? That depends entirely on how the SDF was written, and it varies molecule to molecule.

The label [oxime: C3=N2-O1] in the .gjf title is a human-readable bookmark: "in this particular file, the C=N–O atoms are at positions 3, 2, 1." When you later parse NBO output for the C=N π-bond or the N–O σ* orbital, you know exactly which atom numbers to look for without opening Avogadro.

---
Step 1 — SMARTS pattern

OXIME_PAT = Chem.MolFromSmarts('[C:1]=[N:2]-[O+:3]')

This is a SMARTS query with atom map numbers (:1, :2, :3). The SMARTS encodes the connectivity of the activated protonated oxime:
- [C:1] — any carbon, labelled 1
- =[N:2] — double bond to any nitrogen, labelled 2
- -[O+:3] — single bond to a positively charged oxygen (the [OH2+]), labelled 3

The [O+] is the key fix from earlier — the neutral [OH1] pattern never matched because our molecules are protonated activated oximes (C=N-[OH2+]), not neutral hydroxylamine oximes.

---
Step 2 — substructure match

match = mol.GetSubstructMatch(OXIME_PAT)

GetSubstructMatch returns a tuple of RDKit atom indices (0-based) for the atoms that match the query in order of the map numbers :1, :2, :3. For mol_019_E it returns (2, 1, 0):
- atom index 2 → C (map label 1)
- atom index 1 → N (map label 2)
- atom index 0 → O (map label 3)

---
Step 3 — convert to Gaussian 1-based numbering

ci, ni, oi = (idx + 1 for idx in match)
oxime_label = f"[oxime: C{ci}=N{ni}-O{oi}]"

The coords list is built by iterating enumerate(mol.GetAtoms()), so atom at RDKit index i lands on line i+1 of the coordinate block. Adding 1 converts to Gaussian's 1-based numbering. For mol_019_E: C3=N2-O1.

---
Why this matters for NBO7

When the Gaussian job finishes, the NBO7 output will contain entries like:

     10. BD ( 1) C  3 - N  2        → this is the C=N π bond
     ...
     E2PERT:  BD*(1) N  2 - O  1 /  BD  C  3 - N  2   15.3 kcal/mol

Because the .gjf title already says [oxime: C3=N2-O1], you can write a parser that reads the label, extracts the three atom numbers, and uses them as keys to pull the right NBOs out of the output — without hard-coding any indices. The same parsing logic will work for all 68 molecules because each file carries its own map.

---
The connection to the tests

test_step3_oxime_atom_map_matches_sdf closes the loop: it reads each .gjf back, re-runs the SMARTS match on the SDF molecule, and asserts that the label in the file equals the match result. This catches any future bug where, say, the coordinate ordering and the SMARTS match diverge.
    

