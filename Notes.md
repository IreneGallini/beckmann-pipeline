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