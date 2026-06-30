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

Step 2 — substructure match

match = mol.GetSubstructMatch(OXIME_PAT)

GetSubstructMatch returns a tuple of RDKit atom indices (0-based) for the atoms that match the query in order of the map numbers :1, :2, :3. For mol_019_E it returns (2, 1, 0):
- atom index 2 → C (map label 1)
- atom index 1 → N (map label 2)
- atom index 0 → O (map label 3)

Step 3 — convert to Gaussian 1-based numbering

ci, ni, oi = (idx + 1 for idx in match)
oxime_label = f"[oxime: C{ci}=N{ni}-O{oi}]"

The coords list is built by iterating enumerate(mol.GetAtoms()), so atom at RDKit index i lands on line i+1 of the coordinate block. Adding 1 converts to Gaussian's 1-based numbering. For mol_019_E: C3=N2-O1.

Why this matters for NBO7

When the Gaussian job finishes, the NBO7 output will contain entries like:

     10. BD ( 1) C  3 - N  2        → this is the C=N π bond
     ...
     E2PERT:  BD*(1) N  2 - O  1 /  BD  C  3 - N  2   15.3 kcal/mol

Because the .gjf title already says [oxime: C3=N2-O1], you can write a parser that reads the label, extracts the three atom numbers, and uses them as keys to pull the right NBOs out of the output — without hard-coding any indices. The same parsing logic will work for all 68 molecules because each file carries its own map.

The connection to the tests

test_step3_oxime_atom_map_matches_sdf closes the loop: it reads each .gjf back, re-runs the SMARTS match on the SDF molecule, and asserts that the label in the file equals the match result. This catches any future bug where, say, the coordinate ordering and the SMARTS match diverge.

# Testing HPC steps
Set of 5 and 6 membered substrates with methyl or methoxy substituents in position 4 
- mols: 2, 6, 20, 21

## Implement Orbital Resolved Electron Routing Framework

---

## NBO Output: mol_002_E (first completed DFT run)

**Experimental outcome: F (fragmentation), 100% product B**
**DFT level: wB97XD/6-311+G(d,p), NBO7 single-point on optimised geometry**

### Atom map for mol_002_E

The `.gjf` title line carries: `[oxime: C11=N12-O13]`

This means in the coordinate block (and in all NBO output line references):
- Atom 11 = C (the oxime carbon, C=N)
- Atom 12 = N (the imine nitrogen)
- Atom 13 = O (the protonated leaving group, OH2+)

The two C–C bonds flanking the oxime carbon C11:
- **C6–C11**: aryl bond — connects the aromatic ring (C6 is the ipso-like ring carbon) to the oxime carbon
- **C10–C11**: alkyl bond — the methylene carbon on the other side of the ring

These two bonds are the candidates for migration. In classical Beckmann, the bond anti to the leaving group migrates. In the CN-handoff picture, the bond that donates more strongly into the N–O σ* (and reorganises the virtual manifold) is the one that migrates — or fragments.

### E2PERT key interactions

**Donors into BD\*(1) N12–O13** (the N–O σ\* = the bond being broken by the leaving group):

| Donor | Bond type | E2 (kcal/mol) | E(j)–E(i) (a.u.) | F(i,j) (a.u.) |
|---|---|---|---|---|
| BD(1) C6–C11 | aryl C–C σ | **12.63** | 0.83 | 0.091 |
| BD(1) C10–C11 | alkyl C–C σ | **3.38** | 0.80 | 0.047 |
| CR(1) N12 | N core orbital | 1.60 | 14.39 | 0.138 |

The aryl bond donates ~3.7× more strongly into the breaking N–O bond. A naive E2 analysis would predict aryl migration → rearrangement. The experiment gives fragmentation.

**Donors into BD\*(1) and BD\*(2) C11=N12** (the C=N σ\* and π\* — relevant for CN-handoff):

| Donor | Bond type | Acceptor | E2 (kcal/mol) |
|---|---|---|---|
| BD(2) C6–C7 | aryl π bond | BD\*(2) C11–N12 π\* | **55.20** |
| BD(1) C10–H22/H23 | alkyl C–H σ | BD\*(2) C11–N12 π\* | 7.42 / 7.41 |
| BD(1) C9–C10 | alkyl C–C σ | BD\*(1) C11–N12 σ\* | 5.63 |
| BD(1) C6–C7 | aryl C–C σ | BD\*(1) C11–N12 σ\* | 3.49 |
| BD(1) C6–C11 | aryl C–C σ | BD\*(1) C11–N12 σ\* | 2.80 |

The 55.20 kcal/mol aryl π → C=N π\* interaction is the ground-state aromatic conjugation into the imine. This is expected (resonance), but its magnitude sets the baseline for how much the C=N π\* is stabilised by the aryl side before N–O activation begins.

### Interpretation and the failure of simple E2 analysis

The simple E2 rank (aryl > alkyl into σ\*(N–O)) predicts aryl migration → R.
The experiment is unambiguously F (100%).

This is the core case the CN-handoff model needs to explain. A working hypothesis:

> The aryl π system conjugates so strongly into C=N π\* (55.20 kcal/mol) that as the N–O bond lengthens, the lowest unoccupied orbital reorganises away from σ\*(N–O) / σ\*(C–aryl) toward a CN-like character. The migrating group never builds up the required orbital overlap to complete rearrangement, and fragmentation wins instead.

The parse_nbo.py parser should test whether the ratio of aryl π → C=N π\* vs aryl σ → N–O σ\* (55.20 vs 12.63 here) is a descriptor that separates F from R cases in the benchmark.

### What parse_nbo.py must extract (minimum viable descriptor set)

For each molecule, using the `[oxime: C{ci}=N{ni}-O{oi}]` label from the `.gjf` title to identify atom numbers:

1. `E2_aryl_to_NO_star`: BD(1) C_aryl–C{ci} → BD\*(1) N{ni}–O{oi}
2. `E2_alkyl_to_NO_star`: BD(1) C_alkyl–C{ci} → BD\*(1) N{ni}–O{oi}
3. `E2_aryl_pi_to_CN_pi_star`: BD(2) aryl → BD\*(2) C{ci}–N{ni} (the 55 kcal/mol term)
4. `E2_aryl_to_CN_star`: BD(1) C_aryl–C{ci} → BD\*(1) C{ci}–N{ni}
5. `E2_alkyl_to_CN_star`: BD(1) C_alkyl–C{ci} → BD\*(1) C{ci}–N{ni}
6. Wiberg bond indices for N{ni}–O{oi}, C{ci}–N{ni}, C_aryl–C{ci}, C_alkyl–C{ci} (from BNDIDX)

**Key challenge:** identifying which neighbour of C{ci} is aryl and which is alkyl without hard-coding atom numbers. Two approaches:
- Use RDKit on the SDF to label the two C{ci} neighbours before running NBO — write aryl atom index into the `.gjf` title alongside the oxime label (e.g. `[oxime: C11=N12-O13 | aryl=C6 alkyl=C10]`)
- In the parser, identify them from BNDIDX: the aryl neighbour will have a Wiberg C–C index > 1.3 (aromatic), the alkyl will be near 1.0

The extended label approach is cleaner because it keeps all atom assignments in one place and does not require bond order logic in the parser.

### Atom map consistency across the benchmark

Each molecule in the benchmark has a different atom map because RDKit writes atoms in the order they appear in the SMILES, which varies by molecule. The `[oxime: C{ci}=N{ni}-O{oi}]` label in the `.gjf` title is the anchor that makes cross-molecule comparison possible:

- The parser reads the label, extracts `ci`, `ni`, `oi`
- It then scans the E2PERT table for lines where the acceptor column contains `BD*(1) N{ni} – O{oi}` or `BD*(2) C{ci} – N{ni}`
- This logic is molecule-agnostic; it does not need any hard-coded atom numbers

The test `test_sp_oxime_label_matches_sdf` already verifies the label is correct for all 34 molecules. When parse_nbo.py is added, a corresponding test should verify that the parser extracts a non-null E2 value for at least the two C–C → N–O σ\* entries in every completed log.

