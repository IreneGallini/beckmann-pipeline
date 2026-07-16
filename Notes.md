### Atom mapping
Gaussian numbers atoms by their position in the coordinate block atom 1 is the first line, atom 2 is the second, and so on. NBO7's output (E2PERT donor/acceptor table, bond indices, etc.) refers back to these same numbers. But which atom is the oxime carbon? Which is nitrogen? That depends entirely on how the SDF was written, and it varies molecule to molecule.

The label [oxime: C3=N2-O1] in the .gjf title is a human-readable bookmark: "in this particular file, the C=N–O atoms are at positions 3, 2, 1." When you later parse NBO output for the C=N π-bond or the N–O σ* orbital, you know exactly which atom numbers to look for without opening Avogadro.

---
Step 1 SMARTS pattern

OXIME_PAT = Chem.MolFromSmarts('[C:1]=[N:2]-[O+:3]')

This is a SMARTS query with atom map numbers (:1, :2, :3). The SMARTS encodes the connectivity of the activated protonated oxime:
- [C:1] any carbon, labelled 1
- =[N:2] double bond to any nitrogen, labelled 2
- -[O+:3] single bond to a positively charged oxygen (the [OH2+]), labelled 3

The [O+] is the key fix from earlier the neutral [OH1] pattern never matched because our molecules are protonated activated oximes (C=N-[OH2+]), not neutral hydroxylamine oximes.

Step 2 substructure match

match = mol.GetSubstructMatch(OXIME_PAT)

GetSubstructMatch returns a tuple of RDKit atom indices (0-based) for the atoms that match the query in order of the map numbers :1, :2, :3. For mol_019_E it returns (2, 1, 0):
- atom index 2 → C (map label 1)
- atom index 1 → N (map label 2)
- atom index 0 → O (map label 3)

Step 3 convert to Gaussian 1-based numbering

ci, ni, oi = (idx + 1 for idx in match)
oxime_label = f"[oxime: C{ci}=N{ni}-O{oi}]"

The coords list is built by iterating enumerate(mol.GetAtoms()), so atom at RDKit index i lands on line i+1 of the coordinate block. Adding 1 converts to Gaussian's 1-based numbering. For mol_019_E: C3=N2-O1.

Why this matters for NBO7

When the Gaussian job finishes, the NBO7 output will contain entries like:

     10. BD ( 1) C  3 - N  2        → this is the C=N π bond
     ...
     E2PERT:  BD*(1) N  2 - O  1 /  BD  C  3 - N  2   15.3 kcal/mol

Because the .gjf title already says [oxime: C3=N2-O1], you can write a parser that reads the label, extracts the three atom numbers, and uses them as keys to pull the right NBOs out of the output without hard-coding any indices. The same parsing logic will work for all 68 molecules because each file carries its own map.

The connection to the tests

test_step3_oxime_atom_map_matches_sdf closes the loop: it reads each .gjf back, re-runs the SMARTS match on the SDF molecule, and asserts that the label in the file equals the match result. This catches any future bug where, say, the coordinate ordering and the SMARTS match diverge.

# Testing HPC steps
Set of 5 and 6 membered substrates with methyl or methoxy substituents in position 4 
- mols: 2, 6, 20, 21

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
- **C6–C11**: aryl bond connects the aromatic ring (C6 is the ipso-like ring carbon) to the oxime carbon
- **C10–C11**: alkyl bond the methylene carbon on the other side of the ring

These two bonds are the candidates for migration. In classical Beckmann, the bond anti to the leaving group migrates. In the CN-handoff picture, the bond that donates more strongly into the N–O σ* (and reorganises the virtual manifold) is the one that migrates or fragments.

### E2PERT key interactions

**Donors into BD\*(1) N12–O13** (the N–O σ\* = the bond being broken by the leaving group):

| Donor | Bond type | E2 (kcal/mol) | E(j)–E(i) (a.u.) | F(i,j) (a.u.) |
|---|---|---|---|---|
| BD(1) C6–C11 | aryl C–C σ | **12.63** | 0.83 | 0.091 |
| BD(1) C10–C11 | alkyl C–C σ | **3.38** | 0.80 | 0.047 |
| CR(1) N12 | N core orbital | 1.60 | 14.39 | 0.138 |

The aryl bond donates ~3.7× more strongly into the breaking N–O bond. A naive E2 analysis would predict aryl migration → rearrangement. The experiment gives fragmentation.

**Donors into BD\*(1) and BD\*(2) C11=N12** (the C=N σ\* and π\* relevant for CN-handoff):

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
- Use RDKit on the SDF to label the two C{ci} neighbours before running NBO write aryl atom index into the `.gjf` title alongside the oxime label (e.g. `[oxime: C11=N12-O13 | aryl=C6 alkyl=C10]`)
- In the parser, identify them from BNDIDX: the aryl neighbour will have a Wiberg C–C index > 1.3 (aromatic), the alkyl will be near 1.0

The extended label approach is cleaner because it keeps all atom assignments in one place and does not require bond order logic in the parser.

### Atom map consistency across the benchmark

Each molecule in the benchmark has a different atom map because RDKit writes atoms in the order they appear in the SMILES, which varies by molecule. The `[oxime: C{ci}=N{ni}-O{oi}]` label in the `.gjf` title is the anchor that makes cross-molecule comparison possible:

- The parser reads the label, extracts `ci`, `ni`, `oi`
- It then scans the E2PERT table for lines where the acceptor column contains `BD*(1) N{ni} – O{oi}` or `BD*(2) C{ci} – N{ni}`
- This logic is molecule-agnostic; it does not need any hard-coded atom numbers

The test `test_sp_oxime_label_matches_sdf` already verifies the label is correct for all 34 molecules. When parse_nbo.py is added, a corresponding test should verify that the parser extracts a non-null E2 value for at least the two C–C → N–O σ\* entries in every completed log.

## N–O Bond Scan: mol_002_E (Stage 3 relaxed PES scan)

**Job:** `mol_002_E_scan.gjf`  `wB97XD/6-311+G(d,p) opt=(ModRedundant,MaxCycles=200) pop=nboread geom=checkpoint`
**Scan:** N12–O13 bond stretched from R0 in 4 steps of 0.1 Å (5 points: R0 to R0+0.4 Å)
**NBO keywords:** `E2PERT BNDIDX NBOSUM`
**Status:** Normal termination (20 min wall time, 2h 47min CPU)

### How many NBO analyses ran

Gaussian ran NBO **twice**, not five times:
1. **At R0**: on the geometry read from `_opt.chk` before the scan optimisation starts (the equilibrium DFT geometry)
2. **At R0+0.4 Å**: on the final converged scan geometry (scan point 5)

Intermediate points (R0+0.1, +0.2, +0.3) did not produce separate NBO output. This is a Gaussian behaviour with `opt=(ModRedundant)` + `pop=nboread`: population analysis runs at the initial and final geometries only, not at each intermediate scan point. To get all 5 NBO analyses, we need 5 separate single-point NBO jobs at fixed N–O distances (or 3 additional jobs for the missing intermediate points).

### E2PERT at R0 (equilibrium, from Stage 2 _nbo.log confirmed in scan)

Atom map: C11=N12–O13 | aryl=C6, alkyl=C10

| Donor | Type | Acceptor | E2 (kcal/mol) |
|---|---|---|---|
| BD(1) C6–C11 | aryl C–C σ | BD\*(1) N12–O13 | **12.63** |
| BD(1) C10–C11 | alkyl C–C σ | BD\*(1) N12–O13 | **3.38** |

**Ψ(R0) = 12.63 / 3.38 = 3.74** (aryl dominates → classical prediction: R)

### Full 5-point E2PERT scan

The scan job only ran NBO at R0 and R0+0.4. The 3 intermediate geometries were extracted from the scan log using `scripts/dft/extract_scan_sp.py` and submitted as separate single-point jobs (`_sp2.gjf`, `_sp3.gjf`, `_sp4.gjf`).

| R(N–O) Å | Source | Aryl → acceptor | Alkyl → acceptor | Dominant acceptor |
|---|---|---|---|---|
| 1.6119 (R0) | scan initial + _nbo | **BD\*(1) N12–O13: 12.63** | BD\*(1) N12–O13: 3.38 | σ\*(N–O) |
| 1.7119 (R0+0.1) | sp2 | **BD\*(1) N12–O13: 15.70** | BD\*(1) N12–O13: 4.54 | σ\*(N–O) |
| 1.8119 (R0+0.2) | sp3 | **LP\*(2) N12: 21.80** | LP\*(2) N12: 7.23 | LP\*(N) + σ\*(C–N) |
| 1.9119 (R0+0.3) | sp4 | **LP\*(2) N12: 24.10** | LP\*(2) N12: 9.08 | LP\*(N) + σ\*(C–N) |
| 2.0119 (R0+0.4) | scan final | **LP\*(2) N12: 25.83** | LP\*(2) N12: 11.21 | LP\*(N) + σ\*(C–N) |

Ψ values (E2\_aryl / E2\_alkyl into the dominant N-O channel):

| R(N–O) Å | E2\_aryl (kcal/mol) | E2\_alkyl (kcal/mol) | Ψ |
|---|---|---|---|
| 1.6119 | 12.63 | 3.38 | **3.74** |
| 1.7119 | 15.70 | 4.54 | **3.46** |
| 1.8119 | 21.80 (into LP\*N) | 7.23 | **3.02** |
| 1.9119 | 24.10 | 9.08 | **2.66** |
| 2.0119 | 25.83 | 11.21 | **2.30** |

d/dR ≈ (2.30 − 3.74) / 0.4 = **−3.6 Å⁻¹** (Ψ decreasing as N–O stretches)

### CN-handoff: when and what changes

The σ\*(N–O) acceptor **disappears between R0+0.1 and R0+0.2 Å** (between 1.7119 and 1.8119 Å). Before this crossing:
- Both C–C bonds donate into BD\*(1) N12–O13 (σ\* of the breaking N–O bond)
- Aryl leads Ψ ~ 3.5–3.7×

After the crossing (R ≥ 1.8119 Å):
- σ\*(N–O) is no longer a distinct NBO
- The dominant acceptors are LP\*(2) N12 (N lone-pair antibonding) and BD\*(1/2) C11–N12
- Aryl still leads (Ψ ~ 2.3–3.0) but into a CN-like channel
- The aryl π system (BD(2) C6–C7 → BD\*(2) C11–N12) grows from 46 → 50 kcal/mol across the range

**For mol_002_E (experiment = F):** the equilibrium Ψ predicts aryl migration (R), but the CN-handoff occurs very early (between R0+0.1 and R0+0.2). As the N–O bond stretches, aryl is preferentially stabilising the developing C11–N12 σ\* channel (LP\*N and σ\*CN become the sinks), which may prevent the aryl group from accumulating the bond order needed to complete migration. The alkyl C10–C11 fragmentation pathway wins.

The decreasing Ψ slope (d/dR = −3.6 Å⁻¹) means aryl's relative advantage over alkyl SHRINKS as the reaction proceeds, consistent with the fragmentation outcome.

**Caveat on the Ψ number above:** this early hand calculation used a simplified proxy,
`E2_aryl / E2_alkyl` where both donors feed the *same* acceptor (the N-O σ*). The
preprint's actual Ψ definition (below, now implemented in code) is different: the
denominator `K_frag` is a **sum over every donor** feeding a **different** acceptor
(the developing alkyl/fragmentation-channel antibond `σ*(C_ox-C_alkyl)`), not a single
E2 value into the N-O channel. Don't reuse the 3.74/3.46/... numbers above as if they
were Ψ in the paper's sense — they aren't, they were an earlier ad-hoc approximation.

### CMO analysis: DONE (superseding the "still needed" note below)

CMO output comes directly from `pop=nbo7read` -- no separate `gennbo7`/`.47`-archive
post-processing step is needed (that was only ever necessary as a workaround while we
thought we were stuck with the bundled NBO 3.1; see the "NBO7 setup on Citadel"
section in CLAUDE.md for how `pop=nbo7read` was made to actually work). All 4 test
molecules (002, 006, 020, 021) have complete NBO7+CMO data at all 5 N-O scan points.

---

## Corrected descriptor formulas (per "Ring Size and Substituent Effects in the
## Beckmann Rearrangement", Sections 2.2-2.5) -- implemented in `beckmann/dft/parse_cmo.py`
## and `beckmann/dft/descriptors.py`

The original task description below (the "Implement Orbital Resolved Electron Routing
Framework" section) described Λ, Ψ, and wCNmax only in prose, before the preprint was
available for direct comparison. An earlier version of this codebase implemented Λ as
an **unrestricted max over the whole virtual window** -- that was wrong. This section
documents the corrected, paper-matched formulas actually implemented now.

**Aryl/alkyl role tagging (`beckmann/dft/descriptors.py::get_substituent_map()`):**
Every downstream channel-resolved descriptor needs to know which of the two carbons
bonded to the oxime carbon (`ci`) is the aryl-side substituent (`c_aryl`, migrates in
rearrangement) and which is the alkyl-side substituent (`c_alkyl`, leaves in
fragmentation). This is derived **fresh via RDKit** on every run -- loads
`best_per_substrate.sdf`, reuses `beckmann.analysis.classical.get_oxime_atoms()`
(checks `GetIsAromatic()` on each neighbor of the oxime carbon), and cross-validates
the result against the independently-derived `(ci, ni, oi)` parsed from the
molecule's own `.gjf` title line -- a mismatch raises rather than silently trusting
either source. Deliberately does **not** read the pre-computed `c_aryl_idx`/
`c_allyl_idx` columns already sitting in `classical_rule_results.csv` (even though
they'd give the identical answer) -- kept independent of that CSV/script by design.

**wX^max (generic form, covers `w17max`, `w78max`, and `wcnmax`):** for a target
antibond X, scan every virtual MO from the LUMO up to LUMO+0.4 a.u., take X's
coefficient in that MO's CMO expansion (0 if X doesn't appear in the printed >5%
list -- see caveat below), square it, take the max across the window. Both the
`BD*(1)` (sigma) and `BD*(2)` (pi) components of X are eligible; whichever gives the
larger squared coefficient wins for that MO.
- `w17max` = wX^max for X = `BD*(C{ci}-C{c_aryl})` (rearrangement channel)
- `w78max` = wX^max for X = `BD*(C{ci}-C{c_alkyl})` (fragmentation channel)
- `wcnmax` = wX^max for X = `BD*(C{ci}-N{ni})` (nitrilium/routing channel)

**Per-geometry channel extraction (`data/output/analysis/cmo_channel_extraction.csv`,
`beckmann/dft/parse_cmo.py::compute_channel_weights()`):** the wX^max summary columns
above are a max over the whole virtual window, which throws away *which* MO achieved
that max. Since canonical MOs are energy-ordered and can swap character between scan
points (e.g. the CN channel's leading MO could be MO 49 at one R and a different index
at the next), that identity shift matters for anything that wants to look at what's
actually happening around an extremum -- e.g. an avoided-crossing check (small
eigenvalue gap + character exchange between two nearby virtual MOs right around a
wCNmax extremum), which isn't implemented yet but needs this data preserved to be
built later. This CSV has one row per `(mol, stage, channel)` for `channel` in
`{cn, 17, 78}`, with columns `R_NO, MO_index, epsilon_i_star, coefficient, weight` --
`weight` is exactly the corresponding wX^max value, `coefficient` is the signed value
before squaring, `epsilon_i_star` is that MO's orbital energy (canonical MO
eigenvalue). This is additive to `cmo_descriptors.csv`, not a replacement -- the
existing summary columns are unchanged. `scripts/analysis/summarize_descriptors.py`'s
wCNmax-extremum check backfills `MO_index`/`epsilon_i_star` from this table when an
extremum is found, rather than only reporting a bare yes/no.

**Caveat on "0 if X doesn't appear":** Gaussian's CMO printout only lists contributions
above a 5% threshold ("Leading (> 5%) NBO Contributions to Molecular Orbitals"). A
missing entry means the coefficient's *square* is below 0.05 (|coefficient| < ~0.22),
not necessarily exactly zero. Treating "not printed" as "0" (per the formula's own
wording) is a defensible engineering approximation, but it does mean `w17max`/`w78max`
come out as `None` fairly often in our own 4 test molecules -- when that happens for
`w17max` specifically, `Lambda` hits its `1e-3` floor and swings to large values
(seen up to ~117 for mol_002). This is the formula behaving exactly as specified, not
a parsing bug, but it's worth flagging: it means Lambda's magnitude is quite sensitive
to whether the rearrangement-channel antibond happens to clear the 5% print cutoff at
a given geometry, not a smoothly-varying quantity.

**Lambda (`beckmann/dft/parse_cmo.py::compute_descriptors()`):**
```
Lambda     = max(w78max) / max(max(w17max), 1e-3)
log_lambda = log10(Lambda)
```
A fragmentation-channel-over-rearrangement-channel dominance ratio -- **not** an
unrestricted max (the earlier, wrong version of this code). Floor is on the
denominator only, per the paper.

**Psi (`beckmann/dft/descriptors.py::compute_psi_row()`):**
```
Psi = K_anti / (K_frag + epsilon)
```
- `K_anti` = E2PERT stabilization from donor `BD(C{c_aryl}-C{ci})` into acceptor
  `BD*(N{ni}-O{oi})` (summed if both sigma and pi donor components are present,
  though in practice this bond is a single sigma bond so usually only one row
  matches).
- `K_frag` = **sum** of every E2PERT row whose acceptor is `BD*(C{ci}-C{c_alkyl})`,
  regardless of donor. Verified by hand against `5_s0_Me.log`'s raw E2PERT table
  (16 separate donor rows summing to exactly `K_frag = 21.11` kcal/mol -- see Task 5
  validation below).
- `epsilon = 1e-6` -- **not specified anywhere in the paper or Tetiana's handouts**.
  Only matters when `K_frag ≈ 0`. Flagged for confirmation with Tetiana/Carrie, not
  a settled choice.

**d/dR (`beckmann/dft/descriptors.py::least_squares_slope()`):** ordinary
least-squares slope over the 5-point R(N-O) series (`slope =
sum((R_i-mean(R))*(y_i-mean(y))) / sum((R_i-mean(R))**2)`), computed for Psi,
log10(Lambda), wCNmax, w17max, w78max. The 5-point series uses stages
`["nbo", "sp2", "sp3", "sp4", "scan_2"]` -- `nbo` and `scan_1` are the same R0
geometry (see parse_nbo.py/parse_cmo.py docstrings), so `scan_1` is skipped to avoid
double-counting one point in the regression.

### Task 5 validation against `5_s0_Me.log` (Tetiana's reference log, compound 3 / Me)

Run via `scripts/analysis/validate_reference_descriptors.py`.

**Gotcha found and handled:** this log's route line includes `Stable=Opt`, which
triggers a wavefunction stability re-test and reruns population analysis -- the file
contains **two** full NBO/CMO sections (two `NBO 7.0` banners, two `Normal
termination` lines), not one. Always use the **last** occurrence of each, per the
handout's own warning.

**UPDATE (2026-07-16):** the line originally here said our own pipeline's `.gjf`
files don't use `Stable=Opt`, so this only mattered for reading this one external
file -- no longer true. The rigid-scan architecture's NBO blocks
(`_scan_gjf_rigid()`, `RIGID_SCAN_MIGRATION.md`) also use `Stable=Opt`, and this
exact gotcha (two tables per point, same R, only the last trustworthy) bit
`parse_nbo.py`/`parse_cmo.py` for real during the rigid-scan merge -- see
`RIGID_SCAN_MIGRATION.md`'s "MERGED" section for the full story and the fix
(both files' `parse_log()` now dedupes by R, keeping only the last table).

**Tier 1 (single-geometry, R0 only) -- PASSED:** `wCNmax = 0.457` at MO 32
(coefficient `-0.676`, `BD*(2) C7-N17`), exactly matching the wCNmax handout's worked
example. `K_anti`/`K_frag` sums verified by hand against the raw E2PERT table.
`Psi = 0.478` at R0 -- lower than 1 (fragmentation-channel stabilization exceeds
rearrangement-channel stabilization at equilibrium, even though compound 3
experimentally rearranges). Table 2 only reports `d/dR(Psi)`, not an absolute-value
threshold, and all 4 compounds in Table 2 show *positive* `d/dR(Psi)` regardless of
R/F outcome -- so a sub-1 Psi at R0 may not actually be surprising. Not confirmed
either way (Section 3.2's text wasn't read, out of scope for this task).

**Tier 2 (d/dR slope validation against Table 2) -- BLOCKED:** only `5_s0_Me.log`
(R0) exists in the repo for compound 3. Table 2's `d/dR` values need all 5 scan
points. Not approximated from one point -- the other 4 scan-point logs need to be
requested from Tetiana.

---

## Open issue: the wCNmax "extremum" test conflates two different mechanisms (not yet implemented)

`find_wcnmax_extremum()` in `scripts/analysis/summarize_descriptors.py` (called from
`main()`, feeds the `descriptor_summary.md` table) currently flags a "handoff" purely by
shape: is an interior point of the 5-point wCNmax(R) series a strict local min/max of its
two neighbors? A `depth` column was added later to distinguish a deep dip from a shallow
wobble, but depth alone still doesn't explain *why* a dip happens.

**Manually tracing every virtual MO that contributes to the target `BD*(C{ci}-N{ni})`
antibond (not just the single winning wcnmax MO) across all 6 test molecules
(mol_002/006/014/020/021/029)** shows the MO carrying `wcnmax` always switches identity
at the same scan step (between R0+0.1 and R0+0.2), but two genuinely different things can
happen there:

- **Real near-degenerate mixing** — the outgoing and incoming MOs are *both* printed
  simultaneously in the CMO table at that one geometry (i.e. both still clear NBO's 5%
  print threshold), with a small energy gap between them, and their combined weight is
  roughly conserved (character splits across two orbitals rather than vanishing):
  - `mol_014_Z` (exp = **F**): MO44 (weight 0.0625) + MO45 (weight 0.3612) at R0+0.2,
    ΔE = 0.0088 a.u. — total ≈ 0.42, matching the single-orbital value on either side.
  - `mol_020_E` (exp = **R**): MO48 (weight 0.0538) + MO49 (weight 0.4058) at R0+0.2,
    ΔE = 0.0181 a.u.
  This is what a textbook avoided crossing looks like.
- **A clean relay, no coexistence** — the outgoing MO has already dropped below the 5%
  print threshold by the next geometry, so only the incoming MO is ever visible in the
  table:
  - `mol_002_E` (F) and `mol_006_E` (R): clean relay with **no dip at all** — the
    coefficient magnitude barely changes across the switch (e.g. 002: 0.662 → 0.663;
    006: 0.653 → 0.654), so the scalar max glides through smoothly.
  - `mol_021_E` (R) and `mol_029_E` (R): also a clean relay (no coexisting second MO
    found), but with a small dip anyway (depth 0.016 / 0.010) — likely just numerical
    wobble in exactly where the crossover geometry lands, not a physical mixing event.

The current yes/no (or yes/no + depth) test cannot tell a real avoided crossing apart from
a numerical wobble in a clean relay — both get flagged "yes". And **checked against
experimental outcome, neither version cleanly separates R from F on this 6-molecule set
either way**: real near-degenerate mixing occurs in mol_014_Z (F) *and* mol_020_E (R)
alike; the shallow wobble occurs in two R's (021, 029); no handoff at all occurs in one F
(002) and one R (006). This isn't presented as a fix that will crack R/F prediction — it's
a more mechanistically honest descriptor to bring to the supervisor than the current
scalar-dip proxy.

**Proposed approach (not implemented yet):**

1. In `beckmann/dft/parse_cmo.py`, generalize `max_weight_for_target()` (currently
   returns only the single best-matching MO) to expose the full sorted list of matches, so
   the runner-up MO is available, not just the winner. Keep `max_weight_for_target()`'s own
   return value byte-identical (it's just `matches[0]`) so `wX^max`/`Lambda`/`Psi` are
   unaffected.
2. Add `second_MO_index`, `second_epsilon_i_star`, `second_coefficient`, `second_weight`,
   and `mo_gap` (= `abs(epsilon - second_epsilon)`) columns to
   `cmo_channel_extraction.csv` (`EXTRACTION_FIELDS`), populated from the second-best match
   per `(mol, stage, channel)`.
3. In `summarize_descriptors.py`, add a `classify_handoff()` next to
   `find_wcnmax_extremum()`: find the scan step where `MO_index` changes; if the *previous*
   stage's `MO_index` shows up as `second_MO_index` at that step with a small `mo_gap`
   (empirically < 0.03 a.u. cleanly separates the two confirmed real cases, 0.0088/0.0181,
   from the two wobble cases, which have no coexisting second MO at all) → label
   `"near-degenerate mixing"`; else if `MO_index` changes at all → `"clean handoff"`; else
   → `"no handoff"`. Add this as a new column in `descriptor_summary.md`, additive to the
   existing extremum/depth columns.
4. Regenerate via `python scripts/dft/parse_cmo.py` then
   `python scripts/analysis/summarize_descriptors.py`, and check the result against the
   hand-derived answers above before trusting it on future substrates.

No test currently covers `parse_cmo.py` at all (`tests/test_descriptors.py` only covers
`get_substituent_map`/`least_squares_slope`) — add one for the second-match/gap extraction
when this is implemented.

**IMPORTANT UPDATE (see below): mol_006_E's "clean relay, no dip at all" classification
above was a scan-resolution artifact, not a genuine absence of a dip** — a real, sharp
interior minimum sits at R=1.6608 Å, landing almost exactly on the midpoint between this
grid's R0+0.1 (1.6108 Å) and R0+0.2 (1.7108 Å) samples, so the standard 0.1 Å grid stepped
directly over it without ever landing near it. This calls into question whether the other
five molecules' "no dip"/"clean relay" classifications above are similarly under-sampled
rather than genuinely smooth — none of them have been checked at finer resolution. Treat
every "no handoff"/"no dip" result above as unconfirmed until re-scanned at 0.05 Å (or
finer) before relying on it for anything beyond this document's own working notes.

---

## mol_006_E's missing wCNmax minimum: resolved — it was a scan-resolution problem (2026-07-15)

The supervisor's reference paper reports an interior wCNmax minimum for this compound
(confirmed to be the same substrate as mol_006_E — see below); our standard 5-point,
0.1 Å-step scan showed a smooth, monotonic increase instead
(0.4225 → 0.4264 → 0.4277 → 0.4290 → 0.4303 across R0 to R0+0.4 Å), no minimum anywhere.

**Investigation trail (each ruled out in turn before landing on the real cause):**

1. **Scan architecture** (Gaussian's native internal multi-point walk vs. this session's
   new independent rigid-displacement-per-point architecture, see
   `RIGID_SCAN_MIGRATION.md`) — ruled out. wCNmax is ~identical between the two
   architectures at every point both have (differences ~0.000-0.004, noise-level).
2. **Substrate misidentification** — ruled out. Rigorously confirmed (graph/ring-position
   analysis on both structures, not eyeballing) that mol_006_E is the same compound as the
   supervisor's reference log `5_s0_Me.log` (repo root): the methyl substituent sits
   exactly 3 ring-bonds from the aryl-fusion carbon in both, vs. 2 and 1 bonds for the
   benchmark's other two methyl-indanone positional isomers (mol_009_E, mol_013_E).
3. **Aryl/alkyl channel assignment** (`get_substituent_map()`) — ruled out. Verified
   correct via direct connectivity inspection of the optimized structure.
4. **The LUMO-to-LUMO+0.4 a.u. energy window** (the paper's own method, vs. this
   codebase's unrestricted full-virtual-manifold search, see the `parse_cmo.py` docstring
   fix earlier in this document) — ruled out via
   `scripts/analysis/compare_wcnmax_window.py`: windowed and unrestricted searches give
   byte-identical results (same MO, same weight) at every point for mol_006_E, so the
   unrestricted search isn't reaching past the window to find a different, hidden
   maximum — there's nothing there to hide.
5. **Basis set** (our `wB97XD/6-311+G(d,p)` vs. the ~120-basis-function set implied by
   `5_s0_Me.log`'s `genecp` route line, consistent with `6-31G(d)`) — a strong lead at
   first (our own R0 geometry gives wCNmax=0.4225 vs. 0.457 on the reference geometry, and
   the R0 N-O bond length itself differs by 0.046 Å between the two) — but the supervisor
   confirmed directly that her `oxime_001_scan.gjf` reference file (also in the repo root)
   was only a basis/method **sensitivity test**, and results for this compound type don't
   depend on it. Ruled out by her explicit statement, not further computation on our end.
   (The R0 geometry/bond-length discrepancy itself is still unexplained, separately from
   the missing-minimum question — see `RIGID_SCAN_MIGRATION.md` for the raw numbers if
   this needs revisiting.)
6. **Scan resolution** — confirmed. The supervisor suggested keeping our basis but scanning
   at 0.05 Å steps instead of 0.1 Å (doubling the point density over the same R0 to
   R0+0.4 Å range: 9 points total including R0, vs. 5). Ran via
   `_scan_gjf_rigid(..., step=0.05, n_points=8)` (new optional parameters on the
   rigid-scan architecture, `beckmann/dft/inputs.py`), same
   `wB97XD/6-311+G(d,p)` basis throughout, base geometry reused unchanged from the
   existing (already-converged) Stage 1. Directory:
   `data/output/dft_opt_finescan/mol_006_E_finescan/`.

**Result: a real, sharp interior minimum at R = 1.6608 Å, wCNmax = 0.3260** — roughly
0.10 lower than the neighboring points (0.4264 at R=1.6108, 0.4277 at R=1.7108), landing
almost exactly on the midpoint between those two 0.1 Å-grid samples. The full 9-point
series:

| R (Å) | wCNmax |
|---|---|
| 1.5108 (R0) | 0.4225 |
| 1.5608 | 0.4238 |
| 1.6108 | 0.4264 |
| **1.6608** | **0.3260** |
| 1.7108 | 0.4277 |
| 1.7608 | 0.4277 |
| 1.8108 | 0.4290 |
| 1.8608 | 0.4303 |
| 1.9108 | 0.4303 |

Verified genuine, not an artifact:
- The underlying geometry optimization at this point fully converged to a real stationary
  point (`Optimization completed` / `Stationary point found`, Maximum Force = 0.000008 a.u.,
  far under the 0.00045 threshold) — not a crashed or partially-optimized structure.
- The winning MO (MO 45) is the *same* orbital carrying wCNmax at every neighboring point
  too — this isn't a different, spurious orbital being picked up. Its coefficient just
  drops sharply right at this one geometry (+0.571 vs. ~0.65-0.66 on either side), a real
  dip in that orbital's C{ci}-N{ni} antibond character, not an MO-identity handoff
  artifact (contrast with the near-degenerate-mixing/clean-relay handoffs discussed above,
  which happen at a *different* MO index changing hands — this is the same MO throughout).

**Visualization**: `data/output/analysis/plots/mol006_finescan_wcnmax.png` — see
`scripts/analysis/plot_mol006_finescan.py`.

**Takeaway for the rest of the benchmark**: a 0.1 Å-spaced 5-point scan is not fine enough
to reliably detect a real interior wCNmax minimum if it's narrow — this one would have
been completely invisible without the supervisor's suggested resolution increase. Every
other molecule's "no minimum" result in this document (and in `descriptor_summary.md`)
should be treated as unconfirmed at standard resolution, not as evidence the mechanism
genuinely doesn't occur there.

---

## Implement Orbital Resolved Electron Routing Framework

*(Original task prompt, kept for history. The prose descriptions of Λ/Ψ/wCNmax below
are ambiguous as literal formulas -- an earlier version of the code got Λ wrong as a
result. See "Corrected descriptor formulas" above for the precise, paper-matched
definitions actually implemented now.)*

Goal: Implement orbital resolved electron routing framework: move beyond single point ground state analysis and perform relaxed potential energy surface (PES) scans to capture electronic reorganization preceding bond cleavage. Selective rearrangement is determined by a specific avoided crossing event in the virtual manifold as N-O bond elongates. 

N-OH2 bond stretch: Gaussian will perform an optimization at the initial distance, then increment by 0.1 Å and re-optimize the rest of the molecule for each of the 5 snapshots
- relaxed scan, where the N–O bond is fixed at specific lengths while all other internal coordinates are optimized
- Gaussian **Opt=ModRedundant** keyword
- **The Geometry Section:** After the molecular coordinates, specify the bond to be stretched. 
- Important consideration: when the N-O bond is stretched the since oxygen is protonated, N-OH2 hydrogens should move along with oxygen not be strained away by mistake

Steps:
1. **Relaxed Potential Energy Surface (PES) Scan:** Perform a scan of the N–O bond 
2. **NBO Analysis on Snapshots:** Run NBO7 on the optimized geometry of each scan point using the following command string in your Gaussian input: `$NBO E2PERT BNDIDX NBOSUM CMO $END`.
3. **Data Extraction (Parsing):**
    - **From** **E2PERT** **/** **NBOSUM****:** Extract E(2) values for donor → acceptor interactions involving the activation coordinate (σNO∗​), the rearrangement channel (σC1−C7∗​), and the nitrilium channel (σC7−N17∗​).
    - **From** **CMO****:** Extract the leading NBO contributions (antibonding weights w) for all virtual orbitals within an energy window of **0.4 a.u. above the LUMO**.
4. **Metric Calculation:** Use your code to compute Ψ, logΛ, and wCNmax​ for each scan point.
5. **Differential Response:** Calculate the **least-squares slopes (**d/dR**)** of these metrics with respect to the N–O distance (R) to measure the rate of electronic reorganization

Analysis: Parse the NBO output at each scan point (N-O bond lengths R, R+0.1 A, R+0.2 A, R+0.3 A, R+0.4 A) to calculate these descriptors:
	- Hyperconjugative Competition (Ψ)
	- Frontier Dominance Metric (Λ)
	- CN-weighted Acceptor Response (wCNmax​)
	- Differential Response (d/dR)
Descriptors: key descriptors used in the Beckmann rearrangement study are custom metrics derived from raw NBO7 output. 
- **Hyperconjugative Competition (**Ψ**):** Calculated by taking the ratio of specific E(2) stabilization energies (e.g., the migrating bond feeding the activation coordinate).
- **Frontier Dominance (**Λ**):** A dimensionless measure calculated from the maximum antibonding weights (w) found in the virtual manifold.
- **CN-weighted Acceptor Response (**wCNmax​**):** Derived by parsing the Canonical Molecular Orbital (CMO) analysis to find the highest weight of specific antibonds (like σC1−C7∗​) within a specific energy window above the LUMO.
