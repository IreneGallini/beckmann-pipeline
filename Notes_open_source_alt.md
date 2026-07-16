# Open-source wCNmax prototype (branch: open-source-wcnmax-prototype)

Exploratory only. Nothing here touches `beckmann/dft/` or the validated main-branch
pipeline. Main branch's NBO7-derived numbers remain the trusted ground truth throughout
-- this branch asks "do these open-source alternatives rank/trend the same way," not
"are these numbers more correct."

New code lives entirely in `beckmann_alt/`:
- `geometry.py` -- reference-case geometries/atom maps, reusing
  `beckmann.dft.scan.parse_standard_orientations`/`ATOMIC_SYMBOLS` for the main-pipeline
  case and a local fallback for the external reference log's `nosymm` quirk (see below).
- `pyscf_livvo.py` -- Task 1: PySCF + LIVVO (`pyscf.lo.vvo`).
- `ao_projection.py` -- Task 2: crude AO-projection fallback.
- `compare.py` -- Task 3: both prototypes vs. the trusted NBO7 numbers.

## Reference cases

**mol_002_E** (primary): our own pipeline, level of theory fully known
(`wB97XD/6-311+G(d,p)`, `scrf=(smd,solvent=water)`, charge=1, singlet -- `beckmann/config.py`).
Geometry is the converged Stage-1 DFT geometry (`mol_002_E_opt.log`'s final "Standard
orientation"), the same geometry the trusted NBO7 numbers were computed on. Atom map
(`ci=11, ni=12, oi=13, c_aryl=6, c_alkyl=10`) reused from `beckmann.dft.descriptors`/
`beckmann.dft.scan` conventions, not re-derived.

**5_s0_Me** (secondary, explicitly NOT apples-to-apples): Tetiana's external reference
log (compound 3, "Ring Size and Substituent Effects in the Beckmann Rearrangement",
Table 2). Its actual route line is `wb97xd/genecp scrf=(smd,solvent=water) opt
scf=(Tight,xqc) nosymm`. Two things found before writing any PySCF code:
- **`genecp` is real, not incidental** -- the log's "Pseudopotential Parameters" section
  shows a custom hand-specified ECP applied to every C/N/O center (4 valence electrons
  per carbon, minimal split-valence primitives), not our all-electron `6-311+G(d,p)`.
  Reproducing it exactly would mean transcribing every exponent/coefficient/ECP
  parameter from the log by hand. Not attempted -- this case is run at **our own**
  `6-311+G(d,p)` instead, so any comparison against it is a different-basis check, not
  a validated reproduction of the paper's numbers.
- **`nosymm` changes the geometry block header** from `Standard orientation:` to `Input
  orientation:` -- `beckmann.dft.scan.parse_standard_orientations` (correct and
  unmodified for every log our own pipeline produces, which never uses `nosymm`) finds
  nothing in this file. `beckmann_alt/geometry.py` falls back to a locally-defined
  parser for `Input orientation:` blocks for this one external file's quirk only.

Trusted NBO7 numbers for both cases: wCNmax=0.457 (MO32, compound 3) is the paper's own
worked example, already validated by `scripts/analysis/validate_reference_descriptors.py`
(Tier 1, PASS). mol_002's numbers come straight from the main branch's already-generated
`cmo_channel_extraction.csv` (nbo/R0 stage). Neither was recomputed here.

## What had to change to actually run PySCF (found while building, not assumed upfront)

- **`wb97x-d` is blacklisted by PySCF's high-level dispersion dispatch**
  (`pyscf.scf.dispersion._black_list`). That blacklist exists to stop PySCF silently
  stacking an extra D3/D4 correction on a functional that already has its own dispersion
  built in -- exactly Gaussian's wB97XD. libxc itself recognizes `'wb97x-d'` fine
  (confirmed directly via `libxc.xc_type`). Removed it from the blacklist locally in
  `pyscf_livvo._enable_wb97xd()` rather than substituting PySCF's supported
  `wb97x-d3bj` variant, which uses a different (D3(BJ)) dispersion correction than
  Gaussian's wB97XD -- that substitution would be a bigger, less transparent deviation
  than bypassing an overcautious string check for a functional libxc already implements.
- **SMD's Python API exists but the compiled CDS extension does not.**
  `pyscf.solvent.smd.SMD` and its solvent parameter table (water entry present) are
  real, but the non-electrostatic cavity/dispersion/solvent-structure term requires a
  compiled `libsolvent` extension built with `-DENABLE_SMD=ON`, which this pip wheel
  does not include (`RuntimeError: SMD module is not available`). Building pyscf from
  source for this is a bigger ask than the plain `pip install` that was approved, so
  both prototypes use **ddCOSMO** instead. ddCOSMO is materially different from SMD --
  different cavity construction (fixed atomic-radius spheres vs. SMD's own tessellated
  solvent-accessible surface), and no explicit cavitation/dispersion/solvent-structure
  term at all. Treat the solvent contribution as a rough implicit-solvent stand-in, not
  a reproduction of Gaussian's actual SMD/water calculation.
- **Density fitting enabled for tractable runtime.** A conventional (non-DF) RSH-hybrid
  SCF on mol_002's 358-basis-function system did not finish in a reasonable time on this
  machine. `mf.density_fit()` is a standard, well-controlled approximation to the
  two-electron integrals -- flagged for completeness, not expected to be a significant
  source of divergence, but it is a deviation from Gaussian's conventional SCF.
- **Frozen geometry, no PySCF geometry optimization.** Both prototypes run a
  single-point calculation on the geometry Gaussian already converged to
  (`beckmann_alt/geometry.py`). This isolates "does the localization/projection method
  reproduce similar orbital character" from "does PySCF's optimizer land somewhere
  else," which isn't the question being asked here.

## LIVVO channel-identification rule

For each LIVVO (`pyscf.lo.vvo.vvo`, built from IAOs on the occupied space), compute its
Mulliken population on the two atoms defining the target antibond (`ci+ni` for CN,
`ci+c_aryl` for 17, `ci+c_alkyl` for 78) -- atom roles reused from
`beckmann.dft.descriptors`/`beckmann.analysis.classical`, not re-derived. The channel's
orbital is whichever LIVVO has the largest **combined** population on exactly that atom
pair, subject to a minimum acceptance threshold (`CHANNEL_POP_THRESHOLD = 0.5`,
`pyscf_livvo.py`) so a LIVVO that merely touches one of the two atoms isn't picked by a
low-confidence argmax. `compare.py`'s output reports explicitly whenever a channel is
rejected for falling below this threshold, rather than silently returning a
low-confidence answer.

## Results

All four runs (LIVVO x2, crude AO-projection x2) completed. Trusted NBO7 numbers are
pulled from existing artifacts, not recomputed: mol_002 from the main branch's already-
generated `cmo_channel_extraction.csv` (nbo/R0 stage); 5_s0_Me by re-running the same
already-validated `beckmann.dft.parse_cmo` extraction functions
`validate_reference_descriptors.py` itself uses (Tier 1, already PASS-validated against
the paper's own worked example).

### mol_002_E

| channel | NBO7 (trusted) | LIVVO | crude AO-projection |
|---|---|---|---|
| wCNmax | **0.4356** (MO48) | 0.0990 (MO112) | 0.0159 (MO48) |
| w17max | 0.1722 (MO124) | 0.0990 (MO112) | 0.0377 (MO48) |
| w78max | 0.1037 (MO118) | 0.0748 (MO126) | 0.0366 (MO51) |

NBO7 ranking: **cn > 17 > 78**, and cn is the clear standout (~2.5x w17, ~4.2x w78).

- **LIVVO**: identified the *same* LIVVO (#29) as the best match for both the cn and
  17 channels (combined Mulliken population 0.570 for cn, 0.414 for 17 -- below the
  0.5 acceptance threshold for 17, flagged `REJECTED` in the raw output but the number
  is still reported), so wCNmax and w17max come out numerically identical (0.0990).
  LIVVO does **not** reproduce NBO7's cn-dominant pattern here -- it ties cn with 17
  instead of ranking cn clearly on top. w78max (0.0748, its own distinct LIVVO #9,
  combined population 0.411, also below threshold) is the smallest of the three, which
  at least matches NBO7's ordering of 78 as smallest.
- **Crude AO-projection**: wCNmax's winning canonical MO is **MO48** -- the *exact*
  same MO NBO7 itself identifies as the wCNmax carrier. That's a real, specific
  agreement on orbital identity, not noise. But the *magnitude* ranking is inverted:
  w17max (0.0377) > w78max (0.0366) > wCNmax (0.0159) -- crude AO-projection ranks the
  CN channel **last**, the opposite of NBO7's cn-dominant pattern. w17max also lands on
  MO48 (same MO as wCNmax), so the two channels aren't well separated by this method
  either.

### 5_s0_Me (run at our basis, not the paper's actual GenECP basis -- see caveat above)

| channel | NBO7 (trusted, paper's own basis) | LIVVO (our basis) | crude AO-projection (our basis) |
|---|---|---|---|
| wCNmax | **0.4570** (MO32) | 0.1329 (MO43) | 0.0274 (MO53) |
| w17max | 0.0784 (MO47) | 0.1329 (MO43) | 0.0282 (MO44) |
| w78max | 0.0906 (MO39) | 0.1329 (MO43) | 0.0296 (MO44) |

NBO7 ranking: **cn > 78 > 17**, cn again the clear standout (~5x w78, ~5.8x w17).

- **LIVVO**: complete collision -- **all three channels** mapped to the same LIVVO
  (#9), all with sub-threshold combined population (0.417 / 0.414 / 0.362), giving
  identical wmax (0.1329) for cn, 17, and 78. No ranking information survives at all
  for this molecule -- the channel-identification step failed to distinguish any of
  the three antibonds from each other.
- **Crude AO-projection**: w78max (0.0296) > w17max (0.0282) > wCNmax (0.0274) -- again
  inverts NBO7's cn-dominant ranking, placing cn last, though the values are so close
  together (all within 0.002 of each other) that "ranking" is barely meaningful here.
  w17max and w78max share the same winning MO (44), so again not well separated.

### Plain assessment (not averaged, not cherry-picked)

**Neither prototype reproduces the central signature that makes wCNmax useful in the
current pipeline: that it is clearly, substantially larger than w17max/w78max.** NBO7
shows this pattern strongly in both reference cases (cn 2.5-5.8x the other channels).
LIVVO ties or completely collapses the three channels together in both cases. Crude
AO-projection actively **inverts** the ranking in both cases, placing the CN channel
last rather than first.

The specific, reproducible failure mode in both prototypes is the **channel-
identification step**, not the underlying projection math: NBO7's per-atom-pair Lewis-
structure search produces antibonds that are, by construction, sharply localized on
exactly two atoms. Neither a single-vector Mulliken-population argmax over LIVVOs, nor
a single hand-built p-orbital trial function, reliably finds an equally sharply-
localized two-center object -- both approaches most often land on some orbital with
substantial density spread across more of the ring system (the LIVVO/crude vectors
matching multiple different target channels to the *same* winning MO, in half of the
eight results above, is the direct symptom of this).

One genuinely interesting positive result: crude AO-projection's wCNmax for mol_002
converges on the *identical* canonical MO (MO48) that the trusted NBO7 pipeline
independently identifies -- suggesting the crude trial function does capture something
real about that specific orbital's character, even though it fails to rank its
*magnitude* correctly relative to the other channels.

### Follow-up: does a better channel-identification metric fix it?

The plain assessment above pointed at the channel-identification step (picking *which*
LIVVO represents a given antibond) as the bottleneck, not the projection math. The
natural next thing to try: `identify_channel_livvo`'s scoring rule sums each atom's
*diagonal* Mulliken population (`pop[a] + pop[b]`), which can't see whether atom A's
and atom B's contributions to a given orbital are in-phase or out-of-phase -- two atoms
that both happen to carry a lot of density in a delocalized orbital score identically
to a genuine antibonding combination between them. `identify_channel_livvo_v2`
(`beckmann_alt/pyscf_livvo.py`) replaces this with the **interatomic (off-diagonal)
Mulliken overlap population** between the two atoms -- the standard bonding/
antibonding diagnostic (large negative = antibonding character concentrated
specifically between A and B) -- and picks the LIVVO with the most negative value.
Tested via cached SCF/LIVVO arrays (no new DFT needed) against both reference cases:

**mol_002_E**: v2 *does* separate cn from w17 (previously both landed on livvo#29;
now cn stays on #29, w17 moves to a different orbital, #25) -- real progress on that
specific collision. But w78 then collides with the *new* w17 pick (both land on #25,
where v1 had kept w78 separate). And the resulting ranking (w17 = w78 > cn) still
doesn't reproduce NBO7's cn-dominant pattern -- if anything it's worse than v1 on this
specific point, since v1 at least tied cn for first place; v2 puts cn in last place.

**5_s0_Me**: v2 gives **no improvement at all** -- all three channels still collide,
just onto a different single LIVVO (#26 instead of #9) instead of being separated.

**Conclusion**: a better scoring metric is not the fix. Both scoring rules are
searching the *same* fixed set of ~29-30 LIVVOs (one global SVD rotation of the whole
virtual space against IAOs, computed once per molecule) -- if that set doesn't happen
to contain three orbitals each cleanly localized on a different one of the three target
atom pairs, no post-hoc scoring rule over that fixed set can conjure a fourth,
better-localized option into existence. This points at something more fundamental than
"pick a better argmax criterion": LIVVO's single, non-iterative SVD construction is
structurally different from NBO's actual algorithm (Step 4 in the reference-material
breakdown at the top of this file), which builds a *fresh*, small 2-center subspace
directly from each atom pair's own block of the density matrix, one pair at a time,
subtracting occupancy as it goes. A real fix likely means building something closer to
that iterative per-atom-pair construction directly (a meaningfully larger undertaking
than swapping an identification formula), not choosing a different scoring rule on top
of LIVVO's existing output.

### Second follow-up: the local per-atom-pair construction (`beckmann_alt/pair_nbo.py`)

The diagnosis above pointed at something more fundamental than a scoring formula:
LIVVO draws from one fixed, global, non-iterative ~29-30-orbital set, and no scoring
rule over a fixed set can manufacture a better-localized orbital if the set doesn't
contain one. This module tests the actual implied next step: instead of picking from a
pre-built global set, build a **fresh, local subspace for each atom pair individually**,
directly from that pair's own block of the density matrix (in a Löwdin-orthogonalized
IAO basis) -- structurally the same operation as NBO's real Step 4, just for one
requested pair at a time rather than NBO's full iterative multi-pair deflation across
the whole molecule (see caveat below).

**First pass (single lowest-occupation eigenvector per pair):** immediately fixed the
collision problem completely -- cn/17/78 landed on three *distinct* MOs in both test
cases, something neither LIVVO variant achieved even once across four attempts. Bond
occupations came out ~2.0000 in every channel of both cases, confirming the local
blocks are finding genuine bonding pairs. w17max/w78max for mol_002 landed within 1 MO
index and ~10% of the trusted NBO7 magnitude -- clearly the best result of any method
tried so far on those two channels. But wCNmax stayed badly wrong (0.087 vs trusted
0.436 for mol_002; 0.117 vs 0.457 for 5_s0_Me -- consistently ~4-5x too small in both
cases, and only in the cn channel).

**Why, and the fix:** printing every local eigenvector's occupation for mol_002's cn=
(C11,N12) pair showed *two* low-occupation eigenvectors, not one: `[0.016, 0.219, 0.568,
1.02, ...]`. C=N is a double bond -- a real sigma*/pi* pair should exist locally, and
projecting the *second*-lowest candidate (occupation 0.219) instead of the lowest gave
wmax=0.460 at MO47, matching NBO7's trusted 0.436 at MO48 to within 6%. This is exactly
what the real wX^max definition already does ("both BD*(1) and BD*(2) ... whichever
gives the larger squared coefficient wins," Notes.md) -- the fix (`ANTIBOND_OCC_THRESHOLD`
in `pair_nbo.py`) is to treat every local eigenvector below occupation 1.0 as a
candidate, project all of them, and keep whichever wins, mirroring that rule exactly
rather than assuming the single lowest-occupation eigenvector is always the right one.

**Result after the fix:**

| channel | NBO7 (trusted) | mol_002 (this method) | 5_s0_Me trusted | 5_s0_Me (this method) |
|---|---|---|---|---|
| wCNmax | 0.4356 (MO48) | **0.4599** (MO47) | 0.4570 (MO32) | **0.4602** (MO43) |
| w17max | 0.1722 (MO124) | 0.2009 (MO57) | 0.0784 (MO47) | 0.2729 (MO49) |
| w78max | 0.1037 (MO118) | 0.2282 (MO47) | 0.0906 (MO39) | 0.2334 (MO43) |

**For the first time across every method tried in this exploration, wCNmax is both
correctly ranked as the dominant channel AND numerically close to the trusted value in
both reference cases** (within 5.6% for mol_002, within 0.7% for 5_s0_Me) -- LIVVO
never got closer than ~4x off, crude AO-projection never got closer than ~15x off and
actively ranked cn last. This is the single result this whole exploration was aimed at
reproducing, and it reproduces in both cases, not just one.

**The trade-off, reported plainly, not hidden:** broadening the candidate search to fix
wCNmax reintroduced collisions elsewhere -- w78 now lands on the exact same MO as
wCNmax in *both* cases (MO47 for mol_002, MO43 for 5_s0_Me), and w17/w78 are both
noticeably over-estimated relative to NBO7 (roughly 1.2-3.5x too large), with their
*relative* order (78 vs 17) not matching NBO7's in either case. Widening the candidate
pool per channel made it more likely that unrelated channels' searches converge on the
same broadly-delocalized virtual MO -- the same underlying tension as before (a small
set of "generically effective" virtual MOs competing against genuinely 2-center-local
character), just shifted from the LIVVO-selection step to the candidate-projection step.

**Caveat -- what this still does not reproduce:** NBO's real algorithm is iterative
across the *whole* molecule -- once a pair's bond/antibond is accepted, its occupancy is
subtracted from the density matrix before the next pair is searched, so two target
pairs sharing an atom (cn and w17 both include ci) don't draw on the same undiminished
density there. This module searches every pair completely independently with no
deflation between them. The w78/wCNmax MO collisions above are plausibly a direct
symptom of that missing deflation step -- worth testing directly if this continues.

**Bottom line for the open-source goal, updated**: the local per-atom-pair construction
is a substantial, genuine step forward -- it is the first approach in this exploration
that reproduces wCNmax's central diagnostic signature (dominant, and numerically close
to NBO7) in both reference cases. It is not a finished replacement: w17max/w78max are
still noticeably off and prone to colliding with each other or with wCNmax. The
highest-value next step, if this continues, is implementing the missing deflation
(subtracting each accepted candidate's density contribution before searching the next
pair) rather than trying yet another localization method from scratch -- this
follow-up's own diagnosis suggests that's precisely what's missing now.

## Crude AO-projection caveats (Task 2)

Expected to be the weaker of the two prototypes by construction: it sums every p-type
AO shell on each atom with equal weight (6-311+G(d,p) has several per atom -- core,
valence, diffuse), with no attempt to weight toward the shell that actually dominates a
real valence antibond, and performs no orthogonalization against the occupied space
(unlike the IAO/VVO route, which orthogonalizes by construction). Meant as a rough
same/different-ballpark sanity check, not comparably trustworthy to the LIVVO prototype
-- see `beckmann_alt/ao_projection.py`'s module docstring.
