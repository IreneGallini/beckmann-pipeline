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
- `pair_nbo.py` -- the local per-atom-pair follow-up (see "Second follow-up" below),
  the best-performing method and the only one actively maintained going forward.

**Scope decision (current)**: wCNmax is the actual predictive descriptor this project
needs; w17max/w78max are secondary and were never the target of this prototype effort.
`pair_nbo.py`'s code has been simplified to compute wCNmax only -- the w17/w78 channel
logic and the entire "Third follow-up: iterative deflation" section described below
(which existed specifically to try to fix a w17/w78 collision) have been **removed from
the current code**, since none of that machinery was needed for wCNmax's own result
(confirmed: wCNmax is computed independently of w17/w78 in every version tried, and in
the deflation version specifically, cn was always processed first, so nothing had been
deflated yet when its own value was computed -- the deflation follow-up contributed
nothing to wCNmax's accuracy). The sections below documenting that work are kept as the
honest historical record of what was tried and why -- available via
`git log -- beckmann_alt/pair_nbo.py` for the actual removed code -- not as a
description of `pair_nbo.py`'s current contents.

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

### Third follow-up: iterative deflation (`run_case_deflated` etc. in `beckmann_alt/pair_nbo.py`)

Implements what the previous section identified as missing: after a channel's antibond
is picked, subtract its accepted bond(s)' density from the working matrix before
searching the next channel, so `cn`/`w17`/`w78` (which share atom `ci`) don't all draw
on the same undiminished density. A Plan-agent design review beforehand confirmed the
core deflation formula (`dm_ao - occupation * |bond_vec><bond_vec|`, using the
*measured* occupation rather than hardcoding 2.0) is exact via a trace/electron-count
identity, and flagged two real risks to build diagnostics for: (1) deflating only the
single highest-occupation eigenvector per pair would leave a real gap for the C=N
double bond specifically (two near-2.0 eigenvectors, sigma- and pi-like, not one), and
(2) later channels' local blocks can develop negative eigenvalues (a rank-1-downdate
effect) after deflation -- not hypothetical, needed explicit checking.

**First implementation attempt failed badly, caught immediately by the PSD check the
review recommended.** Deflating every eigenvector above the same threshold used to
define antibond candidates (`ANTIBOND_OCC_THRESHOLD = 1.0`) produced severe PSD
violations -- local occupations down to **-0.56 to -0.60** for later channels, far past
noise. Cause: every local block has an "ambiguous middle" cluster of eigenvectors
(occupation roughly 0.5-1.4) that are not clean, fully-localized 2-electron bonds --
e.g. mol_002's `w17` block: `[0.03, 0.53, 0.84, 0.94, 1.00, 1.17, 1.37, 1.96, 2.00,
2.00]`. Treating everything above 1.0 as "an accepted bond, safe to fully subtract" was
a real over-correction, not a minor rounding issue.

**Fix**: a separate, stricter `BOND_OCC_THRESHOLD = 1.9` for what counts as a genuine
bond to deflate, picked from an empirical pattern that held identically across both
molecules and all three channels -- the ambiguous middle always tops out below ~1.4,
genuine bonds always cluster at ~1.96-2.00, with a clean gap in between every time.
Eigenvectors strictly between the two thresholds are left untouched (neither deflated
nor treated as an antibond candidate). This dropped the PSD violations to
borderline/noise-level (-0.0001 to -0.011, vs. the earlier -0.56 to -0.60) -- not
completely eliminated, but no longer a real correctness problem.

**Results, both orderings, both cases** (baseline = no deflation, from the Second
follow-up):

| | wCNmax | w17max | w78max |
|---|---|---|---|
| **mol_002** trusted | 0.4356 (MO48) | 0.1722 (MO124) | 0.1037 (MO118) |
| baseline (no deflation) | 0.4599 (MO47) | 0.2009 (MO57) | 0.2282 (MO47) |
| deflated, occ-order (cn,w78,w17) | 0.4599 (MO47) | 0.2048 (MO123) | 0.2636 (MO47) |
| deflated, cn-first (cn,w17,w78) | 0.4599 (MO47) | 0.1926 (MO123) | 0.2636 (MO47) |
| **5_s0_Me** trusted | 0.4570 (MO32) | 0.0784 (MO47) | 0.0906 (MO39) |
| baseline (no deflation) | 0.4602 (MO43) | 0.2729 (MO49) | 0.2334 (MO43) |
| deflated, occ-order (cn,w78,w17) | 0.4602 (MO43) | 0.2453 (MO49) | 0.2662 (MO43) |
| deflated, cn-first (cn,w17,w78) | 0.4602 (MO43) | 0.2270 (MO49) | 0.2662 (MO43) |

**Plain verdict: partially resolved, not fixed, and not uniformly.**

- **wCNmax: unaffected either way** (expected and desired -- cn is always processed
  first in both orderings tried, so nothing is deflated before its own search runs).
  Still the strongest result of this whole exploration.
- **mol_002's w17: real improvement.** Moved from colliding with no trusted MO in
  particular (MO57 baseline) to **MO123 in both orderings** -- 1 MO index from trusted
  MO124, and magnitude (0.19-0.20) much closer to trusted (0.172) than baseline's 0.20
  was already close, actually further tightened under cn-first ordering specifically.
- **w78: not fixed in either case.** It still lands on the *exact same MO as wCNmax*
  (MO47 for mol_002, MO43 for 5_s0_Me) in **every** deflated run tried, both orderings,
  both molecules -- identical to the undeflated baseline. Deflating cn's (and w17's)
  bonds did not stop w78's search from converging on the same broadly-effective virtual
  MO cn already claims. If anything its magnitude got *worse* (further from trusted in
  both cases: mol_002 0.228->0.264, 5_s0_Me 0.233->0.266).
- **5_s0_Me's w17: no real improvement.** Stays on the same MO (49) as baseline in both
  orderings; magnitude moves modestly (0.273->0.245/0.227) but remains ~3x the trusted
  0.0784 -- deflation trimmed the overestimate slightly without addressing it.
- **Ordering does measurably change the numbers** (e.g. mol_002 w17: 0.2048 occ-order
  vs 0.1926 cn-first, a genuine ~6% difference from processing w17 before vs. after
  w78), confirming deflation is order-dependent as expected -- but the two orderings
  landed on the *same* winning MO in every channel in both cases, and neither ordering
  resolved the w78/cn collision. In this data, `bond_occupation` for all three channels
  sat within 0.0001 of 2.0000 in both molecules, so decreasing-occupation order was
  effectively a coin flip here (exactly what the design review anticipated) -- it did
  not turn out to be a more informative choice than cn-first in practice, though there
  wasn't a principled reason to expect it to be either.

**What this suggests**: the w78/cn collision looks less like "shared, undeflated
density" and more like a genuine structural preference of this construction --
something about w78's local block (ci, c_alkyl) keeps resolving toward the same
virtual MO cn's block resolves toward, independent of what's been deflated beforehand.
The design review's alternative idea (excluding already-claimed MO indices at the
projection step, rather than touching the density matrix at all) was deliberately not
implemented this pass, but this result makes it a more interesting thing to try next
than a different deflation order or threshold -- it would directly test whether the
collision is really about density-sharing (which deflation targets and mostly didn't
fix) or something else entirely.

## Crude AO-projection caveats (Task 2)

Expected to be the weaker of the two prototypes by construction: it sums every p-type
AO shell on each atom with equal weight (6-311+G(d,p) has several per atom -- core,
valence, diffuse), with no attempt to weight toward the shell that actually dominates a
real valence antibond, and performs no orthogonalization against the occupied space
(unlike the IAO/VVO route, which orthogonalizes by construction). Meant as a rough
same/different-ballpark sanity check, not comparably trustworthy to the LIVVO prototype
-- see `beckmann_alt/ao_projection.py`'s module docstring.

## Fourth check: wCNmax across all 6 test-set molecules

After the decision to scope `pair_nbo.py` down to wCNmax only (see the scope-decision
note near the top of this file), tested against every main-pipeline test-set molecule
with completed NBO7 data (`beckmann.dft.inputs.TEST_IDS`), not just mol_002 -- using
`beckmann_alt.geometry.load_test_set_case()`, which resolves atom maps and geometry
fresh via `resolve_mol_name`/`oxime_atom_map_from_gjf` (the main pipeline's own
utilities) rather than hand-transcribing six atom maps.

| mol | computed wCNmax (MO) | trusted wCNmax (MO) | % diff | MO offset |
|---|---|---|---|---|
| mol_002_E | 0.4599 (47) | 0.4356 (48) | +5.6% | -1 |
| mol_006_E | 0.4473 (43) | 0.4225 (44) | +5.9% | -1 |
| mol_014_Z | 0.4431 (43) | 0.4212 (44) | +5.2% | -1 |
| mol_020_E | 0.4781 (47) | 0.4502 (48) | +6.2% | -1 |
| mol_021_E | 0.4944 (51) | 0.4692 (52) | +5.4% | -1 |
| mol_029_Z | 0.4933 (47) | 0.4665 (48) | +5.7% | -1 |
| 5_s0_Me | 0.4602 (43) | 0.4570 (32) | +0.7% | -11 (different basis, see caveat above -- not directly comparable) |

**The pattern is tight and systematic, not noisy.** Across all 6 of our own-pipeline
molecules (same basis, same functional, same solvent-caveat throughout): the percent
error sits in a narrow 5.2-6.2% band (mean ~5.7%), **always an overestimate, never an
underestimate**, and the winning canonical MO is **exactly 1 index lower than NBO7's own
winning MO in every single case** -- not "close," not "usually," literally 6/6. This
level of consistency across six chemically different substrates (different ring sizes,
substituents, and connectivity) is much more characteristic of a fixed, structural
offset (e.g. a systematic difference in how PySCF's and Gaussian's virtual manifolds are
ordered/counted near the frontier, or a consistent solvent/basis/dispersion-driven
energy shift that happens to move exactly one virtual MO's relative position) than of
six independent, coincidentally-similar numerical errors. Not root-caused yet -- would
need a direct comparison of virtual orbital energies near the frontier between the
Gaussian and PySCF calculations for the same molecule to pin down which piece (solvent
model, dispersion treatment, or something else) produces the shift. If diagnosed, a
fixed 5-6%/1-MO-index systematic offset would be straightforward to calibrate out;
right now it's reported as an observed, reproducible pattern, not yet an explained or
corrected one.

`5_s0_Me` sits apart, as expected given its different-basis caveat -- its % error is
much smaller (0.7%) but its MO offset is much larger (-11), underscoring that magnitude
agreement and orbital-identity agreement are answering different questions here, and
that `5_s0_Me`'s numbers shouldn't be pooled with the other 6 in any summary statistic.
