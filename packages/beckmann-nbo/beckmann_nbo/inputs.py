"""
Prepare Gaussian 16 input files for DFT/NBO analysis.

Two workflows:
  prepare_opt / main_opt — three-stage workflow for the test set (opt + NBO + scan)
  prepare_sp  / main_sp  — single-point NBO directly on AIMNet2 geometry (all 34 molecules)
"""
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from pathlib import Path
from rdkit import Chem

from beckmann_core.constants import CHARGE, MULTIPLICITY
from beckmann_core.geometry import displace_leaving_group
from beckmann_nbo.config import (
    DATA_OUTPUT,
    FUNCTIONAL, BASIS, NPROC, MEM_GB,
    NBO_KEYWORDS, SOLVENT,
)

TEST_IDS  = {"002", "006", "020", "021", "014", "029"}  # original test-set subset, kept for fast iteration
OXIME_PAT = Chem.MolFromSmarts('[C:1]=[N:2]-[O+:3]')


def _all_ids() -> set[str]:
    """Numeric IDs for every substrate in the benchmark set (e.g. '001'..'034'),
    derived from benchmark_meta.json rather than hand-maintained, so scope
    tracks the benchmark set instead of drifting from it."""
    import json
    from beckmann_nbo.config import DATA_INPUT
    meta = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())
    return {key.split("_")[1] for key in meta}


ALL_IDS = _all_ids()

# Three molecules whose canonical dft_opt/{mol}/{mol}_scan.gjf is a crashed
# Stage 3 attempt (see JOB_ISSUES.md, 2026-07-20 entries) but which have one
# or more successful reruns at a different scan step size, living as sibling
# directories under dft_opt_stepscan/. mol_020_E and mol_034_E succeeded at
# *both* listed step sizes -- both are kept and merged into one denser
# per-molecule series (see merge_scan_rows()/build_stage_relabel_map() below)
# rather than picking one as canonical, since the scan's purpose is resolving
# whether wCNmax is monotonic or has an interior minimum, and more points is
# strictly better for that. mol_003_E only has one successful rerun (step04).
STEP_SCAN_SOURCES: dict[str, list[str]] = {
    "mol_003_E": ["mol_003_E_step04"],
    "mol_020_E": ["mol_020_E_step07", "mol_020_E_step04"],
    "mol_034_E": ["mol_034_E_step07", "mol_034_E_step04"],
}


def step_scan_dir() -> Path:
    return DATA_OUTPUT / "dft_opt_stepscan"


def resolve_mol_name(mol_id: str, dft_opt_dir: Path) -> str | None:
    """Find the isomer-suffixed mol dir name for a numeric id (e.g. '014' ->
    'mol_014_Z') -- the AIMNet2-lower-energy isomer isn't always E."""
    matches = sorted(dft_opt_dir.glob(f"mol_{mol_id.zfill(3)}_*"))
    return matches[0].name if matches else None


def build_stage_relabel_map(r_no_values: set) -> dict:
    """Map each unique (rounded) r_no to a fresh 'scan_i' label, i=1.. in
    ascending R order -- used to renumber scan-point rows collected from
    more than one source log (e.g. mol_020_E's step07 + step04 reruns) into
    one coherent series. 'nbo' (R0) rows aren't included -- they pass
    through unrenumbered since there's exactly one baseline point."""
    ordered = sorted(r for r in r_no_values if r is not None)
    return {r: f"scan_{i}" for i, r in enumerate(ordered, start=1)}


def relabel_rows(rows: list[dict], mol_name: str, relabel: dict,
                  stage_key: str = "stage", r_no_key: str = "r_no") -> list[dict]:
    """Rewrite every row's mol field to mol_name; rows whose stage_key is
    'nbo' pass through as-is (relabeled mol only), everything else gets its
    stage_key rewritten via relabel (a build_stage_relabel_map() result,
    keyed by rounded r_no) so scan points collected from multiple source
    logs land as one 'scan_1'..'scan_N' sequence sorted by actual R(N-O).
    r_no_key lets callers point at a differently-cased/named R(N-O) field
    (e.g. parse_cmo.py's channel-extraction rows use 'R_NO', not 'r_no')."""
    out = []
    for row in rows:
        if row[stage_key] == "nbo":
            out.append({**row, "mol": mol_name})
            continue
        key = round(row[r_no_key], 4) if row[r_no_key] is not None else None
        out.append({**row, "mol": mol_name, stage_key: relabel[key]})
    return out


# ── three-stage opt workflow ────────────────────────────────────────────────────

def _opt_gjf(name: str, coords: list[tuple], oxime_label: str, basis: str = BASIS, calcfc: bool = False) -> str:
    """Stage 1: geometry optimisation no NBO block.

    calcfc=True switches the route line to opt=(CalcFC,MaxCycles=300) --
    see _scan_gjf_rigid()'s docstring for the full rationale (fused-ring
    pucker oscillation fix, JOB_ISSUES.md). Used when a molecule's Stage 1
    itself oscillates (e.g. mol_016_E), not just Stage 3 scan points.
    """
    opt_kw = "opt=(CalcFC,MaxCycles=300)" if calcfc else "opt"
    return (
        f"%chk={name}_opt.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{basis} {opt_kw} {SOLVENT}\n"
        f"\n"
        f"{name} opt  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        + "\n".join(
            f"{sym:<3}  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
            for sym, x, y, z in coords
        )
        + "\n\n\n"
    )


def _nbo_gjf(name: str, oxime_label: str, basis: str = BASIS) -> str:
    """Stage 2: NBO7 single-point at DFT geometry.

    pop=nbo7read (not nboread) routes through Gaussian's external-program
    interface (Link 612 -> gaunbo7 -> g16nbo -> nbo7), which is required for
    CMO-based descriptors (Lambda, wCNmax) the bundled NBO 3.1 (pop=nboread)
    doesn't support CMO at all.
    """
    return (
        f"%chk={name}_nbo.chk\n"
        f"%oldchk={name}_opt.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{basis} sp pop=nbo7read geom=checkpoint guess=read {SOLVENT}\n"
        f"\n"
        f"{name} NBO  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        f"\n"
        f"$NBO {NBO_KEYWORDS} $END\n"
        f"\n\n"
    )


def _scan_gjf(name: str, ni: int, oi: int, oxime_label: str) -> str:
    """LEGACY (old internal-walk scan architecture -- see RIGID_SCAN_MIGRATION.md;
    superseded by _scan_gjf_rigid()/prepare_scan_rigid() below, kept for
    reference/rollback, not called by prepare_opt()). Stage 3: relaxed N-O
    bond scan 5 points (R to R+0.4 Å) via Gaussian's own internal multi-point
    walk -- only ran full NBO7 at 2 of the 5 points (R0 and R0+0.4)."""
    return (
        f"%chk={name}_scan.chk\n"
        f"%oldchk={name}_opt.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} opt=(ModRedundant,MaxCycles=200) pop=nbo7read geom=checkpoint guess=read {SOLVENT}\n"
        f"\n"
        f"{name} scan  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
        f"\n"
        f"B {ni} {oi} S 4 0.1\n"
        f"\n"
        f"$NBO {NBO_KEYWORDS} $END\n"
        f"\n\n"
    )


def _scan_gjf_rigid(
    name: str, base_atoms: list, ni: int, oi: int, oxime_label: str,
    basis: str = BASIS, step: float = 0.05, n_points: int = 6, calcfc: bool = False,
) -> str:
    """Stage 3 (rigid-scan architecture): n_points independent points
    (R0+step .. R0+step*n_points), each built from the SAME Stage-1 base
    geometry -- not chained point-to-point, not extracted from an internal
    Gaussian scan walk. Per point: a rigid O-N displacement (pure geometry,
    via displace_leaving_group), then a constrained optimization
    (opt=ModRedundant, bond frozen at the new length, everything else
    relaxes), then a same-checkpoint Link1 continuation into a full
    NBO7/CMO single point. Matches the PI's reference file
    (oxime_001_scan.gjf) -- see JOB_ISSUES.md for the full rationale.

    Defaults (step=0.05, n_points=6) give the standard 7-point series
    (R0..R0+0.30 A, R0 covered separately by Stage 2). This replaced an
    earlier (step=0.1, n_points=4)/R0..R0+0.4 A default after mol_006_E's
    finer 0.05 A rerun showed the coarser 0.1 A grid can step directly over
    a real, sharp interior wCNmax minimum (see Notes.md). The range is
    truncated at R0+0.30 A rather than extended to R0+0.4 A for two
    reasons: mol_006_E's real minimum landed well inside R0+0.30 A, and the
    far end of the range (R0+0.35/+0.4 A) has repeatedly been the least
    stable part of the scan -- JOB_ISSUES.md documents a genuine
    non-converging double-well oscillation at R0+0.4 A for mol_020_E and
    outright crashes/noisy convergence at R0+0.4 A for mol_006_E under two
    earlier architectures.

    NBO keywords deliberately omit her NBOMO=P120 print-window restriction:
    parse_cmo.py was fixed earlier to search the entire virtual manifold
    unconditionally (see its docstring) after discovering real target
    antibonds mixing in above a narrower window for some substrates -- a
    fixed narrow NBOMO range would reintroduce that at the Gaussian-printing
    level, upstream of anything the parser can recover.

    calcfc=True applies opt=(ModRedundant,CalcFC,MaxCycles=300) to every
    point's opt block uniformly, instead of the default opt=(ModRedundant)
    SCF=(Tight,XQC) NoSymm -- see JOB_ISSUES.md's 2026-07-16/18 entry. This
    was originally a per-crashed-point-only patch, which got reverted for
    creating a methodological inconsistency within one molecule's series
    (5 points under one setting, 1 under another). This parameter instead
    applies it molecule-wide -- if any point in a series needs it, rerun
    the whole series with it -- being tested as a candidate default for
    molecules with a known/suspected fused-ring pucker oscillation risk.
    """
    opt_kw = "opt=(ModRedundant,CalcFC,MaxCycles=300)" if calcfc else "opt=(ModRedundant)"
    blocks = []
    for pt in range(1, n_points + 1):
        delta = pt * step
        displaced = displace_leaving_group(base_atoms, ni, oi, delta)
        coord_block = "\n".join(
            f"{sym:<3}  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}" for sym, x, y, z in displaced
        )
        chk = f"{name}_scan_pt{pt}.chk"
        # No %oldchk -- literal coordinates are given directly below (not
        # geom=check), so there's nothing to copy from another checkpoint,
        # matching the reference file (oxime_001_scan.gjf never uses %oldchk).
        blocks.append(
            f"%chk={chk}\n"
            f"%nprocshared={NPROC}\n"
            f"%mem={MEM_GB}GB\n"
            f"#p {FUNCTIONAL}/{basis} {opt_kw} SCF=(Tight,XQC) NoSymm {SOLVENT}\n"
            f"\n"
            f"{name} scan pt{pt} (R0+{delta:.2f}A) rigid O-N displacement then constrained opt  {oxime_label}\n"
            f"\n"
            f"{CHARGE} {MULTIPLICITY}\n"
            f"{coord_block}\n"
            f"\n"
            f"B {ni} {oi} F\n"
        )
        blocks.append(
            f"%chk={chk}\n"
            f"%nprocshared={NPROC}\n"
            f"%mem={MEM_GB}GB\n"
            f"#p {FUNCTIONAL}/{basis} Geom=Check Guess=Read Stable=Opt Pop=NBO7Read Density=Current {SOLVENT}\n"
            f"\n"
            f"{name} SP+NBO7/CMO after constrained opt at R0+{delta:.2f}A  {oxime_label}\n"
            f"\n"
            f"{CHARGE} {MULTIPLICITY}\n"
            f"\n"
            f"$NBO CMO PRINT=2 E2PERT=0.05 BNDIDX $END\n"
        )
    return "\n--Link1--\n".join(blocks) + "\n\n"


def prepare_opt(
    sdf_path: Path,
    outdir: Path,
    test_ids: set[str] = TEST_IDS,
) -> None:
    """Write _opt.gjf and _nbo.gjf for each molecule in test_ids (Stages 1-2
    only). Stage 3 (_scan.gjf, rigid-scan architecture) can't be generated
    here -- _scan_gjf_rigid() needs Stage 1's *converged* geometry, which
    doesn't exist until Stage 1 has actually run on Citadel and its log is
    downloaded. Call prepare_scan_rigid() as a separate step once that's
    done (see its docstring) -- mirrors the same "generate after the
    previous stage completes" pattern extract_scan_sp.py used under the old
    architecture (see RIGID_SCAN_MIGRATION.md)."""
    outdir.mkdir(parents=True, exist_ok=True)

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mols  = [m for m in suppl if m is not None]
    test_mols = [m for m in mols if m.GetProp("_Name").split("_")[1] in test_ids]

    print(f"\n{'Name':<24} {'Atoms':>5}  {'Oxime':>20}  Stage1  Stage2")
    print("-" * 70)

    for mol in test_mols:
        name = mol.GetProp("_Name")
        conf = mol.GetConformer()
        coords = [
            (atom.GetSymbol(), *conf.GetAtomPosition(i))
            for i, atom in enumerate(mol.GetAtoms())
        ]
        match = mol.GetSubstructMatch(OXIME_PAT)
        if match:
            ci, ni, oi = (idx + 1 for idx in match)
            oxime_label = f"[oxime: C{ci}=N{ni}-O{oi}]"
        else:
            oxime_label = "[oxime: not found]"
            print(f"  WARNING: {name} — oxime pattern not found, Stage 3 won't be possible later")

        mol_dir = outdir / name
        mol_dir.mkdir(exist_ok=True)
        (mol_dir / f"{name}_opt.gjf").write_text(_opt_gjf(name, coords, oxime_label))
        (mol_dir / f"{name}_nbo.gjf").write_text(_nbo_gjf(name, oxime_label))
        print(f"  {name:<24} {len(coords):>5}  {oxime_label:>20}   ✓      ✓")

    print(f"\n{len(test_mols)} structures written to {outdir}")
    print("\nSubmit Stage 1 on Citadel via hpc_sync.py:")
    print("  python scripts/dft/hpc_sync.py --mol 002 upload")
    print("  python scripts/dft/hpc_sync.py --mol 002 submit-opt")
    print("  python scripts/dft/hpc_sync.py status")
    print("  python scripts/dft/hpc_sync.py --mol 002 download")
    print("\nOnce Stage 1 shows Normal termination, generate Stage 3 with")
    print("beckmann.dft.inputs.prepare_scan_rigid(), then upload/submit-scan/download it too.")


def main_opt() -> None:
    prepare_opt(
        sdf_path = DATA_OUTPUT / "aimnet_optimized" / "best_per_substrate.sdf",
        outdir   = DATA_OUTPUT / "dft_opt",
    )


def prepare_scan_rigid(mol_dir: Path, name: str, basis: str = BASIS,
                        step: float = 0.05, n_points: int = 6, calcfc: bool = False) -> Path:
    """Stage 3 (rigid-scan architecture) generation, run as a separate step
    AFTER Stage 1 (_opt.gjf) has completed on Citadel and its .log has been
    downloaded to mol_dir -- see prepare_opt()'s docstring for why this can't
    happen upfront like Stages 1-2. Reads the converged geometry from
    {name}_opt.log and writes {name}_scan.gjf via _scan_gjf_rigid().

    step/n_points default to the standard 6-point/0.05 A series (R0..R0+0.30 A)
    -- see _scan_gjf_rigid()'s docstring for the full rationale (mol_006_E's
    missed interior minimum + convergence risk at R0+0.35/+0.4 A).

    calcfc=True applies CalcFC+MaxCycles=300 to all n_points uniformly --
    see _scan_gjf_rigid()'s docstring and JOB_ISSUES.md's 2026-07-16/18 entry."""
    # Local imports: beckmann_nbo.scan imports TEST_IDS/resolve_mol_name from
    # this module at top level, so importing it back at module scope here
    # would be circular (same reason geometry.py was split out).
    from beckmann_nbo.geometry import parse_standard_orientations
    from beckmann_nbo.scan import oxime_atom_map_from_gjf

    ci, ni, oi, oxime_label = oxime_atom_map_from_gjf(mol_dir / f"{name}_opt.gjf")
    lines = (mol_dir / f"{name}_opt.log").read_text().splitlines()
    base_atoms = parse_standard_orientations(lines)[-1][1]

    text = _scan_gjf_rigid(name, base_atoms, ni, oi, oxime_label, basis=basis, step=step, n_points=n_points, calcfc=calcfc)
    out_path = mol_dir / f"{name}_scan.gjf"
    out_path.write_text(text)
    return out_path


# ── single-point NBO workflow ──────────────────────────────────────────────────

def _sp_gjf(name: str, coords: list[tuple], oxime_label: str) -> str:
    """Single-point NBO7 directly on AIMNet2 geometry."""
    header = (
        f"%chk={name}.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{BASIS} sp pop=nbo7read {SOLVENT}\n"
        f"\n"
        f"{name}  {oxime_label}\n"
        f"\n"
        f"{CHARGE} {MULTIPLICITY}\n"
    )
    coord_block = "\n".join(
        f"{sym:<3}  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
        for sym, x, y, z in coords
    )
    return header + coord_block + f"\n$NBO {NBO_KEYWORDS} $END\n\n\n"


def prepare_sp(sdf_path: Path, outdir: Path) -> None:
    """Write single-point NBO .gjf for every molecule in sdf_path."""
    outdir.mkdir(parents=True, exist_ok=True)

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mols  = [m for m in suppl if m is not None]

    print(f"\n{'Name':<24} {'Atoms':>5}  {'Oxime (1-based)':>18}")
    print("-" * 52)

    for mol in mols:
        name = mol.GetProp("_Name")
        conf = mol.GetConformer()
        coords = [
            (atom.GetSymbol(), *conf.GetAtomPosition(i))
            for i, atom in enumerate(mol.GetAtoms())
        ]
        match = mol.GetSubstructMatch(OXIME_PAT)
        if match:
            ci, ni, oi = (idx + 1 for idx in match)
            oxime_label = f"[oxime: C{ci}=N{ni}-O{oi}]"
        else:
            oxime_label = "[oxime: not found]" # TODO: show error message in final UI

        mol_dir = outdir / name
        mol_dir.mkdir(exist_ok=True)
        (mol_dir / f"{name}.gjf").write_text(_sp_gjf(name, coords, oxime_label))
        print(f"  {name:<24} {len(coords):>5}  {oxime_label:>18}")

    print(f"\n{len(mols)} .gjf files → {outdir}")
    print(
        "\nTo submit on Citadel:\n"
        "  python scripts/dft/hpc_sync.py --dir data/output/dft_sp upload\n"
        "  python scripts/dft/hpc_sync.py --dir data/output/dft_sp submit-sp"
    )


def main_sp() -> None:
    prepare_sp(
        sdf_path = DATA_OUTPUT / "aimnet_optimized" / "best_per_substrate.sdf",
        outdir   = DATA_OUTPUT / "dft_sp",
    )
