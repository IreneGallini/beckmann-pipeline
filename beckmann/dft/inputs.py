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

from beckmann.config import (
    DATA_OUTPUT,
    FUNCTIONAL, BASIS, NPROC, MEM_GB, CHARGE, MULTIPLICITY,
    NBO_KEYWORDS, SOLVENT,
)
from beckmann.dft.geometry import displace_leaving_group

TEST_IDS  = {"002", "006", "020", "021", "014", "029"}
OXIME_PAT = Chem.MolFromSmarts('[C:1]=[N:2]-[O+:3]')


def resolve_mol_name(mol_id: str, dft_opt_dir: Path) -> str | None:
    """Find the isomer-suffixed mol dir name for a numeric id (e.g. '014' ->
    'mol_014_Z') -- the AIMNet2-lower-energy isomer isn't always E."""
    matches = sorted(dft_opt_dir.glob(f"mol_{mol_id.zfill(3)}_*"))
    return matches[0].name if matches else None


# ── three-stage opt workflow ────────────────────────────────────────────────────

def _opt_gjf(name: str, coords: list[tuple], oxime_label: str, basis: str = BASIS) -> str:
    """Stage 1: geometry optimisation no NBO block."""
    return (
        f"%chk={name}_opt.chk\n"
        f"%nprocshared={NPROC}\n"
        f"%mem={MEM_GB}GB\n"
        f"#p {FUNCTIONAL}/{basis} opt {SOLVENT}\n"
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
    basis: str = BASIS, step: float = 0.1, n_points: int = 4,
) -> str:
    """Stage 3 (rigid-scan architecture): n_points independent points
    (R0+step .. R0+step*n_points), each built from the SAME Stage-1 base
    geometry -- not chained point-to-point, not extracted from an internal
    Gaussian scan walk. Per point: a rigid O-N displacement (pure geometry,
    via displace_leaving_group), then a constrained optimization
    (opt=ModRedundant, bond frozen at the new length, everything else
    relaxes), then a same-checkpoint Link1 continuation into a full
    NBO7/CMO single point. Matches the supervisor's reference file
    (oxime_001_scan.gjf) -- see JOB_ISSUES.md for the full rationale.

    Defaults (step=0.1, n_points=4) match the standard 5-point series (with
    R0 covered separately by Stage 2). A finer step/more points covers the
    same R0..R0+0.4 A range at higher resolution -- e.g. mol_006_E's
    supervisor-requested 0.05 A / 8-point rerun, to check whether a narrow
    interior wCNmax minimum was being stepped over by the coarser default
    grid (see JOB_ISSUES.md).

    NBO keywords deliberately omit her NBOMO=P120 print-window restriction:
    parse_cmo.py was fixed earlier to search the entire virtual manifold
    unconditionally (see its docstring) after discovering real target
    antibonds mixing in above a narrower window for some substrates -- a
    fixed narrow NBOMO range would reintroduce that at the Gaussian-printing
    level, upstream of anything the parser can recover.
    """
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
            f"#p {FUNCTIONAL}/{basis} opt=(ModRedundant) SCF=(Tight,XQC) NoSymm {SOLVENT}\n"
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
                        step: float = 0.1, n_points: int = 4) -> Path:
    """Stage 3 (rigid-scan architecture) generation, run as a separate step
    AFTER Stage 1 (_opt.gjf) has completed on Citadel and its .log has been
    downloaded to mol_dir -- see prepare_opt()'s docstring for why this can't
    happen upfront like Stages 1-2. Reads the converged geometry from
    {name}_opt.log and writes {name}_scan.gjf via _scan_gjf_rigid().

    step/n_points default to the standard 4-point/0.1 A series; pass a finer
    step for a higher-resolution scan (e.g. step=0.05, n_points=8, as used
    for mol_006_E after its standard-resolution scan missed a real interior
    wCNmax minimum -- see Notes.md)."""
    # Local imports: beckmann.dft.scan imports TEST_IDS/resolve_mol_name from
    # this module at top level, so importing it back at module scope here
    # would be circular (same reason geometry.py was split out).
    from beckmann.dft.geometry import parse_standard_orientations
    from beckmann.dft.scan import oxime_atom_map_from_gjf

    ci, ni, oi, oxime_label = oxime_atom_map_from_gjf(mol_dir / f"{name}_opt.gjf")
    lines = (mol_dir / f"{name}_opt.log").read_text().splitlines()
    base_atoms = parse_standard_orientations(lines)[-1][1]

    text = _scan_gjf_rigid(name, base_atoms, ni, oi, oxime_label, basis=basis, step=step, n_points=n_points)
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
