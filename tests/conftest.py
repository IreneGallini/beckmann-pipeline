"""
Shared fixtures for the Beckmann pipeline test suite.

Each fixture skips automatically if its prerequisite output does not exist,
so the full test suite can be run at any stage of the pipeline.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import json
import pytest
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent


@pytest.fixture(scope="session")
def project_root():
    return PROJECT_ROOT


@pytest.fixture(scope="session")
def molecules_smi_path(project_root):
    p = project_root / "data" / "input" / "molecules.smi"
    if not p.exists():
        pytest.skip("molecules.smi not found — run scripts/00_benchmark_to_oximes.py")
    return p


@pytest.fixture(scope="session")
def benchmark_meta(project_root):
    p = project_root / "data" / "input" / "benchmark_meta.json"
    if not p.exists():
        pytest.skip("benchmark_meta.json not found — run scripts/00_benchmark_to_oximes.py")
    return json.loads(p.read_text())


@pytest.fixture(scope="session")
def conformers_sdf_path(project_root):
    paths = sorted(
        (project_root / "data" / "output" / "conformers").glob("molecules_*/molecules_out.sdf")
    )
    if not paths:
        pytest.skip("No conformers SDF found — run scripts/01_smiles_to_conformers.py")
    return paths[-1]


@pytest.fixture(scope="session")
def aimnet_sdf_path(project_root):
    p = project_root / "data" / "output" / "aimnet_optimized" / "best_aimnet_optimized.sdf"
    if not p.exists():
        pytest.skip("AIMNet SDF not found — run scripts/02_select_and_optimize.py")
    return p


@pytest.fixture(scope="session")
def best_per_substrate_sdf_path(project_root):
    p = project_root / "data" / "output" / "aimnet_optimized" / "best_per_substrate.sdf"
    if not p.exists():
        pytest.skip("best_per_substrate.sdf not found — run scripts/02_select_and_optimize.py")
    return p


@pytest.fixture(scope="session")
def benchmark_csv_path(project_root):
    p = project_root / "data" / "output" / "week1_benchmark_results.csv"
    if not p.exists():
        pytest.skip("Benchmark CSV not found — run scripts/04_extract_dihedrals_and_predict.py")
    return p
