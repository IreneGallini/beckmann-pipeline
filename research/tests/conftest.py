"""
Shared fixtures for research/'s test suite -- these validate the benchmark
batch pipeline's generated artifacts (molecules.smi, benchmark_meta.json,
conformer/AIMNet2 SDFs, TS search outputs), not beckmann-core's pure
functions in isolation (see packages/beckmann-core/tests/ for those). Each
fixture skips automatically if its prerequisite output does not exist.
"""
import json
import sys
from pathlib import Path

import pytest

_RESEARCH_DIR = Path(__file__).resolve().parent.parent
if str(_RESEARCH_DIR) not in sys.path:
    sys.path.insert(0, str(_RESEARCH_DIR))

from beckmann_nbo.config import DATA_INPUT, DATA_OUTPUT, ROOT


@pytest.fixture(scope="session")
def project_root():
    return ROOT


@pytest.fixture(scope="session")
def molecules_smi_path():
    p = DATA_INPUT / "molecules.smi"
    if not p.exists():
        pytest.skip("molecules.smi not found -- run research/benchmark_pipeline/00_benchmark_to_oximes.py")
    return p


@pytest.fixture(scope="session")
def benchmark_meta():
    p = DATA_INPUT / "benchmark_meta.json"
    if not p.exists():
        pytest.skip("benchmark_meta.json not found -- run research/benchmark_pipeline/00_benchmark_to_oximes.py")
    return json.loads(p.read_text())


@pytest.fixture(scope="session")
def conformers_sdf_path():
    paths = sorted((DATA_OUTPUT / "conformers").glob("molecules_*/molecules_out.sdf"))
    if not paths:
        pytest.skip("No conformers SDF found -- run research/benchmark_pipeline/01_smiles_to_conformers.py")
    return paths[-1]


@pytest.fixture(scope="session")
def aimnet_sdf_path():
    p = DATA_OUTPUT / "aimnet_optimized" / "best_aimnet_optimized.sdf"
    if not p.exists():
        pytest.skip("AIMNet SDF not found -- run research/benchmark_pipeline/02_select_and_optimize.py")
    return p


@pytest.fixture(scope="session")
def best_per_substrate_sdf_path():
    p = DATA_OUTPUT / "aimnet_optimized" / "best_per_substrate.sdf"
    if not p.exists():
        pytest.skip("best_per_substrate.sdf not found -- run research/benchmark_pipeline/02_select_and_optimize.py")
    return p
