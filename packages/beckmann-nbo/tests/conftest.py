"""
Shared fixtures for the beckmann-nbo test suite. Each fixture skips
automatically if its prerequisite output does not exist, so the suite can
run at any stage of the (Citadel-dependent) pipeline. Paths come from
beckmann_nbo.config's own DATA_INPUT/DATA_OUTPUT resolution, not a
hardcoded project_root -- this package no longer sits at the repo root.
"""
import json

import pytest

from beckmann_nbo.config import DATA_INPUT, DATA_OUTPUT


@pytest.fixture(scope="session")
def data_input():
    return DATA_INPUT


@pytest.fixture(scope="session")
def data_output():
    return DATA_OUTPUT


@pytest.fixture(scope="session")
def benchmark_meta():
    p = DATA_INPUT / "benchmark_meta.json"
    if not p.exists():
        pytest.skip("benchmark_meta.json not found -- run research/benchmark_pipeline/00_benchmark_to_oximes.py")
    return json.loads(p.read_text())


@pytest.fixture(scope="session")
def best_per_substrate_sdf_path():
    p = DATA_OUTPUT / "aimnet_optimized" / "best_per_substrate.sdf"
    if not p.exists():
        pytest.skip("best_per_substrate.sdf not found -- run research/benchmark_pipeline/02_select_and_optimize.py")
    return p
