"""
Repo-level environment sanity check: everything this pipeline needs is
importable, including pyscf and flask -- both confirmed absent from
beckmann-pipeline's own environment.yml/pyproject.toml (undeclared,
"happens to be installed" dependencies there), so this repo's own
environment.yml must declare them explicitly and this test is what actually
verifies that. Mirrors beckmann-pipeline/tests/test_env.py's pattern.
"""
import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")  # must be set before Auto3D/AIMNet2 import (macOS)


def test_core_imports():
    import ase  # noqa: F401
    import matplotlib  # noqa: F401
    import numpy  # noqa: F401
    import rdkit  # noqa: F401


def test_ml_potential_imports():
    import aimnet  # noqa: F401
    import Auto3D  # noqa: F401


def test_pyscf_imports():
    import pyscf  # noqa: F401


def test_flask_imports():
    import flask  # noqa: F401
