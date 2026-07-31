"""
Verifies the "zero Citadel/Gaussian dependency" hard rule as an actual test,
not just a stated intention:

  1. No file under backend/ references beckmann_nbo (the package that owns
     Citadel/Gaussian code -- hpc.py, inputs.py, etc.), paramiko, or fabric
     -- a static source check.
  2. beckmann_nbo is never loaded into sys.modules after importing the full
     pipeline package.
  3. A full pipeline.predict() run makes no outbound (non-loopback) network
     connection -- if an SSH workflow (or anything else) were ever reached,
     this fails loudly with a clear exception instead of silently
     succeeding or hanging on an unreachable host.

Note: earlier versions of check 2 also asserted paramiko/fabric were absent
from sys.modules globally -- dropped, because that's a property of the
whole Python process, not of this package's own import chain. pysisyphus
(used only by research/ts_ml/, nothing to do with this package) imports
paramiko itself for unrelated reasons, which made that check fail whenever
this test ran in the same pytest session as research/'s tests even though
beckmann-pyscf itself never touched it. "beckmann_nbo was never imported"
is the precise, correct invariant for this repo's package structure: if
beckmann-pyscf never imports beckmann_nbo, it structurally cannot reach any
Citadel code, regardless of what unrelated libraries elsewhere in the
process happen to use.
"""
import socket
from pathlib import Path

import pytest

from beckmann_pyscf import pipeline

_LOOPBACK = {"127.0.0.1", "::1", "localhost"}
_BACKEND_DIR = Path(__file__).resolve().parent.parent


class NetworkAttempted(Exception):
    pass


def test_no_hpc_reference_in_source():
    """Static check: no .py file under backend/ mentions beckmann_nbo,
    paramiko, or fabric."""
    offenders = []
    for py_file in _BACKEND_DIR.rglob("*.py"):
        if py_file == Path(__file__):
            continue
        text = py_file.read_text()
        if "beckmann_nbo" in text or "import paramiko" in text or "import fabric" in text:
            offenders.append(str(py_file.relative_to(_BACKEND_DIR)))
    assert not offenders, f"HPC-related reference found in: {offenders}"


def test_beckmann_nbo_never_imported():
    """After importing the full pipeline package, beckmann_nbo (which owns
    hpc.py) must never have been loaded."""
    import sys
    top_level_names = {m.split(".")[0] for m in sys.modules}
    assert "beckmann_nbo" not in top_level_names


@pytest.mark.slow
def test_predict_makes_no_outbound_network_calls(monkeypatch, tmp_path):
    """Runs a full predict() call with socket.socket.connect patched to
    raise on any non-loopback address. Loopback connections pass through to
    the real implementation unchanged (some libraries' internal
    multiprocessing/IPC legitimately uses localhost sockets; that's not
    what this test is checking for) -- only a genuine outbound attempt
    (e.g. to citadel.chem.cmu.edu) trips it."""
    original_connect = socket.socket.connect

    def guarded_connect(self, address, *a, **kw):
        host = address[0] if isinstance(address, tuple) else address
        if host not in _LOOPBACK:
            raise NetworkAttempted(f"outbound connection attempted to {address!r}")
        return original_connect(self, address, *a, **kw)

    monkeypatch.setattr(socket.socket, "connect", guarded_connect)

    result = pipeline.predict("O=C1CCC2=CC=CC=C21", workdir=tmp_path)  # alpha-tetralone
    assert result["prediction"] in ("R", "F")
    assert len(result["wcnmax_series"]) >= 3
