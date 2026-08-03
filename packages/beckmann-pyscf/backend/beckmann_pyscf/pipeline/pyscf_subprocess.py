"""
Run beckmann_pyscf.engine.pair_nbo.run_from_case() in its own independent OS
process, isolated from any PyTorch/AIMNet2 state in the calling process.

Root cause (confirmed by direct reproduction during a real-substrate smoke
test): PyTorch's bundled OpenMP/MKL runtime (loaded by Auto3D and AIMNet2)
and PySCF's own BLAS/OpenMP runtime conflict when both are resident in the
same process -- PySCF's SCF segfaults on its very first call whenever torch/
AIMNet2 has been used earlier in that process, even with no AIMNet2 model
object still alive at the time. KMP_DUPLICATE_LIB_OK=TRUE (already required
elsewhere in this pipeline to avoid a related libomp duplicate-load abort on
macOS) suppresses the *loud* version of this conflict, not this one -- it
trades a clean abort for unpredictable native memory corruption.

This deliberately uses subprocess.run() to launch a genuinely independent
`python -m ..._pyscf_worker_main` process, NOT multiprocessing.Process
(even with the "spawn" start method). That distinction was verified by direct
reproduction, not assumed: multiprocessing's spawn-context child still
segfaulted PySCF's SCF when launched from a torch-touched parent -- spawn's
child evidently still shares enough OS-level state with the parent (exactly
what isn't clear; both the multiprocessing docs and macOS's own behavior here
are underspecified) that the same conflict survives. A plain subprocess.run()
child, an actually-independent process from the OS's perspective, does not
hit it -- the identical case succeeds and reproduces the expected wCNmax.

wcnmax_pyscf.run_scan_series() interleaves AIMNet2 relaxations and PySCF
single-points 7 times per job -- every PySCF call must go through
run_pyscf_isolated() here rather than calling run_from_case() directly, so
the calling process only ever touches torch/AIMNet2 and each PySCF call runs
in its own clean, genuinely separate process.
"""
import os
import pickle
import subprocess
import sys
import tempfile


def run_pyscf_isolated(case: dict) -> dict:
    """run_from_case(case), executed via `python -m
    beckmann_pyscf.pipeline._pyscf_worker_main` as an independent process
    (subprocess.run, not multiprocessing -- see module docstring for why).
    Raises RuntimeError if the subprocess dies (e.g. segfaults, a negative
    returncode) without producing a result, or if run_from_case() itself
    raised inside the subprocess."""
    label = f"{case.get('name')}/{case.get('stage')}"

    case_fd, case_path = tempfile.mkstemp(suffix=".pkl", prefix="beckmann_pyscf_case_")
    result_path = case_path + ".result"
    try:
        with os.fdopen(case_fd, "wb") as f:
            pickle.dump(case, f)

        env = os.environ.copy()
        # Cap the subprocess's own BLAS/OpenMP thread pool to 1: a real
        # reproduction found PySCF's SCF still segfaulting in a genuinely
        # separate subprocess.run() child, immediately after AIMNet2 use in
        # the parent -- consistent with a thread-oversubscription race
        # (PyTorch's OpenMP pool still resident/idle in the parent competing
        # with PySCF's own multi-threaded BLAS in the child for CPU
        # affinity) rather than a strict same-process conflict. Forcing the
        # child single-threaded is the standard mitigation for exactly this
        # class of cross-library OpenMP conflict.
        for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            env[var] = "1"

        proc = subprocess.run(
            [sys.executable, "-m", "beckmann_pyscf.pipeline._pyscf_worker_main", case_path, result_path],
            capture_output=True, text=True, env=env,
        )
        if proc.returncode != 0:
            crash_note = " (likely a native crash)" if proc.returncode < 0 else ""
            stderr_tail = (proc.stderr or "").strip()[-2000:]
            raise RuntimeError(
                f"PySCF subprocess for {label} exited with code {proc.returncode}"
                f"{crash_note}: {stderr_tail}"
            )
        if not os.path.exists(result_path):
            raise RuntimeError(f"PySCF subprocess for {label} exited 0 but produced no result file")
        with open(result_path, "rb") as f:
            return pickle.load(f)
    finally:
        for p in (case_path, result_path):
            try:
                os.remove(p)
            except OSError:
                pass
