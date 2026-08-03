"""
CLI entrypoint for pyscf_subprocess.run_pyscf_isolated() -- invoked as
`python -m beckmann_pyscf.pipeline._pyscf_worker_main <case_pkl> <result_pkl>`.
Not meant to be imported or run directly otherwise.

Launched via subprocess.run(), not multiprocessing.Process: a real
reproduction found that multiprocessing's spawn-context child (which shares
some OS-level state with the parent that plain subprocess.run's fresh
interpreter does not) still segfaults PySCF's SCF when the parent process
has touched PyTorch/AIMNet2, while an actually-independent process launch
does not. See pyscf_subprocess.py's module docstring for the full story.
"""
import pickle
import sys


def main() -> None:
    case_path, result_path = sys.argv[1], sys.argv[2]
    with open(case_path, "rb") as f:
        case = pickle.load(f)

    from beckmann_pyscf.engine.pair_nbo import run_from_case

    result = run_from_case(case)

    with open(result_path, "wb") as f:
        pickle.dump(result, f)


if __name__ == "__main__":
    main()
