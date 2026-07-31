"""
Run Auto3D conformer generation on the full benchmark's molecules.smi.
Moved here (from beckmann/conformers.py's main()) since it hardcodes the
benchmark's own paths -- generate_conformers() itself is beckmann_core's
reusable function, called unmodified.
"""
from beckmann_core.conformers import generate_conformers
from beckmann_nbo.config import DATA_INPUT, DATA_OUTPUT


def main() -> None:
    input_src  = DATA_INPUT  / "molecules.smi"
    output_dir = DATA_OUTPUT / "conformers"

    sdf_path = generate_conformers(input_src, output_dir)
    print(f"\nDone. Output written to: {sdf_path}")


if __name__ == '__main__':
    main()
