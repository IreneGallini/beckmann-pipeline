"""`beckmann-nbo init` -- interactive .env writer.

Prompts for the same 5 keys documented in the repo-root .env.example
(HPC_HOST, HPC_REMOTE_DIR, G16_PATH, NBOEXE, NBO_WRAPPER_DIR) and writes
them in the KEY=value format beckmann_nbo.hpc.load_config() parses. Refuses
to overwrite an existing .env without confirmation.
"""
from beckmann_nbo.hpc import ENV_FILE

PROMPTS = [
    ("HPC_HOST", "SSH destination, e.g. user@cluster.hostname (or an alias from ~/.ssh/config)", True),
    ("HPC_REMOTE_DIR", "Working directory on the cluster, e.g. ~/beckmann/dft_opt", True),
    ("G16_PATH", "Full path to the g16 executable on the cluster, e.g. /opt/g16/g16", True),
    ("NBOEXE", "Full path to the G16/NBO7 interface binary (optional, e.g. /opt/nbo7/bin/g16nbo.i8.exe)", False),
    ("NBO_WRAPPER_DIR", "Directory with an executable gaunbo7/gaunbo6 (optional, e.g. ~/beckmann/nbo7_bin)", False),
]


def cmd_init(args) -> None:
    if ENV_FILE.exists():
        answer = input(f"{ENV_FILE} already exists. Overwrite? [y/N] ").strip().lower()
        if answer != "y":
            print("Aborted -- .env left unchanged.")
            return

    values: dict[str, str] = {}
    print(f"Writing {ENV_FILE} -- press Enter to leave an optional key blank.\n")
    for key, hint, required in PROMPTS:
        suffix = "" if required else " (optional)"
        while True:
            value = input(f"{key}{suffix} -- {hint}\n{key}= ").strip()
            if value or not required:
                break
            print(f"{key} is required.")
        if value:
            values[key] = value

    lines = [f"{key}={value}" for key, value in values.items()]
    ENV_FILE.write_text("\n".join(lines) + "\n")
    print(f"\nWrote {ENV_FILE}.")

    if "NBO_WRAPPER_DIR" not in values:
        print(
            "\nNBO_WRAPPER_DIR not set -- if pop=nbo7read fails on your cluster because\n"
            "the vendor-installed gaunbo7/gaunbo6 aren't executable by your user, copy\n"
            "them somewhere you own and point NBO_WRAPPER_DIR at it, e.g.:\n"
            "  ssh <HPC_HOST> \"mkdir -p ~/beckmann/nbo7_bin && \\\n"
            "    cp <NBOEXE's directory>/gaunbo7 <NBOEXE's directory>/gaunbo6 ~/beckmann/nbo7_bin/ && \\\n"
            "    chmod +x ~/beckmann/nbo7_bin/gaunbo7 ~/beckmann/nbo7_bin/gaunbo6\"\n"
            "then re-run `beckmann-nbo init` and set NBO_WRAPPER_DIR=~/beckmann/nbo7_bin."
        )

    print("\nNext: run `beckmann-nbo verify` to check the cluster is reachable and NBO7 is set up correctly.")
