"""Tests for the beckmann-nbo CLI layer (cli.py + cli_*.py). Mocks the HPC
layer rather than requiring a live cluster: monkeypatches beckmann_nbo.hpc's
subprocess-backed `run()` chokepoint (used internally by cmd_upload/
cmd_submit_opt/cmd_submit_scan/cmd_download/cmd_status), and each cli_*.py
module's own `load_config`/`require_config` bindings (each does
`from beckmann_nbo.hpc import load_config, ...`, which creates a separate
name in that module's namespace, so the patch target is the importing
module, not beckmann_nbo.hpc itself).
"""
import subprocess
import sys

import pytest

import beckmann_nbo.cli as cli
import beckmann_nbo.cli_init as cli_init
import beckmann_nbo.cli_predict as cli_predict
import beckmann_nbo.cli_recover as cli_recover
import beckmann_nbo.cli_verify as cli_verify
import beckmann_nbo.hpc as hpc

FAKE_CONFIG = {
    "HPC_HOST": "user@example.edu",
    "HPC_REMOTE_DIR": "~/beckmann/dft_opt",
    "G16_PATH": "/opt/g16/g16",
    "NBOEXE": "/opt/nbo7/bin/g16nbo.i8.exe",
    "NBO_WRAPPER_DIR": "~/beckmann/nbo7_bin",
}


# ── argparse wiring ─────────────────────────────────────────────────────────

def test_help_lists_all_six_commands(capsys, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["beckmann-nbo", "--help"])
    with pytest.raises(SystemExit):
        cli.main()
    out = capsys.readouterr().out
    for command in ("init", "verify", "predict", "status", "recover", "report"):
        assert command in out


def test_no_command_exits_nonzero(capsys, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["beckmann-nbo"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code != 0


# ── cli_predict: pure/fast helpers (no HPC, no Auto3D/AIMNet2) ─────────────

def test_sanitize_id_prefixes_and_pads():
    assert cli_predict._sanitize_id("42") == "q42"
    assert cli_predict._sanitize_id("5") == "0q5"  # zfill(3)-equivalent to hpc.mol_dirs' own zfill
    assert cli_predict._sanitize_id("a_b") == "qa-b"  # underscores collapsed


def test_smiles_to_oximes_valid():
    # a simple cyclic ketone from the benchmark set
    oximes = cli_predict.smiles_to_oximes("O=C1CCC2=C1C=CC=C2")
    assert len(oximes) >= 1


def test_smiles_to_oximes_invalid_smiles():
    with pytest.raises(ValueError, match="could not parse SMILES"):
        cli_predict.smiles_to_oximes("not a smiles!!")


def test_smiles_to_oximes_no_ketone():
    with pytest.raises(ValueError, match="no ketone group"):
        cli_predict.smiles_to_oximes("CCO")  # ethanol, no ketone


def test_write_smi_uses_mol_prefix_and_id_in_second_token(tmp_path):
    smi_path = cli_predict.write_smi("O=C1CCC2=C1C=CC=C2", "q42", tmp_path)
    assert smi_path.name == "mol_q42.smi"
    lines = smi_path.read_text().splitlines()
    assert lines
    for line in lines:
        smi, name = line.split()
        tokens = name.split("_")
        assert tokens[0] == "mol"
        assert tokens[1] == "q42"
        assert tokens[2] in ("E", "Z", "noEZ")


# ── cli_recover: wrapper wiring ─────────────────────────────────────────────

def test_recover_calls_run_auto_recovery_with_resolved_mols(tmp_path, monkeypatch):
    dft_opt_dir = tmp_path / "dft_opt"
    (dft_opt_dir / "mol_020_E").mkdir(parents=True)
    (dft_opt_dir / "mol_003_E").mkdir(parents=True)

    monkeypatch.setattr(cli_recover, "load_config", lambda: dict(FAKE_CONFIG))
    monkeypatch.setattr(cli_recover, "require_config", lambda config: None)

    captured = {}
    def fake_run_auto_recovery(mols, dft_opt_dir_arg, config, dry_run=False):
        captured["mols"] = sorted(mols)
        captured["dry_run"] = dry_run
    monkeypatch.setattr(cli_recover, "run_auto_recovery", fake_run_auto_recovery)

    class Args:
        dir = str(dft_opt_dir)
        mol = None
        dry_run = True

    cli_recover.cmd_recover(Args())
    assert captured["mols"] == ["mol_003_E", "mol_020_E"]
    assert captured["dry_run"] is True


# ── cli_verify: dry-run and check-failure reporting ─────────────────────────

def test_ssh_check_dry_run_does_not_call_subprocess(monkeypatch):
    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append(1))
    ok = cli_verify._ssh_check("host", "echo ok", dry_run=True)
    assert ok is True
    assert called == []


def test_verify_exits_nonzero_when_checks_fail(monkeypatch):
    monkeypatch.setattr(cli_verify, "load_config", lambda: dict(FAKE_CONFIG))
    monkeypatch.setattr(cli_verify, "require_config", lambda config: None)

    class FakeResult:
        returncode = 1
        stderr = "connection refused"
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: FakeResult())

    class Args:
        dry_run = False

    with pytest.raises(SystemExit) as exc:
        cli_verify.cmd_verify(Args())
    assert exc.value.code == 1


def test_verify_fix_block_mentions_nbo_wrapper_dir():
    block = cli_verify._fix_block(FAKE_CONFIG)
    assert FAKE_CONFIG["NBO_WRAPPER_DIR"] in block
    assert "gaunbo7" in block and "gaunbo6" in block


# ── cli_init: refuses to overwrite without confirmation ─────────────────────

def test_init_aborts_on_no_confirmation(tmp_path, monkeypatch, capsys):
    env_file = tmp_path / ".env"
    env_file.write_text("HPC_HOST=already-here\n")
    monkeypatch.setattr(cli_init, "ENV_FILE", env_file)
    monkeypatch.setattr("builtins.input", lambda *_args: "n")

    class Args:
        pass

    cli_init.cmd_init(Args())
    assert env_file.read_text() == "HPC_HOST=already-here\n"
    assert "Aborted" in capsys.readouterr().out
