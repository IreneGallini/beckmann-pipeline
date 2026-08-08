"""
CLI-layer tests: argparse wiring, the validate_smiles() fast-fail
short-circuit, and pipeline/plot.py's plot_wcnmax(). Fast subset only
exercises the CLI surface (argument parsing, output-path construction,
error handling) -- correctness of the underlying AIMNet2/PySCF pipeline
itself is already covered by test_pipeline.py, not re-tested here except
for the one slow end-to-end CLI test, which exists to catch CLI-layer
wiring bugs (wrong path joins, wrong flag defaults) a pure-pipeline test
can't see.
"""
import argparse
import csv
import sys

import pytest

from beckmann_pyscf import cli
from beckmann_pyscf.pipeline import plot
from test_pipeline import BENCHMARK_CASES


def _fake_args(**kw):
    return argparse.Namespace(**kw)


class TestArgparseSmoke:
    def test_predict_parses_required_flags(self, monkeypatch):
        captured = {}
        monkeypatch.setattr(sys, "argv", ["beckmann-pyscf", "predict", "--smiles", "CCO", "--name", "foo", "--plot"])
        monkeypatch.setattr("beckmann_pyscf.cli.cmd_predict", lambda args: captured.setdefault("args", args))
        cli.main()
        args = captured["args"]
        assert args.smiles == "CCO"
        assert args.name == "foo"
        assert args.plot is True

    def test_predict_neither_smiles_nor_csv_exits_1(self, monkeypatch, capsys):
        """--smiles is no longer required at the argparse level (--csv is a
        valid alternative), so this is now cmd_predict's own validation
        (exit 1), not argparse's (exit 2)."""
        monkeypatch.setattr(sys, "argv", ["beckmann-pyscf", "predict"])
        with pytest.raises(SystemExit) as exc_info:
            cli.main()
        assert exc_info.value.code == 1
        assert "ERROR" in capsys.readouterr().err

    def test_predict_parses_csv_flag(self, monkeypatch):
        captured = {}
        monkeypatch.setattr(sys, "argv", ["beckmann-pyscf", "predict", "--csv", "molecules.csv", "--plot"])
        monkeypatch.setattr("beckmann_pyscf.cli.cmd_predict", lambda args: captured.setdefault("args", args))
        cli.main()
        args = captured["args"]
        assert args.csv == "molecules.csv"
        assert args.smiles is None

    def test_scan_missing_sdf_exits_2(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["beckmann-pyscf", "scan"])
        with pytest.raises(SystemExit) as exc_info:
            cli.main()
        assert exc_info.value.code == 2

    def test_optimize_missing_conformers_sdf_exits_2(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["beckmann-pyscf", "optimize"])
        with pytest.raises(SystemExit) as exc_info:
            cli.main()
        assert exc_info.value.code == 2

    def test_scan_parses_atom_overrides_and_window(self, monkeypatch):
        captured = {}
        monkeypatch.setattr(sys, "argv", [
            "beckmann-pyscf", "scan", "--sdf", "x.sdf",
            "--ci", "1", "--ni", "2", "--oi", "3", "--c-aryl", "4", "--c-alkyl", "5",
            "--r-min", "1.4", "--r-max", "1.9", "--r-step", "0.1",
        ])
        monkeypatch.setattr("beckmann_pyscf.cli.cmd_scan", lambda args: captured.setdefault("args", args))
        cli.main()
        args = captured["args"]
        assert (args.ci, args.ni, args.oi, args.c_aryl, args.c_alkyl) == (1, 2, 3, 4, 5)
        assert (args.r_min, args.r_max, args.r_step) == (1.4, 1.9, 0.1)


class TestValidateSmilesShortCircuit:
    def test_predict_bad_smiles_exits_before_conformers(self, monkeypatch, capsys):
        from beckmann_pyscf import cli_predict

        def boom(*a, **kw):
            raise AssertionError("run_conformers should never be called for a bad SMILES")

        monkeypatch.setattr(cli_predict, "run_conformers", boom)
        args = _fake_args(smiles="not a smiles", name="q", csv=None, out=None, plot=False)
        with pytest.raises(SystemExit) as exc_info:
            cli_predict.cmd_predict(args)
        assert exc_info.value.code == 1
        assert "ERROR" in capsys.readouterr().err

    def test_conformers_non_ketone_smiles_exits_before_run(self, monkeypatch, capsys, tmp_path):
        from beckmann_pyscf import cli_conformers

        def boom(*a, **kw):
            raise AssertionError("run_conformers should never be called for a non-ketone SMILES")

        monkeypatch.setattr(cli_conformers, "run_conformers", boom)
        args = _fake_args(smiles="CCO", name="q", out=str(tmp_path / "out"))
        with pytest.raises(SystemExit) as exc_info:
            cli_conformers.cmd_conformers(args)
        assert exc_info.value.code == 1
        assert "ERROR" in capsys.readouterr().err


class TestCsvBatch:
    def test_missing_csv_file_exits_1(self, monkeypatch, capsys, tmp_path):
        from beckmann_pyscf import cli_predict

        args = _fake_args(smiles=None, name="query", csv=str(tmp_path / "nope.csv"), out=None, plot=False)
        with pytest.raises(SystemExit) as exc_info:
            cli_predict.cmd_predict(args)
        assert exc_info.value.code == 1
        assert "ERROR" in capsys.readouterr().err

    def test_csv_missing_columns_exits_1(self, monkeypatch, capsys, tmp_path):
        from beckmann_pyscf import cli_predict

        csv_path = tmp_path / "bad.csv"
        csv_path.write_text("foo,bar\n1,2\n")
        args = _fake_args(smiles=None, name="query", csv=str(csv_path), out=None, plot=False)
        with pytest.raises(SystemExit) as exc_info:
            cli_predict.cmd_predict(args)
        assert exc_info.value.code == 1
        assert "id" in capsys.readouterr().err

    def test_csv_skips_bad_row_and_continues(self, monkeypatch, tmp_path, capsys):
        from beckmann_pyscf import cli_predict

        csv_path = tmp_path / "mols.csv"
        csv_path.write_text("id,SMILES\nmol_a,BAD\nmol_b,GOOD\n")

        def fake_predict_one(smiles, name, out, plot):
            return smiles == "GOOD"

        monkeypatch.setattr(cli_predict, "_predict_one", fake_predict_one)
        args = _fake_args(smiles=None, name="query", csv=str(csv_path), out=str(tmp_path / "out"), plot=False)
        cli_predict.cmd_predict(args)  # should not raise -- one bad row doesn't abort the batch
        err = capsys.readouterr().err
        assert "SKIP mol_a" in err


class TestPlotWcnmax:
    def _series(self):
        return [
            {"stage": "nbo", "R_NO": 1.50, "weight": 0.45, "MO_index": 39},
            {"stage": "scan_1", "R_NO": 1.55, "weight": 0.42, "MO_index": 39},
            {"stage": "scan_2", "R_NO": 1.60, "weight": 0.41, "MO_index": 39},
        ]

    def test_writes_png(self, tmp_path):
        out_path = tmp_path / "out.png"
        plot.plot_wcnmax(self._series(), None, out_path, name="test")
        assert out_path.exists()
        assert out_path.stat().st_size > 0

    def test_with_minimum(self, tmp_path):
        out_path = tmp_path / "out.png"
        minimum = {"R_star": 1.60, "w_star": 0.41, "MO_index": 39}
        plot.plot_wcnmax(self._series(), minimum, out_path, name="test")
        assert out_path.exists()
        assert out_path.stat().st_size > 0

    def test_empty_series_noop(self, tmp_path):
        out_path = tmp_path / "out.png"
        plot.plot_wcnmax([], None, out_path, name="test")
        assert not out_path.exists()


@pytest.mark.slow
def test_predict_cli_end_to_end(monkeypatch, tmp_path):
    """Exercises CLI-layer wiring (path joins, flag defaults) end-to-end,
    not re-validating AIMNet2/PySCF correctness (test_pipeline.py already
    covers that)."""
    mol_id = sorted(BENCHMARK_CASES)[0]
    smiles, _ = BENCHMARK_CASES[mol_id]
    out_dir = tmp_path / "run"
    monkeypatch.setattr(sys, "argv", [
        "beckmann-pyscf", "predict", "--smiles", smiles, "--name", mol_id,
        "--out", str(out_dir), "--plot",
    ])
    cli.main()

    assert (out_dir / "optimized.sdf").exists()
    assert (out_dir / "wcnmax_series.csv").exists()
    assert (out_dir / "summary.txt").exists()
    assert (out_dir / "wcnmax_vs_rno.png").exists()

    with open(out_dir / "wcnmax_series.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 7
