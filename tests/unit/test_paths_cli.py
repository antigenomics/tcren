"""tcren.paths reference resolution + a CLI smoke test (no arda, no network)."""

from __future__ import annotations

import gzip

import pytest

import tcren.paths as paths


def test_reference_local_hit(tmp_path, monkeypatch):
    monkeypatch.setenv("TCREN_DATA_DIR", str(tmp_path))
    nat = tmp_path / "Native2026"
    nat.mkdir()
    with gzip.open(nat / "1ao7.pdb.gz", "wt") as fh:
        fh.write("ATOM\n")
    assert paths.reference_structure_path("1ao7") == nat / "1ao7.pdb.gz"


def test_reference_missing_raises_without_hf(tmp_path, monkeypatch):
    monkeypatch.setenv("TCREN_DATA_DIR", str(tmp_path))
    monkeypatch.setattr(paths, "_fetch_reference_from_hf", lambda pid, folder=paths.NATIVE2026: None)
    with pytest.raises(FileNotFoundError, match="9zzz"):
        paths.reference_structure_path("9zzz")


def test_cli_info_runs():
    from typer.testing import CliRunner

    from tcren.cli import app

    result = CliRunner().invoke(app, ["info"])
    assert result.exit_code == 0
    assert "tcren" in result.stdout
