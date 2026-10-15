"""`tcren features --chunk` must produce exactly what the whole-set run produces.

Chunking exists for memory: a parsed complex is about a megabyte, so a twenty-thousand-structure
corpus does not fit in one pass. Three separate shell scripts in the benchmark repository grew
their own chunk-split-concatenate loops around this command for that reason, and a loop that
changes the numbers is worse than no loop, so the identity is asserted here rather than assumed.

The property that must NOT be traded for smaller chunks is the batched annotation: one arda call
and one MHC search per chunk, never one per structure.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from tcren.cli import app

_ASSETS = Path(__file__).resolve().parents[1] / "assets"


@pytest.fixture(scope="module")
def structure_list(tmp_path_factory) -> Path:
    pytest.importorskip("arda")
    paths = sorted(_ASSETS.rglob("*.pdb"))[:5]
    if len(paths) < 3:
        pytest.skip("need at least three structure assets to chunk")
    manifest = tmp_path_factory.mktemp("feat") / "structures.txt"
    manifest.write_text("\n".join(str(p) for p in paths) + "\n")
    return manifest


def _run(runner: CliRunner, manifest: Path, out: Path, *extra: str):
    res = runner.invoke(app, ["features", "-s", str(manifest), "-o", str(out),
                              "-i", "placement", "--no-metadata", *extra])
    assert res.exit_code == 0, res.output
    return pl.read_csv(out, separator="\t", infer_schema_length=None).sort("complex.id")


def test_chunked_equals_whole_set(structure_list, tmp_path):
    runner = CliRunner()
    whole = _run(runner, structure_list, tmp_path / "whole.tsv")
    chunked = _run(runner, structure_list, tmp_path / "chunked.tsv", "--chunk", "2")
    assert whole.columns == chunked.columns
    assert whole.equals(chunked), "chunking changed the descriptor table"


def test_resume_is_not_a_supported_scientific_run(tmp_path):
    result = CliRunner().invoke(app, ["features", "-s", str(tmp_path), "--resume"])
    assert result.exit_code == 2
    assert "No such option" in result.output


def test_failed_features_exit_nonzero_and_keep_diagnostics(tmp_path, monkeypatch):
    import tcren.cli
    import tcren.recognition

    monkeypatch.setattr(tcren.cli, "iter_structures", lambda *a, **kw: iter([("failed", object())]))
    monkeypatch.setattr(tcren.recognition, "recognition_table", lambda *a, **kw: [
        {"complex.id": "failed", "error": "missing reference"}])
    out = tmp_path / "failed.tsv"
    result = CliRunner().invoke(app, ["features", "-s", str(tmp_path),
                                      "-o", str(out), "--no-metadata"])
    assert result.exit_code == 1
    assert pl.read_csv(out, separator="\t")["error"][0] == "missing reference"


def test_parts_are_cleaned_up_without_resume(structure_list, tmp_path):
    runner = CliRunner()
    out = tmp_path / "clean.tsv"
    _run(runner, structure_list, out, "--chunk", "2")
    assert not (out.parent / f"{out.stem}.parts").exists(), "part files left behind without --resume"
