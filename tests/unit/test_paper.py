"""Unit tests for the natcompsci2022 paper helpers (no network)."""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

from tcren.paper import compare, copy_external_inputs, copy_legacy_results
from tcren.paper.helpers import _read_any


def _write(path: Path, df: pl.DataFrame, sep=","):
    df.write_csv(path, separator=sep)


def test_compare_pass(tmp_path):
    a = pl.DataFrame({"k": [1, 2, 3], "v": [0.0, 1.0, 2.0]})
    b = pl.DataFrame({"k": [1, 2, 3], "v": [0.0, 1.0, 2.0 + 1e-9]})
    _write(tmp_path / "a.csv", a)
    _write(tmp_path / "b.csv", b)
    res = compare(tmp_path / "a.csv", tmp_path / "b.csv", keys=["k"])
    assert res["status"] == "pass"
    assert res["matched"] == 3 and res["only_old"] == 0


def test_compare_detects_value_and_key_diffs(tmp_path):
    _write(tmp_path / "a.csv", pl.DataFrame({"k": [1, 2], "v": [0.0, 1.0]}))
    _write(tmp_path / "b.csv", pl.DataFrame({"k": [1, 2], "v": [0.0, 9.0]}))
    assert compare(tmp_path / "a.csv", tmp_path / "b.csv", keys=["k"])["status"] == "FAIL"
    _write(tmp_path / "c.csv", pl.DataFrame({"k": [1, 3], "v": [0.0, 1.0]}))
    r = compare(tmp_path / "a.csv", tmp_path / "c.csv", keys=["k"])
    assert r["only_old"] == 1 and r["only_new"] == 1 and r["status"] == "FAIL"


def test_read_any_handles_gzip_and_tsv(tmp_path):
    df = pl.DataFrame({"a": [1], "b": ["x"]})
    gz = tmp_path / "t.tsv.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("a\tb\n1\tx\n")
    out = _read_any(gz)
    assert out.columns == ["a", "b"] and out["b"][0] == "x"


def test_copy_external_inputs_gzips(tmp_path):
    # External published inputs (e.g. Birnbaum) go to the shared data dir, gzipped.
    repo_data = tmp_path / "data"
    repo_data.mkdir()
    (repo_data / "Birnbaum.tsv").write_text("peptide\tx\nAAA\t1\n")
    dest = tmp_path / "out"
    n = copy_external_inputs(dest, repo_data=repo_data)
    assert n >= 1
    assert (dest / "Birnbaum.tsv.gz").exists()


def test_copy_legacy_results_routes_to_data_legacy(tmp_path):
    # Legacy baselines (old TCRen matrix, source_data) go to data_legacy/, not the inputs.
    repo_data = tmp_path / "data"
    (repo_data / "source_data").mkdir(parents=True)
    (repo_data / "TCRen_potential.csv").write_text("residue.aa.from,residue.aa.to,TCRen\nA,A,0.1\n")
    (repo_data / "source_data" / "fig1.csv").write_text("a,b\n1,2\n")
    paper_dir = tmp_path / "natcompsci2022"
    copy_legacy_results(paper_dir=paper_dir, repo_data=repo_data)
    assert (paper_dir / "data_legacy" / "TCRen_potential.csv.gz").exists()
    assert (paper_dir / "data_legacy" / "source_data" / "fig1.csv.gz").exists()
    assert _read_any(paper_dir / "data_legacy" / "source_data" / "fig1.csv.gz")["a"][0] == 1
