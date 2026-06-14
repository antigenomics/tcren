"""Unit tests for the natcompsci2022 paper helpers (no network)."""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

from tcren.paper import compare, copy_paper_data
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


def test_copy_paper_data_gzips(tmp_path):
    repo_data = tmp_path / "data"
    (repo_data).mkdir()
    (repo_data / "Birnbaum.tsv").write_text("peptide\tx\nAAA\t1\n")
    (repo_data / "source_data").mkdir()
    (repo_data / "source_data" / "fig1.csv").write_text("a,b\n1,2\n")
    dest = tmp_path / "out"
    n = copy_paper_data(dest, repo_data=repo_data)
    assert n >= 2
    assert (dest / "Birnbaum.tsv.gz").exists()
    assert (dest / "source_data" / "fig1.csv.gz").exists()
    # round-trips
    assert _read_any(dest / "source_data" / "fig1.csv.gz")["a"][0] == 1
