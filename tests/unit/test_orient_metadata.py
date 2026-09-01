"""Orient metadata: the format `orient` writes is the format `superimpose` reads, and the shipped
database carries its metadata inside the package (an installed wheel has no repo ``data/``)."""

from __future__ import annotations

import json

import polars as pl
import pytest

from tcren.docking.pipeline import _write_metadata
from tcren.docking.superimpose import _matching_ids, _metadata_path


def test_orient_metadata_json_is_what_superimpose_reads(tmp_path):
    db = tmp_path / "Canonical2026"; db.mkdir()
    _write_metadata(pl.DataFrame([{"pdb.id": "1ao7", "status": "ok",
                                   "mhc.class": "MHCI", "species": "Human"}]),
                    db / "orient_metadata.json")
    assert _matching_ids(db, "MHCI", "Human") == ["1ao7"]


def test_csv_suffix_still_writes_csv(tmp_path):
    path = tmp_path / "orient_metadata.csv"
    _write_metadata(pl.DataFrame([{"pdb.id": "1ao7", "status": "ok"}]), path)
    assert path.read_text().startswith("pdb.id,status")


def test_shipped_database_metadata_rides_in_the_package(tmp_path):
    # No metadata beside the structures — the wheel's copy must answer for Canonical2026.
    db = tmp_path / "Canonical2026"; db.mkdir()
    records = json.loads(_metadata_path(db).read_text())
    assert len(records) == 374
    assert {"pdb.id", "status", "mhc.class", "species"} <= set(records[0])


def test_unknown_database_says_how_to_build_one(tmp_path):
    db = tmp_path / "MyRefs"; db.mkdir()
    with pytest.raises(FileNotFoundError, match="tcren orient"):
        _metadata_path(db)
