"""Unit tests for the native database accessor and bootstrap helpers (no network)."""

from __future__ import annotations

import json

import polars as pl
import pytest

from tcren.native import NativeDatabase, default_native_root
from tcren.native.database import CIF_SUFFIX


def _make_fake_db(root):
    (root / "cif").mkdir(parents=True)
    for pid in ("1ao7", "5xyz"):
        (root / "cif" / f"{pid}{CIF_SUFFIX}").write_text("# fake cif\n")
    pl.DataFrame(
        {"pdb_id": ["1ao7"], "tcr_type": ["Alpha"], "cdr3_sequences": ["CAVF"]}
    ).write_csv(root / "tcr_chain_data.tsv", separator="\t")
    pl.DataFrame(
        {"PDB_ID": ["1ao7"], "TCR_complex": ["CLASSI"], "Epitope": ["LLF"]}
    ).write_csv(root / "tcr_complexes_data.tsv", separator="\t")
    return NativeDatabase(root)


def test_custom_root_via_argument(tmp_path):
    db = _make_fake_db(tmp_path / "mydb")
    assert db.root == tmp_path / "mydb"
    assert db.is_present()
    assert sorted(db.pdb_ids()) == ["1ao7", "5xyz"]
    assert db.cif_for("1ao7").name == f"1ao7{CIF_SUFFIX}"


def test_default_root_env(monkeypatch, tmp_path):
    monkeypatch.setenv("TCREN_NATIVE_DIR", str(tmp_path / "envdb"))
    assert default_native_root() == tmp_path / "envdb"
    monkeypatch.delenv("TCREN_NATIVE_DIR")
    assert default_native_root().name == "native"


def test_missing_cif_raises(tmp_path):
    db = _make_fake_db(tmp_path / "db")
    with pytest.raises(FileNotFoundError):
        db.cif_for("9zzz")


def test_tables_load(tmp_path):
    db = _make_fake_db(tmp_path / "db")
    assert db.complex_data["Epitope"][0] == "LLF"
    assert db.chain_data["tcr_type"][0] == "Alpha"


def test_version_roundtrip(tmp_path):
    db = _make_fake_db(tmp_path / "db")
    assert db.version() == {}  # not yet written
    db.version_path.write_text(json.dumps({"n_cif": 2}))
    assert db.version()["n_cif"] == 2


def test_not_present_when_tables_missing(tmp_path):
    root = tmp_path / "empty"
    (root / "cif").mkdir(parents=True)
    assert not NativeDatabase(root).is_present()
