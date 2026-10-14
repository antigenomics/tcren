"""Concurrency guard for the on-demand MHC mmseqs index build.

tcren is routinely fanned out one process/thread per sample (SLURM array / Nextflow) against a
shared ``data/mhc_cache``. ``reference.reference_db`` must serialize the ``createdb`` +
``createindex`` build (via ``arda._locking.build_lock``) so concurrent callers never build twice
or search a half-written index. mmseqs itself is mocked here so the test is fast and deterministic;
it fails if the locking is dropped.
"""

from __future__ import annotations

import threading
import time
from pathlib import Path

import pytest

pytest.importorskip("arda.mmseqs")  # reference_db builds via arda's mmseqs + lock; skip if absent (CI)

from tcren.mhc import reference


@pytest.fixture(autouse=True)
def local_reference(tmp_path, monkeypatch):
    fasta = tmp_path / "reference.fasta"
    fasta.write_text(">synthetic\nACDEFGHIK\n")
    monkeypatch.setattr(reference, "reference_fasta", lambda: fasta)


def _write_db(prefix):
    Path(str(prefix)).write_bytes(b"seq\0")
    Path(str(prefix) + ".index").write_text("0\t0\t4\n")
    Path(str(prefix) + ".dbtype").write_bytes(b"\x01\0\0\0")


def test_reference_db_serializes_concurrent_builds(tmp_path, monkeypatch):
    import arda.mmseqs as mmseqs

    db = tmp_path / "alleles_db"
    calls = {"createdb": 0, "createindex": 0}
    counter_lock = threading.Lock()

    def fake_createdb(fasta, out, dbtype=1):
        with counter_lock:
            calls["createdb"] += 1
        time.sleep(0.2)  # widen the race window so an unguarded build would double-run
        _write_db(out)
        _write_db(str(out) + "_h")

    def fake_run(args, check=True):
        with counter_lock:
            calls["createindex"] += 1
        _write_db(str(args[1]) + ".idx")

    monkeypatch.setattr(mmseqs, "createdb", fake_createdb)
    monkeypatch.setattr(mmseqs, "run", fake_run)

    results: list[Path] = []
    res_lock = threading.Lock()

    def worker():
        r = reference.reference_db(cache_dir=tmp_path)
        with res_lock:
            results.append(r)

    threads = [threading.Thread(target=worker) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    # The build ran exactly once despite 8 concurrent callers, and all got the same db path.
    assert calls["createdb"] == 1, f"build was not serialized: createdb ran {calls['createdb']}x"
    assert calls["createindex"] == 1
    assert {str(r) for r in results} == {str(db)}


def test_reference_db_second_call_is_a_cache_hit(tmp_path, monkeypatch):
    import arda.mmseqs as mmseqs

    calls = {"n": 0}

    def fake_createdb(fasta, out, dbtype=1):
        calls["n"] += 1
        _write_db(out)
        _write_db(str(out) + "_h")

    def fake_run(args, check=True):
        _write_db(str(args[1]) + ".idx")

    monkeypatch.setattr(mmseqs, "createdb", fake_createdb)
    monkeypatch.setattr(mmseqs, "run", fake_run)

    reference.reference_db(cache_dir=tmp_path)
    reference.reference_db(cache_dir=tmp_path)  # fresh cache → no rebuild
    assert calls["n"] == 1


@pytest.mark.parametrize("broken", ["missing", "truncated"])
def test_reference_db_repairs_incomplete_index_despite_fresh_markers(tmp_path, monkeypatch, broken):
    import arda.mmseqs as mmseqs

    db = tmp_path / "alleles_db"
    for suffix in ("", "_h", ".idx"):
        _write_db(str(db) + suffix)
    if broken == "missing":
        Path(str(db) + ".idx").unlink()
    else:
        Path(str(db) + ".idx").write_bytes(b"x")
    calls = []

    def createdb(fasta, out, dbtype):
        calls.append(out)
        _write_db(out)
        _write_db(str(out) + "_h")

    monkeypatch.setattr(mmseqs, "createdb", createdb)
    monkeypatch.setattr(mmseqs, "run", lambda args: _write_db(str(args[1]) + ".idx"))
    assert reference.reference_db(cache_dir=tmp_path) == db
    assert len(calls) == 1
    assert Path(str(db) + ".idx").read_bytes() == b"seq\0"
    backup, = tmp_path.glob(".alleles_db-backup-*")
    assert (backup / "alleles_db.idx.index").read_text() == "0\t0\t4\n"
    if broken == "truncated":
        assert (backup / "alleles_db.idx").read_bytes() == b"x"


def test_reference_db_does_not_publish_an_incomplete_build(tmp_path, monkeypatch):
    import arda.mmseqs as mmseqs

    def createdb(fasta, out, dbtype):
        _write_db(out)
        _write_db(str(out) + "_h")

    monkeypatch.setattr(mmseqs, "createdb", createdb)
    monkeypatch.setattr(mmseqs, "run", lambda args: None)
    with pytest.raises(RuntimeError, match="incomplete"):
        reference.reference_db(cache_dir=tmp_path)
    assert not (tmp_path / "alleles_db.dbtype").exists()
