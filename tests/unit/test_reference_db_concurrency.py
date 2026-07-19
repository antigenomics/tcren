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


def test_reference_db_serializes_concurrent_builds(tmp_path, monkeypatch):
    import arda.mmseqs as mmseqs

    db = tmp_path / "alleles_db"
    calls = {"createdb": 0, "createindex": 0}
    counter_lock = threading.Lock()

    def fake_createdb(fasta, out, dbtype=1):
        with counter_lock:
            calls["createdb"] += 1
        time.sleep(0.2)  # widen the race window so an unguarded build would double-run
        Path(str(out) + ".dbtype").write_text("db")  # createdb marker

    def fake_run(args, check=True):
        with counter_lock:
            calls["createindex"] += 1
        Path(str(args[1]) + ".idx.dbtype").write_text("idx")  # createindex marker (written last)

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
        Path(str(out) + ".dbtype").write_text("db")

    def fake_run(args, check=True):
        Path(str(args[1]) + ".idx.dbtype").write_text("idx")

    monkeypatch.setattr(mmseqs, "createdb", fake_createdb)
    monkeypatch.setattr(mmseqs, "run", fake_run)

    reference.reference_db(cache_dir=tmp_path)
    reference.reference_db(cache_dir=tmp_path)  # fresh cache → no rebuild
    assert calls["n"] == 1
