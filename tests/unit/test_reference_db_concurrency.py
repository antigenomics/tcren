"""Concurrency guard for the on-demand MHC mmseqs index build.

tcren is routinely fanned out one process/thread per sample (SLURM array / Nextflow) against a
shared ``data/mhc_cache``. ``reference.reference_db`` must serialize the ``createdb``
publication (via ``arda._locking.build_lock``) so concurrent callers never build twice
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
    monkeypatch.setattr("arda.mmseqs.version", lambda: "test-version")


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
    assert calls["createindex"] == 0
    assert len(results) == 8
    assert len(set(results)) == 1
    assert results[0].name == db.name
    assert results[0].parent.parent == tmp_path / "generations"


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
        Path(str(db) + "_h").unlink()
    else:
        Path(str(db) + "_h").write_bytes(b"x")
    calls = []

    def createdb(fasta, out, dbtype):
        calls.append(out)
        _write_db(out)
        _write_db(str(out) + "_h")

    monkeypatch.setattr(mmseqs, "createdb", createdb)
    monkeypatch.setattr(mmseqs, "run", lambda args: _write_db(str(args[1]) + ".idx"))
    new = reference.reference_db(cache_dir=tmp_path)
    assert new != db
    assert len(calls) == 1
    assert Path(str(new) + "_h").read_bytes() == b"seq\0"
    assert not Path(str(new) + ".idx").exists()
    assert (tmp_path / "alleles_db.idx.index").read_text() == "0\t0\t4\n"
    if broken == "truncated":
        assert Path(str(db) + "_h").read_bytes() == b"x"


def test_reference_db_does_not_publish_an_incomplete_build(tmp_path, monkeypatch):
    import arda.mmseqs as mmseqs

    def createdb(fasta, out, dbtype):
        _write_db(out)

    monkeypatch.setattr(mmseqs, "createdb", createdb)
    monkeypatch.setattr(mmseqs, "run", lambda args: None)
    with pytest.raises(RuntimeError, match="incomplete"):
        reference.reference_db(cache_dir=tmp_path)
    assert not (tmp_path / "alleles_db.dbtype").exists()


def test_rebuild_preserves_old_readers_and_detects_same_size_corruption(tmp_path, monkeypatch):
    import os
    import arda.mmseqs as mmseqs

    def createdb(fasta, out, dbtype):
        _write_db(out)
        _write_db(str(out) + "_h")

    monkeypatch.setattr(mmseqs, "createdb", createdb)
    monkeypatch.setattr(mmseqs, "run", lambda args: _write_db(str(args[1]) + ".idx"))
    old = reference.reference_db(cache_dir=tmp_path)
    original = {p.name: p.read_bytes() for p in old.parent.iterdir() if p.is_file()}
    fasta = reference.reference_fasta()
    stat = fasta.stat()
    fasta.write_text(">synthetic\nKKKKKKKKK\n")
    os.utime(fasta, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    new = reference.reference_db(cache_dir=tmp_path)
    assert new != old
    assert {p.name: p.read_bytes() for p in old.parent.iterdir() if p.is_file()} == original
    new.write_bytes(b"bad\0")  # same size and valid offsets, but wrong contents
    repaired = reference.reference_db(cache_dir=tmp_path)
    assert repaired not in (old, new)
    assert repaired.read_bytes() == b"seq\0"
    assert new.read_bytes() == b"bad\0"


def test_mmseqs_version_change_selects_new_generation(tmp_path, monkeypatch):
    import arda.mmseqs as mmseqs

    def createdb(fasta, out, dbtype):
        _write_db(out)
        _write_db(str(out) + "_h")

    monkeypatch.setattr(mmseqs, "createdb", createdb)
    monkeypatch.setattr(mmseqs, "run", lambda args: _write_db(str(args[1]) + ".idx"))
    old = reference.reference_db(cache_dir=tmp_path)
    monkeypatch.setattr(mmseqs, "version", lambda: "another-version")
    assert reference.reference_db(cache_dir=tmp_path) != old


def test_separate_processes_share_one_complete_generation(tmp_path):
    import subprocess
    import sys

    driver = tmp_path / 'build_probe.py'
    driver.write_text('''import sys,time
from pathlib import Path
import arda.mmseqs as mmseqs
from tcren.mhc import reference
root=Path(sys.argv[1])
reference.reference_fasta=lambda: root/'reference.fasta'
mmseqs.version=lambda: 'synthetic-version'
def write(prefix):
    Path(str(prefix)).write_bytes(b'seq\\0')
    Path(str(prefix)+'.index').write_text('0\\t0\\t4\\n')
    Path(str(prefix)+'.dbtype').write_bytes(b'\\x01\\0\\0\\0')
def create(fasta,out,dbtype):
    with (root/'builds').open('a') as log: log.write('build\\n')
    time.sleep(.1)
    write(out);write(str(out)+'_h')
mmseqs.createdb=create
mmseqs.run=lambda args: write(str(args[1])+'.idx')
print(reference.reference_db(root))
''')
    workers = [subprocess.Popen([sys.executable, str(driver), str(tmp_path)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
               for _ in range(3)]
    try:
        outputs = [worker.communicate(timeout=20) for worker in workers]
        assert all(worker.returncode == 0 for worker in workers), outputs
        assert len({stdout.strip() for stdout, _ in outputs}) == 1
        assert (tmp_path / 'builds').read_text() == 'build\n'
    finally:
        for worker in workers:
            if worker.poll() is None:
                worker.kill()
            worker.wait()
