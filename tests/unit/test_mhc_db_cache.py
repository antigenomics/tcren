"""Compiled+indexed mmseqs reference DB cache for MHC mapping."""

from __future__ import annotations

import pytest

pytest.importorskip("arda")

from tcren.mhc import reference  # noqa: E402


def _have_reference() -> bool:
    try:
        reference.reference_fasta()
        return True
    except FileNotFoundError:
        return False


@pytest.mark.skipif(not _have_reference(), reason="MHC reference not built (tcren build-mhc-ref)")
def test_reference_db_builds_and_indexes(tmp_path):
    db = reference.reference_db(cache_dir=tmp_path)
    assert db == tmp_path / "alleles_db"
    assert db.with_name("alleles_db.dbtype").exists()        # createdb
    assert db.with_name("alleles_db.idx.dbtype").exists()    # createindex (the speedup)


@pytest.mark.skipif(not _have_reference(), reason="MHC reference not built")
def test_reference_db_is_cached(tmp_path):
    # Second call must not rebuild (marker mtime unchanged).
    db = reference.reference_db(cache_dir=tmp_path)
    m0 = db.with_name("alleles_db.idx.dbtype").stat().st_mtime
    reference.reference_db(cache_dir=tmp_path)
    assert db.with_name("alleles_db.idx.dbtype").stat().st_mtime == m0
