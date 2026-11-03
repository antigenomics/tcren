"""Small compiled MMseqs sequence databases, without dense prefilter indexes."""

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
def test_reference_db_builds_without_dense_prefilter_index(tmp_path):
    db = reference.reference_db(cache_dir=tmp_path)
    assert db.name == "alleles_db"
    assert db.parent.parent == tmp_path / "generations"
    assert db.with_name("alleles_db.dbtype").exists()        # createdb
    assert db.with_name("alleles_db_h.dbtype").exists()
    assert not list(db.parent.glob("*.idx*"))


@pytest.mark.skipif(not _have_reference(), reason="MHC reference not built")
def test_reference_db_is_cached(tmp_path):
    # Second call must not rebuild (marker mtime unchanged).
    db = reference.reference_db(cache_dir=tmp_path)
    m0 = db.with_name("alleles_db.dbtype").stat().st_mtime
    reference.reference_db(cache_dir=tmp_path)
    assert db.with_name("alleles_db.dbtype").stat().st_mtime == m0


def test_bundled_panel_distinguishes_b2702_and_b2705(tmp_path, monkeypatch):
    import arda.mmseqs as mmseqs
    from tcren.mhc.imgt import _read_fasta
    from tcren.mhc.mapper import _best_hits

    fasta = reference._BUNDLED_DIR / 'alleles.aa.fasta'
    sequences = {header.split('|')[0]: sequence for header, sequence in _read_fasta(fasta)}
    names = ['HLA-B*27:02', 'HLA-B*27:05']
    assert sequences[names[0]] != sequences[names[1]]
    monkeypatch.setattr(reference, 'reference_fasta', lambda: fasta)
    db = reference.reference_db(cache_dir=tmp_path / 'db')
    query, output = tmp_path / 'queries.fasta', tmp_path / 'hits.tsv'
    query.write_text(''.join(f'>{i}\n{sequences[name]}\n' for i, name in enumerate(names)))
    mmseqs.easy_search(query, db, output, tmp_path / 'work', search_type=1, max_seqs=50, threads=1)
    best = _best_hits(output)
    assert [best[str(i)]['target'].split('|')[0] for i in range(2)] == names
    assert all(best[str(i)]['pident'] == 100 for i in range(2))
    assert not list(db.parent.glob('*.idx*'))
    assert sum(p.stat().st_size for p in db.parent.iterdir() if p.is_file()) < 1_000_000
