"""Batch dispatch regressions without sequence-search subprocesses."""
from pathlib import Path


def test_archive_is_annotated_in_one_batch(monkeypatch):
    from tcren import annotation, structure
    from tcren.annotation import batch
    from tcren.structure import Structure

    structs = [Structure("one", []), Structure("two", [])]
    seen = []
    monkeypatch.setattr(structure, "iter_structures", lambda *a, **k: iter(zip(("one", "two"), structs)))
    def annotate(items, **kwargs):
        assert kwargs["threads"] == 1
        seen.append(list(items))
        return [{"human": {}, "mouse": {}} for _ in items]
    monkeypatch.setattr(batch, "annotate_batch", annotate)
    def classify(s, **kwargs):
        assert kwargs["precomputed_records"] == {"human": {}, "mouse": {}}
        assert kwargs["organism"] == "mouse"
    monkeypatch.setattr(annotation, "classify_chains", classify)
    assert list(batch.iter_typed(Path("models.tar.gz"), organism="mouse")) == structs
    assert seen == [structs]


def test_scoring_full_ids_preserves_raw_stems_in_files_and_archives(tmp_path, monkeypatch):
    import tarfile

    import polars as pl
    from typer.testing import CliRunner
    from tcren import annotation, mhc, pipeline
    from tcren.annotation import batch
    from tcren.cli import app

    source_dir = tmp_path / "structures"
    source_dir.mkdir()
    fixture = Path(__file__).resolve().parents[1] / "assets/pdb/1ao7.pdb"
    names = ["first_full_stem", "second_full_stem"]
    for name in names:
        (source_dir / f"{name}.pdb").write_bytes(fixture.read_bytes())
    archive = tmp_path / "structures.tar.gz"
    with tarfile.open(archive, "w:gz") as stream:
        for path in sorted(source_dir.iterdir()):
            stream.add(path, arcname=path.name)
    calls = []
    monkeypatch.setattr(batch, "annotate_batch", lambda items, **kw: calls.append(len(items)) or [{} for _ in items])
    monkeypatch.setattr(annotation, "classify_chains", lambda *a, **kw: None)
    monkeypatch.setattr(mhc, "annotate_mhc_batch", lambda *a, **kw: None)
    monkeypatch.setattr(pipeline, "run", lambda structure, **kw: structure)
    monkeypatch.setattr(pipeline, "score_row", lambda s: {"pdb.id": s.pdb_id, "Phi_total": 1.0})
    for source in (source_dir, archive):
        out = tmp_path / "scores.csv"
        result = CliRunner().invoke(app, ["scoring", "-s", str(source), "-o", str(out),
                                         "--no-superimpose", "--full-ids"])
        assert result.exit_code == 0, result.output
        assert pl.read_csv(out)["pdb.id"].to_list() == names
    result = CliRunner().invoke(app, ["scoring", "-s", str(archive), "-o", str(out),
                                     "--no-superimpose"])
    assert result.exit_code == 0, result.output
    assert pl.read_csv(out)["pdb.id"].to_list() == ["first", "second"]
    assert calls == [2, 2, 2]


def test_all_structure_queries_share_one_multicore_call_per_reference(monkeypatch):
    """Never put per-structure MMseqs searches into a Python thread/process pool."""
    import os
    import threading
    from tcren import runtime
    from tcren.annotation import batch
    from tcren.structure import parse_structure

    fixture = Path(__file__).resolve().parents[1] / 'assets/pdb/1ao7.pdb'
    structures = [parse_structure(fixture, pdb_id=name) for name in ('one', 'two')]
    calls=[]
    def search(backend, pairs, organism, threads):
        calls.append((organism, list(pairs), threads, os.environ['OMP_THREAD_LIMIT'], threading.get_ident()))
        return [{} for _ in pairs]
    monkeypatch.setattr(runtime, 'available_cpus', lambda: 2)
    monkeypatch.setattr(batch, '_annotate', search)
    result = batch.annotate_batch(structures, arda=object(), threads=2)
    assert len(result) == 2
    assert [c[0] for c in calls] == ['human', 'mouse']
    expected=[(f'{i}|{c.chain_id}', c.sequence()) for i,s in enumerate(structures) for c in s.chains if c.sequence()]
    assert all(c[1] == expected and c[2:] == (2, '2', threading.get_ident()) for c in calls)
