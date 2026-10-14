"""Batch dispatch regressions without sequence-search subprocesses."""
from pathlib import Path


def test_archive_is_annotated_in_one_batch(monkeypatch):
    from tcren import annotation, structure
    from tcren.annotation import batch
    from tcren.structure import Structure

    structs = [Structure("one", []), Structure("two", [])]
    seen = []
    monkeypatch.setattr(structure, "iter_structures", lambda *a, **k: iter(zip(("one", "two"), structs)))
    def annotate(items):
        seen.append(list(items))
        return [{"human": {}, "mouse": {}} for _ in items]
    monkeypatch.setattr(batch, "annotate_batch", annotate)
    def classify(s, **kwargs):
        assert kwargs["precomputed_records"] == {"human": {}, "mouse": {}}
        assert kwargs["organism"] == "mouse"
    monkeypatch.setattr(annotation, "classify_chains", classify)
    assert list(batch.iter_typed(Path("models.tar.gz"), organism="mouse")) == structs
    assert seen == [structs]
