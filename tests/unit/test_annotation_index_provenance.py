"""Searched-index identity survives batching and rejects changes during a search."""
import hashlib
from types import SimpleNamespace

import pytest

from tcren.annotation import batch
from tcren.structure.model import Atom, Chain, Residue, Structure


def structure():
    return Structure("query", [Chain("A", [Residue(0, 1, "", "A", "ALA", (Atom("CA", "C", (0, 0, 0)),))])])


@pytest.mark.parametrize("mutate", [False, True])
def test_annotation_records_searched_index_and_rejects_midsearch_changes(tmp_path, monkeypatch, mutate):
    index = tmp_path / "db"
    index.write_bytes(b"first generation")

    def identity(organism, seqtype="aa"):
        return {"database_path": str(index), "files_sha256": {
            "db": hashlib.sha256(index.read_bytes()).hexdigest()}}

    def annotate(*args):
        if mutate:
            index.write_bytes(b"another generation")
        return [{"locus": "TRA"}]

    monkeypatch.setattr(batch, "_annotate", annotate)
    item = structure()
    backend = SimpleNamespace(reference_index_identity=identity)
    if mutate:
        with pytest.raises(RuntimeError, match="index.*changed"):
            batch.annotate_batch([item], backend, organisms=("human",))
    else:
        batch.annotate_batch([item], backend, organisms=("human",))
        assert item.annotation_provenance["arda"]["human"]["files_sha256"]["db"] == (
            "f144c659e319d8bf3b6c590cf4f0286b3d4d6b9109113c60c807e1468384165a")


def test_index_provenance_keeps_split_payloads_and_rejects_mixed_batch(tmp_path):
    from tcren import provenance

    for name, data in (("db.0", b"one"), ("db.1", b"two"), ("db.index", b"index"),
                       ("db_h", b"headers"), ("unrelated", b"ignore")):
        (tmp_path / name).write_bytes(data)
    identity = provenance.searched_index_identity(tmp_path / "db")
    assert set(identity["files_sha256"]) == {"db.0", "db.1", "db.index", "db_h"}
    first, second = structure(), structure()
    first.annotation_provenance["mhc"] = identity
    second.annotation_provenance["mhc"] = {**identity, "database_path": "another-index"}
    with pytest.raises(RuntimeError, match="mixed annotation indexes"):
        provenance.annotation_indexes([first, second])


def test_feature_based_recognition_preserves_searched_index_identity(tmp_path):
    import polars as pl
    from tcren.cli import _score_feature_table
    from tcren.provenance import read, stamp

    source, destination = tmp_path / "features.tsv", tmp_path / "scores.tsv"
    pl.DataFrame({"complex.id": ["synthetic"]}).write_csv(source, separator="\t")
    identity = {"arda": {"human": {"files_sha256": {"db": "a-recorded-hash"}}}}
    stamp(source, command="synthetic input", extra={"searched_indexes": identity})
    _score_feature_table(source, destination)
    assert read(destination)["searched_indexes"] == identity
