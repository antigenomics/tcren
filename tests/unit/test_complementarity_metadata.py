import tarfile
from pathlib import Path

import polars as pl
import pytest

from tcren.complementarity.command import _structures
from tcren.complementarity.inputs import metadata_inputs


PDB = Path(__file__).parents[1] / "assets/pdb/1ao7.pdb"


def test_text_metadata_preserves_numeric_looking_ids(tmp_path):
    path = tmp_path / "input.tsv"
    path.write_text(f"id\tstructure_path\n001\t{PDB}\n1\t{PDB}\n")
    assert [r.identity for r in metadata_inputs(path)] == ["001", "1"]


def test_metadata_formats_relative_paths_and_explicit_identity(tmp_path):
    (tmp_path / "model.pdb").write_bytes(PDB.read_bytes())
    table = pl.DataFrame({"id": ["custom:alpha"], "structure_path": ["model.pdb"], "outcome": [None]})
    targets = []
    for name in ["input.tsv.gz", "input.parquet"]:
        path = tmp_path / name
        if name.endswith("parquet"):
            table.write_parquet(path)
        else:
            import gzip
            path.write_bytes(gzip.compress(table.write_csv(separator="\t").encode()))
        rows = metadata_inputs(path)
        assert rows[0].identity == "custom:alpha"
        assert rows[0].path == (tmp_path / "model.pdb").resolve()
        structure = next(_structures(None, tmp_path, rows))
        assert structure.pdb_id == "custom:alpha"
        assert structure.annotation_provenance["visualization_input"]["metadata_row"] == 1
        targets.append(rows)
    assert targets[0] == targets[1]


def test_metadata_archive_aliases_and_missing_members(tmp_path):
    archive = tmp_path / "models.tar"
    with tarfile.open(archive, "w") as stream:
        stream.add(PDB, arcname="nested/native.pdb")
    path = tmp_path / "input.tsv"
    pl.DataFrame({"id": ["one", "two"], "structure_path": ["models.tar"]*2,
                  "member": ["nested/native.pdb"]*2}).write_csv(path, separator="\t")
    rows = metadata_inputs(path, member_column="member")
    assert [s.pdb_id for s in _structures(None, tmp_path, rows)] == ["one", "two"]
    pl.DataFrame({"id": ["missing"], "structure_path": ["models.tar"],
                  "member": ["not-present.pdb"]}).write_csv(path, separator="\t")
    with pytest.raises(ValueError, match="members missing"):
        list(_structures(None, tmp_path, metadata_inputs(path, member_column="member")))


@pytest.mark.parametrize("identities", [["same", "same"], ["", "two"], [None, "two"], [" white", "two"]])
def test_metadata_identity_errors_are_explicit(tmp_path, identities):
    path = tmp_path / "input.tsv"
    pl.DataFrame({"id": identities, "structure_path": [str(PDB)]*2}).write_csv(path, separator="\t")
    with pytest.raises(ValueError, match="identity"):
        metadata_inputs(path)
