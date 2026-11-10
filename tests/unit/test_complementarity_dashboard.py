import gzip
import json

import pytest

from tcren.complementarity.dashboard import write_dashboard
from tcren.complementarity.frame import CanonicalFrame
from tcren.complementarity.scene import build_scene
from tcren.structure.model import Atom, Chain, Residue, Structure
import numpy as np
import polars as pl


def empty_scene(identity):
    return build_scene(Structure(identity, []),
                       CanonicalFrame(np.eye(3), np.zeros(3), "ref", "apex"))


def test_streamed_scene_identity_is_never_a_path(tmp_path):
    scene = empty_scene("../unsafe</script>")
    out = tmp_path / "browser"
    path = write_dashboard(iter([scene]), out)
    assert path.is_file()
    assert json.loads(gzip.decompress((out / "scene-0.json.gz").read_bytes())) == scene
    assert json.loads((out / "index.json").read_text())[0]["id"] == scene["id"]
    assert "../unsafe" not in path.read_text()
    with pytest.raises(FileExistsError):
        write_dashboard([], out)


def test_nan_export_fails_explicitly(tmp_path):
    scene = empty_scene("bad")
    scene["cutoff"] = float("nan")
    with pytest.raises(ValueError):
        write_dashboard([scene], tmp_path / "browser")


def test_parquet_partitions_preserve_typed_empty_tables(tmp_path):
    out = tmp_path / "browser"
    write_dashboard((empty_scene(f"structure-{i}") for i in range(9)), out)
    manifest = json.loads((out / "tables/manifest.json").read_text())
    assert len(manifest["files"]) == 10
    for item in manifest["files"]:
        table = pl.read_parquet(out / "tables" / item["file"])
        assert table.height == 0
        assert table.schema["structure_id"] == pl.String
    atoms = pl.scan_parquet(out / "tables/atoms/*.parquet").collect()
    assert atoms.schema["x"] == pl.Float64
    assert atoms.schema["contact_count"] == pl.Int64


def test_json_and_parquet_share_counts_and_ca_coordinate_owner(tmp_path):
    chains = []
    for index, x in enumerate((0., 3., 30.)):
        residue = Residue(0, 10 + index, "", "G", "GLY",
                          (Atom("CA", "C", np.array([x, 0., 0.])),))
        chains.append(Chain(str(index), [residue], "PEPTIDE"))
    scene = build_scene(Structure("molecule", chains),
                        CanonicalFrame(np.eye(3), np.zeros(3), "ref", "apex"))
    out = tmp_path / "browser"
    write_dashboard([scene], out)
    atoms = pl.read_parquet(out / "tables/atoms/part-00000.parquet")
    residues = pl.read_parquet(out / "tables/residues/part-00000.parquet")
    edges = pl.read_parquet(out / "tables/residue_contacts/part-00000.parquet")
    assert atoms["contact_count"].to_list() == [1, 1, 0]
    assert residues["contact_count"].to_list() == [1, 1, 0]
    assert edges["contact_count"].to_list() == [1]
    joined = residues.join(atoms.select("structure_id", pl.col("id").alias("ca_atom"), "x"),
                           on=["structure_id", "ca_atom"], validate="1:1")
    assert joined["x"].to_list() == [0., 3., 30.]
    assert "x" not in residues.columns
