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
    assert len(manifest["files"]) == 12
    for item in manifest["files"]:
        table = pl.read_parquet(out / "tables" / item["file"])
        expected = (8 if '00000' in item['file'] else 1) if item['table'] == 'structures' else 0
        assert table.height == expected
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
    from polars.testing import assert_frame_equal

    restored = pl.read_csv(out / "tables/atoms/part-00000.tsv.gz", separator="\t", schema_overrides=atoms.schema)
    assert_frame_equal(atoms, restored)
    summary = pl.read_parquet(out / "summary.parquet")
    summary_tsv = pl.read_csv(out / "summary.tsv.gz", separator="\t", schema_overrides=summary.schema)
    assert_frame_equal(summary, summary_tsv)
    report = json.loads((out / "summary.json").read_text())
    assert report["totals"]["atoms"] == 3
    assert report["totals"]["atom_contacts"] == 1
    from tcren.complementarity.summary import summarize_exports

    checked = summarize_exports(out)
    assert checked["row_and_incidence_checks"] == "passed"
    assert checked["totals"]["atoms"] == 3
    modified = atoms.with_columns(pl.lit(0).alias("contact_count"))
    modified.write_parquet(out / "tables/atoms/part-00000.parquet")
    with pytest.raises(ValueError, match="incidence differs"):
        summarize_exports(out)


def test_duplicate_structure_ids_rejected(tmp_path):
    with pytest.raises(ValueError, match="duplicate structure ID"):
        write_dashboard([empty_scene("same"), empty_scene("same")], tmp_path / "browser")


def test_summary_checks_peptide_contacts_against_residue_roles(tmp_path):
    from tcren.complementarity.summary import summarize_exports

    out = tmp_path / "browser"
    write_dashboard([empty_scene("empty")], out)
    assert summarize_exports(out)["zero_tcr_peptide_interfaces"] == 1
    summary = pl.read_parquet(out / "summary.parquet").with_columns(
        pl.lit(1, dtype=pl.Int64).alias("tcr_peptide_atom_contacts"))
    summary.write_parquet(out / "summary.parquet")
    with pytest.raises(ValueError, match="TCR–peptide contacts differ"):
        summarize_exports(out)
