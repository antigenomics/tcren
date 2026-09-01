"""gzip / tar.gz structure inputs + identifier resolution (no mmseqs)."""

from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import pytest

from tcren.structure import (
    iter_structures,
    parse_structure,
    resolve_structure_ids,
    structure_id_from_path,
)

_ASSET = Path(__file__).resolve().parents[1] / "assets" / "cgene" / "1ao7_full.pdb"


@pytest.mark.parametrize(
    "name,expected",
    [
        ("1ao7.pdb", "1ao7"),
        ("1ao7.pdb.gz", "1ao7"),
        ("4x5w_renumbered.cif", "4x5w"),
        ("4x5w_renumbered.cif.gz", "4x5w"),
        ("6uk4_TCRpMHCmodels_polyV.pdb", "6uk4"),
    ],
)
def test_structure_id_from_path(name, expected):
    assert structure_id_from_path(name) == expected


def test_parse_gzipped_pdb_matches_plain(tmp_path):
    plain = parse_structure(_ASSET, pdb_id="1ao7")
    gz = tmp_path / "1ao7.pdb.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(_ASSET.read_text())
    out = parse_structure(gz)
    assert out.pdb_id == "1ao7"  # id resolved from filename, .gz + .pdb stripped
    assert [c.chain_id for c in out.chains] == [c.chain_id for c in plain.chains]
    assert sum(len(c.residues) for c in out.chains) == sum(len(c.residues) for c in plain.chains)


def test_iter_structures_dir_and_targz(tmp_path):
    # A directory with a plain and a gzipped copy.
    d = tmp_path / "structs"
    d.mkdir()
    (d / "1ao7.pdb").write_text(_ASSET.read_text())
    with gzip.open(d / "2xyz.pdb.gz", "wt") as fh:
        fh.write(_ASSET.read_text())
    got = dict(iter_structures(d, importer=parse_structure))
    assert set(got) == {"1ao7", "2xyz"}

    # A .tar.gz archive of both is streamed and parsed.
    tgz = tmp_path / "batch.tar.gz"
    with tarfile.open(tgz, "w:gz") as tar:
        tar.add(d / "1ao7.pdb", arcname="1ao7.pdb")
        tar.add(d / "2xyz.pdb.gz", arcname="nested/2xyz.pdb.gz")
    ids = {pid for pid, _ in iter_structures(tgz, importer=parse_structure)}
    assert ids == {"1ao7", "2xyz"}


def test_resolve_sources_splits_commas_and_repeats():
    from tcren.structure.io import resolve_sources

    assert resolve_sources("a.pdb") == ["a.pdb"]
    assert resolve_sources("a.pdb,b.pdb.gz") == ["a.pdb", "b.pdb.gz"]
    assert resolve_sources(["a.pdb,b.pdb", "c.cif"]) == ["a.pdb", "b.pdb", "c.cif"]


def test_structure_paths_glob_and_manifest(tmp_path):
    """A glob and a .txt manifest must resolve to the same files as a plain listing."""
    from tcren.structure.io import structure_paths

    for name in ("a.pdb", "b.pdb.gz", "notes.md"):
        (tmp_path / name).write_text("x")
    assert [p.name for p in structure_paths(str(tmp_path / "*"))] == ["a.pdb", "b.pdb.gz"]

    # relative manifest entries resolve against the manifest's own directory, not the CWD
    (tmp_path / "models.txt").write_text("# a comment\na.pdb\n\nb.pdb.gz\n")
    assert structure_paths(tmp_path / "models.txt") == [tmp_path / "a.pdb", tmp_path / "b.pdb.gz"]


def test_structure_ids_fall_back_to_stems_when_the_prefix_collides():
    """The PDB-id prefix is lossy: keep it only while it stays unique over the set.

    ``VDJdb_Model_603_min.pdb`` and ``VDJdb_Model_604_min.pdb`` both start ``VDJdb``, so the
    per-file rule collapsed 73 distinct models onto 7 identifiers and silently destroyed the
    rows downstream. Deciding per SET keeps the convenience where it is unambiguous.
    """
    rcsb = ["4x5w_renumbered.cif", "1ao7.pdb.gz", "6uk4_TCRpMHCmodels.pdb"]
    assert list(resolve_structure_ids(rcsb).values()) == ["4x5w", "1ao7", "6uk4"]

    colliding = ["VDJdb_Model_603_min.pdb", "VDJdb_Model_604_min.pdb"]
    assert list(resolve_structure_ids(colliding).values()) == [
        "VDJdb_Model_603_min", "VDJdb_Model_604_min"]

    # <hash>_<epitope>_<mhc> decoys keep the bare hash while the hashes stay unique
    decoys = ["aaa_GILGFVFTL_HLA-A02.pdb", "bbb_GILGFVFTL_HLA-A02.pdb"]
    assert list(resolve_structure_ids(decoys).values()) == ["aaa", "bbb"]


def test_iter_structures_never_yields_a_duplicate_id(tmp_path):
    """End to end: two files whose prefixes collide must come back as two distinct ids."""
    body = _ASSET.read_bytes()
    for n in ("Model_1_min.pdb", "Model_2_min.pdb"):
        (tmp_path / n).write_bytes(body)
    ids = sorted(i for i, _ in iter_structures(tmp_path, importer=parse_structure))
    assert ids == ["Model_1_min", "Model_2_min"]


def test_appledouble_sidecars_are_skipped(tmp_path):
    """Tarring a structure set on macOS writes ``._x.pdb`` beside every ``x.pdb``.

    The sidecar carries a binary resource fork under the extension of the file it shadows,
    so a parser reaches a decode error rather than a structure; the corpus recompute lost
    two sets of 2,000+ models to exactly this.
    """
    from tcren.structure.io import is_structure_file

    assert is_structure_file("4x5w.pdb")
    assert not is_structure_file("._4x5w.pdb")
    assert not is_structure_file("positives/._4x5w.pdb.gz")

    d = tmp_path / "set"
    d.mkdir()
    (d / "1ao7.pdb").write_bytes(_ASSET.read_bytes())
    (d / "._1ao7.pdb").write_bytes(b"\x00\x05\x16\x07\xa3Mac OS X resource fork")
    assert [pid for pid, _ in iter_structures(d, importer=parse_structure)] == ["1ao7"]

    tgz = tmp_path / "set.tar.gz"
    with tarfile.open(tgz, "w:gz") as tar:
        for p in sorted(d.iterdir()):
            tar.add(p, arcname=p.name)
    assert [pid for pid, _ in iter_structures(tgz, importer=parse_structure)] == ["1ao7"]
