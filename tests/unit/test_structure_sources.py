"""gzip / tar.gz structure inputs + identifier resolution (no mmseqs)."""

from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import pytest

from tcren.structure import (
    iter_structures,
    parse_structure,
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
