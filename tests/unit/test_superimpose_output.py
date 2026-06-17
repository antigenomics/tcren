"""superimpose `-o` resolution: directory vs single file + extension/flag validation."""

from __future__ import annotations

import pytest

from tcren.orient.pipeline import _output_target


def test_directory_output_returns_none():
    assert _output_target("oriented", 5, mmcif=False, compress=False) is None
    assert _output_target("some/dir", 1, mmcif=True, compress=True) is None


def test_single_file_extension_must_match_flags():
    assert _output_target("x.pdb", 1, mmcif=False, compress=False).name == "x.pdb"
    assert _output_target("x.cif", 1, mmcif=True, compress=False).name == "x.cif"
    assert _output_target("x.cif.gz", 1, mmcif=True, compress=True).name == "x.cif.gz"
    assert _output_target("x.pdb.gz", 1, mmcif=False, compress=True).name == "x.pdb.gz"


def test_mismatched_extension_rejected():
    with pytest.raises(ValueError, match="--mmCIF"):
        _output_target("x.pdb", 1, mmcif=True, compress=False)
    with pytest.raises(ValueError, match="gz"):
        _output_target("x.cif", 1, mmcif=True, compress=True)
    with pytest.raises(ValueError, match="gz"):
        _output_target("x.pdb.gz", 1, mmcif=False, compress=False)


def test_single_file_with_many_inputs_rejected():
    with pytest.raises(ValueError, match="single file but 3"):
        _output_target("x.pdb", 3, mmcif=False, compress=False)
