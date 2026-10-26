"""Historical source/cluster/calibration adapters preserve established calculations."""
import numpy as np
import polars as pl
import pytest

from tcren.calibration import epitope_weights
from tcren.potential.redundancy import _cluster, cluster_table
from tcren.structure.repair import restore_modified_peptide


def test_clusters_reuse_historical_membership_and_linkage():
    contacts = pl.DataFrame({"pdb.id": ["a", "b", "c"], "chain.type.from": ["TRA", "TRB", "TRD"]})
    markup = pl.DataFrame({"pdb.id": ["a", "b", "c"], "cdr3a": ["AA", "AA", "GG"],
                           "cdr3b": ["TT", "TT", "GG"], "peptide": ["AAAA", "AAAA", "GGGG"]})
    expected_ids, expected = _cluster(markup.head(2), 6., ("cdr3a", "cdr3b", "peptide"), "complete")
    result = cluster_table(contacts, markup)
    assert result.filter(pl.col("alpha_beta"))["pdb.id"].to_list() == expected_ids
    np.testing.assert_array_equal(result.filter(pl.col("alpha_beta"))["cluster"], expected)
    assert result.filter(pl.col("pdb.id") == "c")["cluster"].item() is None
    assert epitope_weights(markup.rename({"peptide": "epitope"})) == {"a": .5, "b": .5, "c": 1.}


def line(index, xyz, chain="C", residue="ALA", record="ATOM"):
    return f"{record:<6}{index:5d}  CA  {residue} {chain}{index:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 90.00           C"


def test_deposited_restore_preserves_nonpeptide_and_rigid_identity():
    coordinates = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    native = [line(i+1, x) for i, x in enumerate(coordinates)]
    other = line(8, [2., 2., 2.], chain="A")
    native.append(other)
    unrelated = line(9, [3., 3., 3.], residue="MSE", record="HETATM")
    native.append(unrelated)
    source = [line(i+1, x+[5., -2., 1.], chain="X") for i, x in enumerate(coordinates)]
    source.append(line(4, [5., -2., 2.], chain="X", residue="CIR", record="HETATM"))
    restored, record = restore_modified_peptide(native, source)
    assert unrelated in restored
    assert other in restored and record["matched_standard_ca"] == 3
    assert record["restored_modified_atoms"] == 1
    modified = next(s for s in restored.splitlines() if s.startswith("HETATM"))
    assert modified[21] == "C"
    np.testing.assert_allclose([float(modified[k:k+8]) for k in (30, 38, 46)], [0., 0., 1.], atol=1e-3)
    with pytest.raises(ValueError, match="modified"):
        restore_modified_peptide(native, source[:-1])


@pytest.mark.parametrize("coordinates", [[[0., 0., 0.]], [[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]])
def test_restore_rejects_underdetermined_orientation(coordinates):
    native = [line(i+1, x) for i, x in enumerate(coordinates)]
    with pytest.raises(ValueError, match="noncollinear"):
        restore_modified_peptide(native, native)
