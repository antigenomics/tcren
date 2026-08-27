"""metadata.tsv: the table that ships beside a structure set."""
import polars as pl
import pytest

from tcren.metadata import find_metadata, join_metadata, read_metadata


def _set(tmp_path, rows, name="metadata.tsv"):
    d = tmp_path / "myset"
    d.mkdir(exist_ok=True)
    pl.DataFrame(rows).write_csv(d / name, separator="\t")
    return d


def test_find_and_read(tmp_path):
    d = _set(tmp_path, {"id": ["a", "b"], "iptm": [0.8, 0.4]})
    assert find_metadata(d) == d / "metadata.tsv"
    assert read_metadata(d).height == 2


def test_found_from_a_file_inside_the_set(tmp_path):
    d = _set(tmp_path, {"id": ["a"], "iptm": [0.8]})
    (d / "a.pdb").write_text("")
    assert find_metadata(d / "a.pdb") == d / "metadata.tsv"


def test_absent_metadata_is_not_an_error(tmp_path):
    assert read_metadata(tmp_path) is None
    t = pl.DataFrame({"complex.id": ["a"], "burial": [1.0]})
    assert join_metadata(t, tmp_path).equals(t)      # unchanged, so callers need no branch


def test_missing_id_column_raises(tmp_path):
    d = _set(tmp_path, {"stem": ["a"], "iptm": [0.8]})
    with pytest.raises(ValueError, match="no 'id' column"):
        read_metadata(d)


def test_duplicate_ids_raise(tmp_path):
    """The defect this guard exists for: vdjdb_binder_benchmark shipped 1,089 rows keyed on a bare
    TCR hash with only 1,068 distinct values, so 566 negatives silently failed to join."""
    d = _set(tmp_path, {"id": ["a", "a"], "y": [1, 0]})
    with pytest.raises(ValueError, match="duplicate ids"):
        read_metadata(d)


def test_join_is_left_and_keyed_on_complex_id(tmp_path):
    d = _set(tmp_path, {"id": ["a", "b"], "y": [1, 0], "iptm": [0.8, 0.4]})
    t = pl.DataFrame({"complex.id": ["a", "c"], "burial": [1.0, 2.0]})
    j = join_metadata(t, d)
    assert j.height == 2                       # left join keeps the unmatched row
    assert j["y"].to_list() == [1, None]
    assert j["iptm"].to_list() == [0.8, None]


def test_column_subset(tmp_path):
    d = _set(tmp_path, {"id": ["a"], "y": [1], "iptm": [0.8], "plddt": [90.0]})
    j = join_metadata(pl.DataFrame({"complex.id": ["a"]}), d, columns=("iptm",))
    assert "iptm" in j.columns and "plddt" not in j.columns


def test_clashing_column_is_prefixed_not_dropped(tmp_path):
    d = _set(tmp_path, {"id": ["a"], "iptm": [0.8]})
    t = pl.DataFrame({"complex.id": ["a"], "iptm": [0.1]})
    j = join_metadata(t, d)
    assert j["iptm"][0] == 0.1 and j["meta.iptm"][0] == 0.8
