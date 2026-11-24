import pickle
import numpy as np
import pytest
from tcren.validation.tcrmodel2 import retain_generated_inputs


def test_retained_inputs_preserve_producer_bytes_and_reject_mismatch(tmp_path):
    features = tmp_path / "features.pkl"
    msa = tmp_path / "msa.pkl"
    out = tmp_path / "retained"
    features.write_bytes(
        pickle.dumps(
            {
                "aatype": np.zeros(3, dtype=np.int32),
                "template_aatype": np.zeros((4, 3), dtype=np.int32),
            },
            protocol=4,
        )
    )
    msa.write_bytes(pickle.dumps({"msa": np.zeros((4096, 3), dtype=np.int32)}, protocol=4))
    record = retain_generated_inputs(features, msa, out)
    assert record["tokens"] == 3
    assert (out / "producer-features.pkl").read_bytes() == features.read_bytes()
    assert (out / "producer-msa.pkl").read_bytes() == msa.read_bytes()
    msa.write_bytes(pickle.dumps({"msa": np.zeros((4096, 4), dtype=np.int32)}, protocol=4))
    with pytest.raises(ValueError, match="dimensions differ"):
        retain_generated_inputs(features, msa, tmp_path / "invalid")
    assert not (tmp_path / "invalid").exists()
    msa.write_bytes((out / "producer-msa.pkl").read_bytes())
    (out / "producer-features.pkl").write_bytes(b"previous retained bytes")
    with pytest.raises(ValueError, match="Refuse to replace"):
        retain_generated_inputs(features, msa, out)
    assert (out / "producer-features.pkl").read_bytes() == b"previous retained bytes"


def test_feature_array_identity_distinguishes_values_shape_and_dtype():
    from tcren.validation.tcrmodel2 import compare_generated_feature_arrays
    import numpy as np

    left = {"msa": np.array([[1, 2]], dtype=np.int32), "template": np.array([b"public"])}
    assert compare_generated_feature_arrays(left, left)["identical"]
    changed = {**left, "msa": np.array([[1, 3]], dtype=np.int32)}
    result = compare_generated_feature_arrays(left, changed)
    assert not result["identical"]
    assert result["arrays"][0]["maximum_absolute_difference"] == 1
    typed = compare_generated_feature_arrays(left, {**left, "msa": left["msa"].astype(np.int64)})
    assert typed["arrays"][0]["values_equal"] and not typed["identical"]
    shaped = compare_generated_feature_arrays(left, {"msa": np.array([1, 2])})
    assert not shaped["identical"] and shaped["arrays"][1]["missing_from"] == "right"


def test_msa_inventory_separates_permutation_from_pairing_changes():
    from tcren.validation.tcrmodel2 import compare_msa_row_inventory
    import numpy as np

    left = np.array([[1, 7], [2, 8]])
    moved = compare_msa_row_inventory(left, left[::-1], chain_ids=[1.0, 2.0])
    assert moved["full"]["multiset_equal"] and not moved["full"]["order_equal"]
    paired = compare_msa_row_inventory(left, np.array([[1, 8], [2, 7]]), chain_ids=[1, 2])
    assert not paired["full"]["multiset_equal"]
    assert all(row["multiset_equal"] for row in paired["chains"].values())
