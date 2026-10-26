"""Small checks for grouped profile evaluation, with no benchmark data."""

import numpy as np
from scipy.special import expit
from tcren.contacts.profile import connected_groups, fit_predict, separation


def test_profile_grouping_and_train_only_prediction():
    groups = connected_groups(["AAA", "AAA", "BBB", "CCC"], ["XXX", "YYY", "YYY", "ZZZ"])
    assert groups[0] == groups[1] == groups[2] != groups[3]
    x = np.array([[0.0, np.nan], [1.0, 1.0], [2.0, 0.0], [3.0, 1.0], [4.0, np.nan], [0.0, 0.0]])
    y = np.array([0, 0, 1, 1, 0, 1])
    train, test = np.arange(4), np.arange(4, 6)
    p, contributions, params = fit_predict(x, y, train, test, 0.1)
    flipped = y.copy()
    flipped[test] = 1 - flipped[test]
    p2, _, _ = fit_predict(x, flipped, train, test, 0.1)
    np.testing.assert_array_equal(p, p2)
    np.testing.assert_allclose(p, expit(params["intercept"] + contributions.sum(axis=1)))
    assert np.isfinite(p).all() and ((p > 0) & (p < 1)).all()
    assert params["imputation"] == [1.5, 2 / 3]
    assert separation([0, 0, 1, 1], [1, 1, 1, 1])["rank_biserial"] == 0
    assert separation([0, 0, 1, 1], [0, 0, 1, 1])["rank_biserial"] == 1


def test_feature_families_keep_id_order():
    import polars as pl
    from tcren.contacts.profile import ordered_matrix

    models = pl.DataFrame({"id": ["b", "a", "c"]})
    features = pl.DataFrame({"id": ["c", "b", "a"], "atoms": [30, 20, 10]})
    np.testing.assert_array_equal(ordered_matrix(models, features, ["atoms"]).ravel(), [20, 10, 30])


def test_contact_states_keep_endpoint_masses_and_gaps():
    import polars as pl
    from tcren.contacts.profile import contact_states, ordered_matrix

    prof = pl.DataFrame(
        dict(
            id=["zero", "partial", "full"],
            feature=["tcr:TRA:4"] * 3,
            contacted_atoms=[0, 2, 4],
            occupation=[0.0, 0.5, 1.0],
            fully_contacted=[0, 0, 1],
        )
    )
    wide = contact_states(prof).pivot(on="feature", index="id", values="value")
    names = ["tcr:TRA:4:" + s for s in ["engaged", "fully_contacted", "partial_occupation"]]
    x = ordered_matrix(pl.DataFrame(dict(id=["zero", "partial", "full", "gap"])), wide, names)
    np.testing.assert_array_equal(x[:3], [[0, 0, 0], [1, 0, 0.5], [1, 1, 0]])
    assert np.isnan(x[3]).all()


def test_digest_supports_python310_and_streams_all_blocks(tmp_path, monkeypatch):
    import hashlib
    from tcren.contacts.profile import digest

    monkeypatch.delattr(hashlib, "file_digest", raising=False)
    payload = b"contact-reference" * 70000
    path = tmp_path / "reference.bin"
    path.write_bytes(payload)
    assert digest(str(path)) == hashlib.sha256(payload).hexdigest()
