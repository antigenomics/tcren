"""Check exact tied AP and physical atom identity in additive contact scores."""

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score
from tcren.contacts.additive import average_precision_matrix, contact_counts, training_percentile


def test_ties_counts_and_train_only_transform():
    y = np.array([1, 0, 1, 0, 1, 0])
    scores = np.array([[1, 0], [1, 1], [2, 2], [3, 2], [3, 3], [0, 3]], float)
    np.testing.assert_allclose(
        average_precision_matrix(y, scores), [average_precision_score(y, s) for s in scores.T]
    )
    c = pl.DataFrame(
        {
            "pdb.id": ["x"] * 3,
            "interface": ["tcr_peptide"] * 3,
            "role.from": ["tcr"] * 3,
            "role.to": ["peptide"] * 3,
            "dist": [3.0] * 3,
            "chain.id.from": ["A", "A", "B"],
            "chain.type.from": ["TRA", "TRA", "TRB"],
            "residue.index.from": [1, 1, 1],
            "atom.from": ["CA", "CB", "CA"],
            "chain.id.to": ["P"] * 3,
            "residue.index.to": [1] * 3,
            "atom.to": ["O", "O", "O"],
        }
    )
    r = contact_counts(c, pl.DataFrame({"pdb.id": ["x", "zero"]}))
    assert r.select("max_pairs", "max_atoms", "peptide_atoms").rows() == [(2, 2, 1), (0, 0, 0)]
    assert r["pairs_per_peptide_atom"].to_list() == [2.0, 0.0]
    np.testing.assert_allclose(training_percentile([0, 1, 1, 2], [-1, 1, 3]), [0, 0.5, 1])


def test_weight_selection_does_not_read_held_out_labels():
    from tcren.contacts.additive import FEATURES, candidate_spec, scores_for, select

    models = pl.DataFrame(
        {
            "y": [0, 1, 0, 1, 0, 1],
            "cohort": ["E"] * 6,
            "iptm": [0.2, 0.8, 0.3, 0.7, 0.4, 0.6],
            "betti": [0.1] * 6,
            **{name: [0.0, 2.0, 1.0, 3.0, 2.0, 4.0] for name in FEATURES},
        }
    )
    specs = candidate_spec()
    scores = scores_for(models, specs)
    train = np.arange(4)
    before, ap = select(models, scores, specs, train)
    changed = models.with_columns(pl.Series("y", [0, 1, 0, 1, 1, 0]))
    after, changed_ap = select(changed, scores, specs, train)
    assert before == after
    np.testing.assert_array_equal(ap, changed_ap)
