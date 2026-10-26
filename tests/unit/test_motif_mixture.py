"""Exact Fisher arithmetic and the three-proportion score."""

import numpy as np
import polars as pl
import pytest
from scipy.stats import combine_pvalues

from tcren.contacts.motif_mixture import (
    fisher_summary,
    empirical_survival,
    simplex_weights,
    linear_score,
    fit_weights,
)


def test_fisher_statistic_and_mock_tail():
    data = pl.DataFrame({"id": ["a", "a", "b"], "tail_probability": [0.1, 0.2, 1.0]})
    result = fisher_summary(data).sort("id")
    expected = combine_pvalues([0.1, 0.2])
    assert result["fisher_statistic"][0] == pytest.approx(expected.statistic)
    assert result["fisher_independent_p"][0] == pytest.approx(expected.pvalue)
    assert result["fisher_statistic"][1] == 0
    np.testing.assert_allclose(empirical_survival([0, 3, 5], [0, 2, 4], [1, 1, 1]), [1, 0.5, 0.25])
    np.testing.assert_allclose(empirical_survival([3], [0, 2, 4], [1, 1, 2]), [0.6])


def test_three_weights_and_bounded_score():
    weights = simplex_weights(10)
    np.testing.assert_allclose(weights.sum(1), 1)
    assert len(weights) == 66 and (weights >= 0).all()
    x = np.array([[0.9, 0.1, 0.1], [0.8, 0.2, 0.2], [0.2, 0.8, 0.8], [0.1, 0.9, 0.9]])
    y = np.array([1, 1, 0, 0])
    selected, ap = fit_weights(x, y, np.array(["e"] * 4), threads=1, weights=weights)
    batch = linear_score(x, weights)
    for i in range(len(weights)):
        np.testing.assert_array_equal(batch[:, i], linear_score(x, weights[i]))
    np.testing.assert_array_equal(selected, [1, 0, 0])
    assert ap.max() == 1
    np.testing.assert_allclose(linear_score([[0, 0, 0], [1, 1, 1]], [0.5, 0.3, 0.2]), [0, 1])
    with pytest.raises(ValueError, match="sum to one"):
        linear_score(x, [1, 1, 1])


def test_fisher_cli_retains_all_controls_and_uses_model_leave_out(tmp_path):
    import json
    from tcren.contacts.residue_surprise import fit_reference
    from tcren.contacts.profile import digest
    from tcren.contacts.motif_mixture import fisher_command

    def rows(ids, counts):
        result = []
        for model, k in zip(ids, counts):
            for rank, (chain, region) in enumerate(
                [(c, r) for c in ["TRA", "TRB"] for r in ["CDR1", "CDR2", "CDR3"]]
                + [("PEPTIDE", "PEPTIDE")]
            ):
                result.append(
                    dict(
                        id=model,
                        chain_id=chain,
                        residue_index=rank + 1,
                        chain=chain,
                        region=region,
                        aa="A",
                        available=5,
                        complete=True,
                        occupied=int(k > 0),
                        pairs=k + int(k > 0),
                        cycles=k,
                    )
                )
        return pl.DataFrame(result)

    control = rows(["a", "b", "c"], [0, 1, 2])
    query = rows(["q"], [1])
    control.write_parquet(tmp_path / "controls.parquet")
    query.write_parquet(tmp_path / "query.parquet")
    resource = fit_reference(control)
    resource.update(reference_label=0, input_sha256=digest(tmp_path / "controls.parquet"))
    (tmp_path / "reference.json").write_text(json.dumps(resource))
    pl.DataFrame({"id": ["a", "b", "c"], "weight": [1.0, 1.0, 1.0]}).write_csv(
        tmp_path / "reference.models.tsv", separator="\t"
    )
    fisher_command(
        features=tmp_path / "query.parquet",
        reference=tmp_path / "reference.json",
        controls=tmp_path / "controls.parquet",
        output=tmp_path / "scores",
        statistic="cycles",
    )
    scored = pl.read_csv(tmp_path / "scores/scores.tsv", separator="\t")
    assert scored["mock_fisher_p"].item() == 0.75
    assert scored["enrichment"].item() == 0.25
    assert pl.read_parquet(tmp_path / "scores/control_fisher.parquet").height == 3
    pl.DataFrame(
        {"id": ["a", "b", "c"], "weight": [1.0] * 3, "PEPTIDE": ["A", "B", "B"]}
    ).write_csv(tmp_path / "reference.models.tsv", separator="\t")
    pl.DataFrame({"id": ["q"], "PEPTIDE": ["A"]}).write_csv(
        tmp_path / "query.identities.tsv", separator="\t"
    )
    fisher_command(
        features=tmp_path / "query.parquet",
        reference=tmp_path / "reference.json",
        controls=tmp_path / "controls.parquet",
        output=tmp_path / "conditioned",
        statistic="cycles",
        condition_epitope=True,
    )
    conditional = pl.read_csv(tmp_path / "conditioned/scores.tsv", separator="\t")
    assert conditional["mock_fisher_p"].item() == 0.5
    assert conditional["null_models"].item() == 1
    assert conditional["epitope_specific_null"].item() is True
    assert pl.read_parquet(tmp_path / "conditioned/control_fisher.parquet").height == 3
    positive = dict(resource, reference_label=1)
    (tmp_path / "positive.json").write_text(json.dumps(positive))
    fisher_command(
        features=tmp_path / "query.parquet",
        reference=tmp_path / "reference.json",
        controls=tmp_path / "controls.parquet",
        output=tmp_path / "ratio",
        statistic="cycles",
        positive_reference=tmp_path / "positive.json",
    )
    ratio = pl.read_csv(tmp_path / "ratio/scores.tsv", separator="\t")
    assert ratio["fisher_statistic"].item() == 0
    assert (
        pl.read_parquet(tmp_path / "ratio/residues.parquet")["evidence_bits"].to_list() == [0.0] * 7
    )


def test_diagnostic_isolates_weight_shift_on_same_rows(tmp_path):
    import json
    from tcren.contacts.motif_mixture import diagnose_command

    rows = pl.DataFrame(
        {
            "id": ["a", "b", "c", "d"],
            "epitope": ["e"] * 4,
            "y": [1, 0, 1, 0],
            "iptm": [0.9, 0.8, 0.9, 0.8],
            "betti": [0.1] * 4,
            "enrichment": [0.2] * 4,
            "prepared_scored": [0.9, 0.8, 0.9, 0.8],
            "mixture_oof": [0.9, 0.8, 0.5, 0.45],
            "mixture_fullfit": [0.9, 0.8, 0.9, 0.8],
        }
    )
    rows.write_parquet(tmp_path / "pred.parquet")
    pl.DataFrame({"id": rows["id"], "fold": [0, 0, 1, 1]}).write_parquet(tmp_path / "folds.parquet")
    (tmp_path / "weights.json").write_text(
        json.dumps({"weights": [1, 0, 0], "imputation": [0.85, 0.1, 0.2]})
    )
    (tmp_path / "fits.json").write_text(
        json.dumps([{"fold": 0, "weights": [1, 0, 0]}, {"fold": 1, "weights": [0.5, 0.5, 0]}])
    )
    diagnose_command(
        predictions=tmp_path / "pred.parquet",
        folds=tmp_path / "folds.parquet",
        fold_weights=tmp_path / "fits.json",
        resource=tmp_path / "weights.json",
        output=tmp_path / "out",
    )
    mean = pl.read_csv(tmp_path / "out/fold_mean.tsv", separator="\t")
    assert mean["AP"].to_list() == [1.0] * 3
    shifts = pl.read_csv(tmp_path / "out/scale_shifts.tsv", separator="\t")
    np.testing.assert_allclose(shifts["common_mean_shift"], [0, -0.375])


def test_global_fit_needs_no_folds_or_residue_inventory(tmp_path):
    import json
    from tcren.contacts.motif_mixture import fit_command, score_command

    ids = ["a", "b", "c", "d", "v", "ambiguous"]
    i = np.array([0.9, 0.8, 0.2, 0.1, 0.7, 0.6])
    b = 1 - i
    pl.DataFrame(
        dict(
            id=ids,
            dataset=["tcrvdb"] * 4 + ["vdjdb", "tcrvdb"],
            stratum=["all"] * 6,
            epitope=["e"] * 6,
            y=[1, 1, 0, 0, 1, None],
            iptm=i,
            prepared=0.55 * i + 0.45 * b,
        )
    ).write_parquet(tmp_path / "base.parquet")
    pl.DataFrame(
        {"pdb.id": ids, "complex_cyclo_frac": b, "receptor_contacts_available": [True] * 6}
    ).write_csv(tmp_path / "graph.tsv", separator="\t")
    pl.DataFrame({"id": ids, "enrichment": b}).write_csv(tmp_path / "motifs.tsv", separator="\t")
    fit_command(
        baseline=tmp_path / "base.parquet",
        topology=[tmp_path / "graph.tsv"],
        motifs=tmp_path / "motifs.tsv",
        residues=None,
        folds=None,
        threads=1,
        output=tmp_path / "fit",
    )
    r = json.loads((tmp_path / "fit/weights.json").read_text())
    assert r["weights"] == [1.0, 0.0, 0.0] and r["training_n"] == 4
    assert r["assessment"].startswith("No folds;")
    assert json.loads((tmp_path / "fit/fold_weights.json").read_text()) == []
    score_command(
        features=tmp_path / "fit/predictions.parquet",
        resource=tmp_path / "fit/weights.json",
        output=tmp_path / "scores.tsv",
    )
    out = pl.read_csv(tmp_path / "scores.tsv", separator="\t")
    assert out.height == 6 and out["mixture_oof"].null_count() == 6
    np.testing.assert_allclose(out["score"], out["iptm"])


def test_cohort_fit_uses_only_its_labels_and_dispatches_frozen_weights(tmp_path):
    import json
    from tcren.contacts.motif_mixture import fit_cohorts_command, score_command

    data = []
    for k, (dataset, stratum) in enumerate(
        [("tcrvdb", "all"), ("vdjdb", "template_present"), ("vdjdb", "template_absent")]
    ):
        for j, y in enumerate([1, 1, 0, 0]):
            x = [0.9, 0.8, 0.2, 0.1][j]
            v = [1 - x] * 3
            v[k] = x
            data.append(
                dict(
                    id=f"{k}-{j}",
                    dataset=dataset,
                    stratum=stratum,
                    epitope="e",
                    y=y,
                    iptm=v[0],
                    betti=v[1],
                    enrichment=v[2],
                )
            )
    pl.DataFrame(data).write_parquet(tmp_path / "components.parquet")
    fit_cohorts_command(
        features=tmp_path / "components.parquet", output=tmp_path / "fit", threads=1
    )
    model = json.loads((tmp_path / "fit/weights.json").read_text())
    assert [g["weights"] for g in model["groups"]] == [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
    assert [g["training_n"] for g in model["groups"]] == [4] * 3
    score_command(
        features=tmp_path / "fit/predictions.parquet",
        resource=tmp_path / "fit/weights.json",
        output=tmp_path / "scores.tsv",
    )
    s = pl.read_csv(tmp_path / "scores.tsv", separator="\t")
    np.testing.assert_allclose(s["score"], s["mixture_fullfit"])
    assert (s.filter(pl.col("y") == 1)["score"] >= 0.8).all()
    bad = pl.DataFrame(data).with_columns(pl.lit("unknown").alias("dataset"))
    bad.write_parquet(tmp_path / "bad.parquet")
    with pytest.raises(ValueError, match="exactly one"):
        score_command(
            features=tmp_path / "bad.parquet",
            resource=tmp_path / "fit/weights.json",
            output=tmp_path / "bad.tsv",
        )
