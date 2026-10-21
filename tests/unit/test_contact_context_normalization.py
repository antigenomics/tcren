"""Exact reference factor identities, neutral unseen support and no hidden omission."""
import copy

import numpy as np
import polars as pl
import pytest

from tcren.contacts.context_normalization import (
    VERSION, apply_normalized, fit_normalizer, fit_partition, normalized_edges, validate_normalizer,
)
from tcren.contacts.plausibility import VERSION as PRIOR_VERSION
from tcren.contacts.weighted import KEYS, STRATA


def edges(tails=(.2, .8, .4), counts=(2, 2, 3)):
    n = len(tails)
    return pl.DataFrame({"pdb.id": ["x"]*n, "chain.id.from": ["T"]*n,
        "residue.index.from": list(range(n)), "chain.id.to": ["P"]*n,
        "residue.index.to": [1]*n, "minimum_tail": list(tails),
        "n_distance": list(counts), **{key: [0]*n for key in STRATA[1:]}})


def normalizer(reference):
    return dict(version=VERSION, strata=STRATA, mhc_class="I", marginal_resource_version=PRIOR_VERSION,
        marginal_resource_sha256="abc", records=fit_partition(reference).to_dicts())


def test_reference_partition_identity_constant_cancellation_and_unseen_missing():
    reference = edges()
    model = normalizer(reference)
    validate_normalizer(model, "abc")
    for record in model["records"]:
        assert record["mean_reference_normalized_factor"] == pytest.approx(1., abs=2e-15)
    result = normalized_edges(reference, model)
    for _, part in result.group_by(STRATA):
        assert np.exp(-part["normalized_increment"].to_numpy()).mean() == pytest.approx(1., abs=2e-15)
    constant = edges((.3, .3, .3))
    np.testing.assert_allclose(normalized_edges(constant, normalizer(constant))["normalized_increment"], 0., atol=2e-15)
    query = edges((.5, None), (99, 2))
    result = normalized_edges(query, model)
    assert result["normalized_increment"].to_list() == [0., None]
    assert result["normalizer_status"].to_list() == ["unseen_multiplicity", "unavailable_geometry"]
    with pytest.raises(ValueError, match="different marginal"):
        validate_normalizer(model, "changed")
    bad = copy.deepcopy(model)
    bad["records"][0]["logZ"] = .01
    with pytest.raises(ValueError, match="finite logZ"):
        validate_normalizer(bad, "abc")
    with pytest.raises(ValueError, match="unique directed"):
        fit_partition(pl.concat([reference, reference]))


def test_directed_join_exact_sum_framework_and_missing_geometry():
    reference = edges()
    increments = normalized_edges(reference, normalizer(reference))
    table = reference.select(KEYS).with_columns(pl.Series("original_contribution", [-2., 3., 4.]),
        pl.Series("geometry_status", ["supported_cdr_geometry", "supported_cdr_geometry", "unassessed_framework"]))
    result, summary = apply_normalized(table, increments, {"original": 5.})
    assert summary["normalized"] == pytest.approx(5.+result["normalized_increment"].sum())
    assert result["normalized_increment"][2] == 0.
    assert result["normalized_contribution"][2] == 4.
    # Wrong direction is unavailable rather than borrowing the reversed edge.
    reversed_edges = increments.with_columns(pl.col("chain.id.to").alias("chain.id.from"), pl.col("chain.id.from").alias("chain.id.to"))
    result, summary = apply_normalized(table, reversed_edges, {})
    assert summary["normalized"] is None and summary["normalized_missing_geometry"] == 2
    result, summary = apply_normalized(table.clear(), increments.clear(), {})
    assert summary["normalized"] == 0.
    with pytest.raises(ValueError):
        apply_normalized(table, pl.concat([increments, increments]), {})


def test_fit_explicit_native_class_and_directed_v2_dependencies():
    from test_contact_observations import complex_fixture
    from tcren.contacts.observations import contact_observations
    from tcren.contacts.plausibility import fit_resource
    obs = contact_observations(complex_fixture(), interface="tcr_peptide")
    ref = pl.concat([obs, obs.with_columns((pl.col("dist")*.9).alias("dist"))])
    prior = fit_resource(ref, obs)
    model = fit_normalizer(obs, prior, "abc", "I")
    validate_normalizer(model, "abc")
    assert model["reference_counts"]["n_edges"] == 1
    with pytest.raises(ValueError, match="no native"):
        fit_normalizer(obs, prior, "abc", "II")
    with pytest.raises(ValueError, match="v2"):
        fit_normalizer(obs, dict(prior, version="contact-marginals-1"), "abc", "I")


def test_normalizer_cli_export_and_comparison_keep_declared_arms(tmp_path, monkeypatch):
    import json
    from test_contact_observations import complex_fixture
    import tcren.annotation
    import tcren.annotation.batch
    import tcren.mhc
    import tcren.structure
    from tcren.contacts.context_normalization import normalizer_command
    from tcren.contacts.observations import contact_observations
    from tcren.contacts.penalty_export import penalty_export
    from tcren.contacts.plausibility import ARMS, comparison_command, fit_resource
    structure = complex_fixture()
    obs = contact_observations(structure, interface="tcr_peptide")
    ref = pl.concat([obs, obs.with_columns((pl.col("dist")*.9).alias("dist"))])
    resource, native, model = [tmp_path/name for name in ["prior.json", "native.parquet", "normalizer.json"]]
    resource.write_text(json.dumps(fit_resource(ref, obs)))
    obs.write_parquet(native)
    normalizer_command(native, resource, "I", model, tmp_path/"reference_checks.tsv")
    source = tmp_path/"synthetic.pdb"
    source.write_text("synthetic parser input")
    monkeypatch.setattr(tcren.structure, "parse_structure", lambda *a, **k: structure)
    monkeypatch.setattr(tcren.annotation.batch, "annotate_batch", lambda structures, **k: [{} for _ in structures])
    monkeypatch.setattr(tcren.annotation, "classify_chains", lambda *a, **k: None)
    monkeypatch.setattr(tcren.mhc, "annotate_mhc_batch", lambda *a, **k: None)
    output = tmp_path/"scores"
    penalty_export([str(source)], resource, output, normalizer_path=model)
    scored = pl.read_csv(output/"scores.tsv", separator="\t")
    assert set(["original", *ARMS, "normalized"]) <= set(scored.columns)
    assert scored["normalized_missing_geometry"].sum() == 0
    assert json.loads((output/"provenance.json").read_text())["normalizer_sha256"]
    labels = tmp_path/"labels.tsv"
    pl.DataFrame(dict(id=[structure.pdb_id], dataset=["d"], task=["ranking"], cohort=["e"],
        mhc_class=["I"], template_coverage=["covered"], y=[1], iptm=[.8])).write_csv(labels, separator="\t")
    comparison_command(labels, output/"scores.tsv", tmp_path/"comparison")
    metrics = pl.read_csv(tmp_path/"comparison/metrics.tsv", separator="\t")
    assert set(metrics["arm"]) == {"original", *ARMS, "normalized", "iptm"}
    assert metrics["n_excluded"].sum() == 0


def test_missing_potential_entry_does_not_require_geometry():
    reference = edges((.2, None), (2, 2))
    increments = normalized_edges(reference, normalizer(edges()))
    table = reference.select(KEYS).with_columns(
        pl.Series("original_contribution", [-2., None]),
        pl.lit("supported_cdr_geometry").alias("geometry_status"))
    result, summary = apply_normalized(table, increments, {})
    assert summary["normalized_missing_geometry"] == 0
    assert summary["normalized"] == result["normalized_contribution"][0]
    assert result["normalized_contribution"][1] is None
