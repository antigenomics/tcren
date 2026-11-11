"""Frozen marginal updating and the explicit finite penalty menu."""
import numpy as np
import polars as pl
from scipy.stats import t


def test_predictive_update_and_penalties():
    from tcren.contacts.plausibility import combine_stats, predictive_compatibility, penalty_contributions
    a = {"n": 3, "mean": 2., "m2": 2.}
    b = {"n": 1, "mean": 4., "m2": 0.}
    merged = combine_stats(a, b)
    assert merged == {"n": 4, "mean": 2.5, "m2": 5.}
    x = np.array([2.5, 4.])
    expected = 2*t.sf(np.abs(x-2.5)/np.sqrt((5/3)*(1+1/4)), 3)
    np.testing.assert_allclose(predictive_compatibility(x, merged), expected)
    e, u = np.array([-2., 3., 0.]), np.array([1., .5, 0.])
    out = penalty_contributions(e, u, 3.)
    np.testing.assert_array_equal(out["fixed"], [-1., 3.5, 0.])
    np.testing.assert_array_equal(out["matrix"], [1., 4.5, 0.])
    np.testing.assert_array_equal(out["flip"], [2., 0., 0.])
    np.testing.assert_array_equal(out["flip_favourable"], [2., 3., 0.])
    for name in ("fixed", "matrix", "flip_favourable"):
        assert np.all(out[name] >= e)
    assert np.isnan(penalty_contributions(np.array([np.nan]), np.array([1.]), 3.)["fixed"]).all()


def test_marginals_follow_biological_direction_and_missing_angles_stay_missing():
    from tcren.contacts.plausibility import marginal_samples
    table = pl.DataFrame({"pdb.id": ["s"], "chain.id.from": ["A"], "chain.id.to": ["B"],
        "residue.index.from": [0], "residue.index.to": [0], "chem.from": ["A:CA"], "chem.to": ["G:N"],
        "role.from": ["tcr"], "role.to": ["peptide"],
        "dist": [3.], "omega_rad": [None], "theta.from_rad": [.3], "theta.to_rad": [.7],
        "phi.from_rad": [.2], "phi.to_rad": [None]}, infer_schema_length=None)
    reverse = table.rename({c: c.replace(".from", ".swap").replace(".to", ".from").replace(".swap", ".to") for c in table.columns})
    def values(t):
        return sorted((r["key"], r["pool"], r["x"]) for r in marginal_samples(t).to_dicts())
    assert values(table) == values(reverse)
    assert len(values(table)) == 4
    keys = {key for key, _, _ in values(table)}
    assert keys == {"distance|A:CA|G:N", "theta|A|G|binder", "theta|A|G|target", "phi|A|G|binder"}
    # Biological exchange changes residue keys, even for the same physical atoms.
    exchanged = table.with_columns(pl.lit("peptide").alias("role.from"), pl.lit("tcr").alias("role.to"))
    changed = values(exchanged)
    assert [v for v in values(table) if v[0].startswith("distance|")] == [v for v in changed if v[0].startswith("distance|")]
    assert {k for k, _, _ in changed if not k.startswith("distance|")} == {
        "theta|G|A|binder", "theta|G|A|target", "phi|G|A|target"}
    omega = table.with_columns(pl.lit(.4).alias("omega_rad"))
    swapped_omega = exchanged.with_columns(pl.lit(.4).alias("omega_rad"))
    assert "omega|A|G" in {k for k, _, _ in values(omega)}
    assert "omega|G|A" in {k for k, _, _ in values(swapped_omega)}
    # Same-amino-acid endpoints must remain distinct in angular backoff, too.
    equal = marginal_samples(table.with_columns(pl.lit("A:N").alias("chem.to")))
    theta = equal.filter(pl.col("key").str.starts_with("theta|"))
    assert set(theta["key"]) == {"theta|A|A|binder", "theta|A|A|target"}
    assert set(theta["pool"]) == {"theta|binder", "theta|target"}
    import pytest
    with pytest.raises(ValueError, match="requires role"):
        marginal_samples(table.drop("role.from", "role.to"))
    with pytest.raises(ValueError, match="biological binder/target"):
        marginal_samples(table.with_columns(pl.lit(None).cast(pl.String).alias("role.from")))
    with pytest.raises(ValueError, match="interface disagrees"):
        marginal_samples(table.with_columns(pl.lit("peptide_mhc").alias("interface")))
    # The external deposit schema explicitly declares its physical from side binder.
    external = table.drop("role.from", "role.to").rename({
        "chain.id.from": "source_chain.from", "residue.index.from": "source_resid.from",
        "chain.id.to": "source_chain.to", "residue.index.to": "source_resid.to"})
    assert values(external) == values(table)
    assert values(external.rename({c: c.replace(".from", ".swap").replace(".to", ".from").replace(".swap", ".to") for c in external.columns})) != values(table)
    for binder, target in (("tcr", "mhc"), ("peptide", "mhc")):
        assert values(table.with_columns(pl.lit(binder).alias("role.from"), pl.lit(target).alias("role.to"))) == values(table)
    # Recover a supported orientation even if the first contact atom has no frame.
    later = table.with_columns(pl.lit(.4).alias("omega_rad"))
    mixed = pl.concat([table, later], how="vertical_relaxed")
    def fitted(t):
        return sorted((r["key"], r["pool"], r["x"]) for r in marginal_samples(t, fitting=True).to_dicts())
    assert fitted(mixed) == fitted(mixed.reverse())


def test_residue_penalty_conservation_and_unchanged_baseline():
    from test_contact_observations import complex_fixture
    from tcren.contacts.observations import contact_observations
    from tcren.contacts.plausibility import fit_resource, score_contact_penalties
    from tcren.contactmap import ContactMap
    from tcren.potential import tcren2
    structure = complex_fixture()
    observations = contact_observations(structure, interface="tcr_peptide")
    # Independent reference coordinates with several distances give nonzero variance.
    reference = pl.concat([observations, observations.with_columns(
        (pl.col("dist")*.95).alias("dist"), pl.lit("second").alias("pdb.id"))])
    resource = fit_resource(reference, observations)
    import pytest
    from tcren.contacts.plausibility import atom_weights
    assert resource["version"] == "contact-marginals-2"
    with pytest.raises(ValueError, match="refit historical v1"):
        atom_weights(observations, dict(resource, version="contact-marginals-1"))
    before = repr(resource)
    contacts = ContactMap.from_structure(structure).interface("tcr_peptide")
    edges, scores = score_contact_penalties(contacts, observations, tcren2(), resource)
    assert repr(resource) == before
    assert edges.height == 1
    for name in ("fixed", "matrix", "flip", "flip_favourable"):
        assert scores[name] == edges[name+"_contribution"].sum()
        np.testing.assert_allclose(edges["original_contribution"].to_numpy()+edges[name+"_increment"].to_numpy(),
                                   edges[name+"_contribution"].to_numpy())
    assert scores["fixed"] >= scores["original"]
    empty_edges, empty_scores = score_contact_penalties(contacts.clear(), observations.clear(), tcren2(), resource)
    assert empty_edges.height == 0
    assert all(empty_scores[arm] == 0 for arm in ("original", "fixed", "matrix", "flip", "flip_favourable"))
    # Reversing physical endpoints preserves the weight and every resulting sum.
    reverse = observations.rename({c: c.replace(".from", ".swap").replace(".to", ".from").replace(".swap", ".to") for c in observations.columns})
    other, other_scores = score_contact_penalties(contacts, reverse, tcren2(), resource)
    np.testing.assert_allclose(other["unnaturalness"].to_numpy(), edges["unnaturalness"].to_numpy())
    assert scores == other_scores


def test_cli_comparison_uses_identical_rows_for_every_arm(tmp_path):
    from tcren.contacts.plausibility import ARMS, comparison_command
    labels = pl.DataFrame(dict(id=["a", "b", "c"], dataset=["d"]*3, task=["ranking"]*3,
        cohort=["e"]*3, mhc_class=["I"]*3, template_coverage=["scarce"]*3, y=[0, 1, 1], iptm=[.1, .9, .8]))
    scores = pl.DataFrame(dict(id=["a", "b", "c"], potential=["p"]*3, mhc_class=["I"]*3,
        original=[1., -1., -2.], missing_potential_entries=[0]*3,
        **{arm: ([1., -1., None] if arm=="fixed" else [1., -1., -2.]) for arm in ARMS}))
    lp, sp, output = tmp_path/"labels.tsv", tmp_path/"scores.tsv", tmp_path/"out"
    labels.write_csv(lp, separator="\t")
    scores.write_csv(sp, separator="\t")
    comparison_command(lp, sp, output)
    metrics = pl.read_csv(output/"metrics.tsv", separator="\t")
    assert metrics["n_used"].unique().to_list() == [2]
    assert metrics["n_excluded"].unique().to_list() == [1]
    assert set(metrics["arm"]) == {"original", "iptm", *ARMS}
    surface_scores = scores.with_columns(pl.lit(0.).alias("fixed"),
        pl.lit(0.).alias("surface_raw_sum"), pl.lit(0.).alias("surface_context_sum"),
        pl.lit(None).cast(pl.Float64).alias("surface_context_mean"))
    surface_scores.write_csv(sp, separator="\t")
    comparison_command(lp, sp, tmp_path/"surface_compare")
    surface_metrics = pl.read_csv(tmp_path/"surface_compare/metrics.tsv", separator="\t")
    assert surface_metrics["n_used"].unique().to_list() == [3]
    assert {"surface_raw_sum", "surface_context_sum"} <= set(surface_metrics["arm"])
    chemical_scores = surface_scores.with_columns(
        pl.Series("chemical_core", [1., -1., None]),
        pl.Series("atom_count", [1., -1., -2.]))
    chemical_scores.write_csv(sp, separator="\t")
    comparison_command(lp, sp, tmp_path/"chemical_compare")
    chemical_metrics = pl.read_csv(tmp_path/"chemical_compare/metrics.tsv", separator="\t")
    assert chemical_metrics["n_used"].unique().to_list() == [2]
    assert {"chemical_core", "atom_count", "original"} <= set(chemical_metrics["arm"])
    import pytest
    scores.with_columns(pl.lit("II").alias("mhc_class")).write_csv(sp, separator="\t")
    with pytest.raises(ValueError, match="MHC class"):
        comparison_command(lp, sp, tmp_path/"mismatch")


def test_bounded_export_writes_complete_readable_artifacts(tmp_path, monkeypatch):
    import gzip
    import json
    import pytest
    from test_contact_observations import complex_fixture
    import tcren.annotation
    import tcren.annotation.batch
    import tcren.mhc
    import tcren.structure
    from tcren.contacts.observations import contact_observations
    from tcren.contacts.plausibility import fit_resource
    from tcren.contacts.penalty_export import penalty_export
    structure = complex_fixture()
    observations = contact_observations(structure, interface="tcr_peptide")
    reference = pl.concat([observations, observations.with_columns((pl.col("dist")*.95).alias("dist"))])
    resource = tmp_path/"resource.json"
    resource.write_text(json.dumps(fit_resource(reference, observations)))
    source = tmp_path/"synthetic.pdb"
    source.write_text("synthetic parser input")
    monkeypatch.setattr(tcren.structure, "parse_structure", lambda *a, **k: structure)
    monkeypatch.setattr(tcren.annotation.batch, "annotate_batch", lambda structures, **k: [{} for _ in structures])
    monkeypatch.setattr(tcren.annotation, "classify_chains", lambda *a, **k: None)
    monkeypatch.setattr(tcren.mhc, "annotate_mhc_batch", lambda *a, **k: None)
    output = tmp_path/"out"
    penalty_export([str(source)], resource, output, chunk=2)
    assert json.loads((output/"provenance.json").read_text())["no_scoring_time_fitting"]
    assert pl.read_csv(output/"scores.tsv", separator="\t").height == 2
    from tcren.contacts.chemical_core import derive_chemical_core
    import tcren.contacts.chemical_core as core
    counts = pl.DataFrame(dict(aa_from=["Y"], atom_from=["OH"],
                               aa_to=["P"], atom_to=["O"], count=[1]))
    frozen_core = derive_chemical_core(counts, counts.head(0))
    chemical_resource = tmp_path/"chemical.json"
    chemical_resource.write_text(json.dumps(frozen_core))
    calls = []
    def chemical_score(contacts, atoms, frozen, total, *, mhc_class):
        assert frozen == frozen_core
        assert mhc_class == "I" and atoms.equals(observations)
        calls.append(total)
        return contacts.with_columns(pl.lit(.25).alias("chemical_core_increment")), {
            "chemical_core": total+.25, "chemical_core_increment": .25}
    monkeypatch.setattr(core, "score_chemical_core", chemical_score)
    chemical_output = tmp_path/"chemical_out"
    penalty_export([str(source)], resource, chemical_output, chunk=2,
                   threads=1, chemical_core_path=chemical_resource)
    before = pl.read_csv(output/"scores.tsv", separator="\t")
    chemical = pl.read_csv(chemical_output/"scores.tsv", separator="\t")
    assert chemical.select(before.columns).equals(before)
    assert calls == before["original"].to_list()
    assert chemical["chemical_core_increment"].to_list() == [.25, .25]
    from tcren.exports import _sha256
    assert json.loads((chemical_output/"provenance.json").read_text())["chemical_core_sha256"] == _sha256(chemical_resource)
    with pytest.raises(ValueError, match="requires score export"):
        penalty_export([str(source)], resource, tmp_path/"chemical_geometry",
                       chemical_core_path=chemical_resource, geometry_only=True)
    for value in ({}, {"version": "malformed"}):
        chemical_resource.write_text(json.dumps(value))
        invalid_output = tmp_path/"invalid_chemical"
        with pytest.raises(ValueError, match="chemical-core"):
            penalty_export([str(source)], resource, invalid_output,
                           chemical_core_path=chemical_resource)
        assert not invalid_output.exists()
    surface_output = tmp_path/"surface"
    penalty_export([str(source)], resource, surface_output, chunk=2, surface_context=True)
    original = pl.read_csv(output/"scores.tsv", separator="\t")
    added = pl.read_csv(surface_output/"scores.tsv", separator="\t")
    assert original.equals(added.select(original.columns))
    assert added["surface_context_sum"].is_finite().all()
    assert pl.read_parquet(surface_output/"surface/atoms-000000.parquet").height == 10
    boundary_output = tmp_path/"boundary"
    penalty_export([str(source)], resource, boundary_output, chunk=2,
                   surface_context=True, surface_ownership="accessible_boundary")
    boundary = pl.read_csv(boundary_output/"scores.tsv", separator="\t")
    assert original.equals(boundary.select(original.columns))
    assert pl.read_parquet(boundary_output/"surface/atoms-000000.parquet").equals(
        pl.read_parquet(surface_output/"surface/atoms-000000.parquet"))
    edges = pl.read_parquet(boundary_output/"surface/edges-000000.parquet")
    assert edges.filter((pl.col("binder_area_A2")>0) & (pl.col("isolated_sasa_A2.to")<=0)).is_empty()
    assert json.loads((boundary_output/"provenance.json").read_text())["surface_context"]["ownership"] == "accessible_boundary"
    diagnostic_output = tmp_path/"diagnostics"
    penalty_export([str(source)], resource, diagnostic_output, chunk=2, geometry_only=True)
    assert not (diagnostic_output/"scores.tsv").exists()
    assert pl.read_csv(diagnostic_output/"geometry.tsv.gz", separator="\t").height == 1
    assert pl.read_csv(diagnostic_output/"geometry_models.tsv", separator="\t").height == 1
    with gzip.open(output/"contacts.tsv.gz", "rt") as handle:
        assert len(handle.read().splitlines()) == 3


def test_geometry_diagnostics_reproduce_weights_and_keep_all_atom_chemistry():
    from test_contact_observations import complex_fixture
    from tcren.contacts.observations import contact_observations
    from tcren.contacts.plausibility import fit_resource, atom_weights
    from tcren.contacts.mechanism import geometry_diagnostics
    obs = contact_observations(complex_fixture(), interface="tcr_peptide")
    ref = pl.concat([obs, obs.with_columns((pl.col("dist")*.95).alias("dist"))])
    resource = fit_resource(ref, obs)
    result = geometry_diagnostics(obs, resource)
    assert result.height == 1
    np.testing.assert_allclose(result["audit_unnaturalness"].item(), atom_weights(obs, resource)["unnaturalness"].max())
    assert result["n_atom_pairs"].item() == obs.height
    assert result["n_distinct_distances"].item() <= obs.height
    assert result["winner_family"].item() in {"distance", "omega", "theta", "phi"}
    assert result["chemical_supported"].item() == bool(obs.select(pl.any_horizontal([pl.col("is_"+t) for t in ("salt_bridge", "hydrogen_bond", "cation_pi", "stacking", "aromatic", "hydrophobic", "polar")]).any()).item())
    assert geometry_diagnostics(obs.clear(), resource).height == 0

    reverse = obs.rename({c: c.replace(".from", ".swap").replace(".to", ".from").replace(".swap", ".to") for c in obs.columns})
    reversed_result = geometry_diagnostics(reverse, resource)
    for column in ["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to", "winner_chem.from", "winner_chem.to"]:
        assert result[column].to_list() == reversed_result[column].to_list()
