"""Matched populations and independent strata are scientific invariants."""
import numpy as np
import polars as pl
import pytest

from tcren.associations import associate


def test_matched_metrics_strata_and_missingness():
    rows = []
    features = []
    for coverage in ("covered", "scarce"):
        for cohort in range(3):
            for i in range(4):
                key = f"{coverage}-{cohort}-{i}"
                rows.append(dict(id=key, dataset="example", task="receptor", cohort=str(cohort),
                                 mhc_class="I", template_coverage=coverage, y=i // 2,
                                 iptm=float(i)/3, score=float(3-i) if i != 0 else None))
                features.append(dict(id=key, relief=float(cohort)))
    result = associate(pl.DataFrame(rows), pl.DataFrame(features), {"score": 1}, ["relief"])
    metrics = result["metrics"]
    auc = metrics.filter(pl.col("metric") == "roc_auc")
    assert auc.height == 6
    assert auc["n_used"].to_list() == [3] * 6
    assert auc["n_excluded"].to_list() == [1] * 6
    assert auc["n_positive"].to_list() == [2] * 6
    assert auc["value"].to_list() == [0.] * 6
    assert auc["iptm_value"].to_list() == [1.] * 6
    assert auc["delta"].to_list() == [-1.] * 6
    assert set(result["associations"]["template_coverage"]) == {"covered", "scarce"}
    assert set(result["associations"]["status"]) == {"constant"}
    with pytest.raises(ValueError, match="duplicate"):
        associate(pl.DataFrame(rows), pl.concat([pl.DataFrame(features)] * 2), {"score": 1}, ["relief"])


def test_ties_invalid_labels_and_missing_feature_accounting():
    rows = pl.DataFrame(dict(id=["a", "b", "c", "d"], dataset=["d"]*4,
        task=["receptor"]*4, cohort=["p"]*4, mhc_class=["II"]*4,
        template_coverage=["unknown"]*4, y=[0, 0, 1, 1], iptm=[.1,.2,.3,.4], score=[1.]*4))
    features = pl.DataFrame(dict(id=["a", "b"], relief=[2., 2.]))
    result = associate(rows, features, {"score": -1}, ["relief"])
    assert result["metrics"].filter(pl.col("metric") == "roc_auc")["value"][0] == .5
    assert result["features"]["n_feature"][0] == 2
    assert result["features"]["n_missing"][0] == 2
    assert result["associations"]["n_cohorts"][0] == 1
    with pytest.raises(ValueError, match="binary"):
        associate(rows.with_columns(pl.lit(2).alias("y")), features, {"score": 1}, ["relief"])
    with pytest.raises(ValueError, match="direction"):
        associate(rows, features, {"score": 0}, ["relief"])


def test_pmhc_features_do_not_depend_on_receptor():
    from test_exports import typed_structure
    from tcren.topology.pmhc import pmhc_features
    from tcren.structure.model import Structure

    s = typed_structure()
    expected = pmhc_features(s)
    stripped = Structure(s.pdb_id, [c for c in s.chains if c.chain_type not in ("TRA", "TRB")])
    actual = pmhc_features(stripped)
    assert actual.keys() == expected.keys()
    for k in actual:
        if isinstance(actual[k], float):
            assert actual[k] == pytest.approx(expected[k], nan_ok=True)
        else:
            assert actual[k] == expected[k]
    assert actual["peptide_sasa_A2"] > actual["sidechain_sasa_A2"] > 0
    assert np.isfinite(actual["sasa_hydropathy"])
    # Explicit hydrogens must not alter the heavy-atom exposure convention.
    from dataclasses import replace
    from tcren.structure.model import Atom

    for chain in stripped.chains:
        residue = chain.residues[0]
        chain.residues[0] = replace(residue, atoms=residue.atoms + (
            Atom("H", "H", residue.atoms[0].coord + np.array([0., 0., 20.])),))
    hydrogenated = pmhc_features(stripped)
    for key in ("peptide_sasa_A2", "sidechain_sasa_A2", "sasa_hydropathy"):
        assert hydrogenated[key] == pytest.approx(actual[key])


def test_association_cli_writes_provenance_and_refuses_overwrite(tmp_path):
    import json
    from typer.testing import CliRunner
    from tcren.cli import app

    table = pl.DataFrame(dict(id=["a", "b"], dataset=["d"]*2, task=["receptor"]*2,
        cohort=["p"]*2, mhc_class=["II"]*2, template_coverage=["unknown"]*2,
        y=[0, 1], iptm=[.2, .8], score=[2., 1.]))
    table.write_csv(tmp_path / "results.csv")
    pl.DataFrame(dict(id=["a", "b"], relief=[1., 1.])).write_csv(tmp_path / "features.csv")
    (tmp_path / "scores.json").write_text(json.dumps({"score": -1}))
    args = ["associate", "--results", str(tmp_path / "results.csv"), "--features",
            str(tmp_path / "features.csv"), "--scores", str(tmp_path / "scores.json"),
            "--feature-columns", "relief", "-o", str(tmp_path / "out")]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    assert (tmp_path / "out/provenance.json").is_file()
    assert CliRunner().invoke(app, args).exit_code != 0
