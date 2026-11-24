"""The two measured cohort formulas and zero-contact screen policy."""

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from tcren.cli import app
from tcren.native_ranking import rank_native_screen


def _z(values):
    values = np.asarray(values)
    return (values - values.mean()) / values.std()


def _inputs():
    return pl.DataFrame({
        "complex.id": ["a", "b", "c", "d"],
        "unit": ["clone"] * 4,
        "iptm": [.2, .45, .65, .8],
        "dPhi_peptide_sum": [5., 1., -2., 0.],
        "d_charge": [.4, .8, .6, .1],
        "coverage": [.2, .6, .4, 1.],
        "n_atom_contacts_tp": [0, 5, 15, 20],
    })


def test_peptide_score_keeps_recorded_weights_and_zero_contact_reason():
    t = _inputs()
    p = t["iptm"].to_numpy()
    base = _z(np.log(p / (1 - p))) / 2 - _z(t["dPhi_peptide_sum"]) / 2
    expected = .8 * _z(base) - .1 * _z(t["d_charge"]) + .1 * _z(t["coverage"])
    r = rank_native_screen(t, task="epitope", ranking_unit="unit", min_cohort=1)
    np.testing.assert_allclose(r["workflow_score"], expected)
    assert r["screen_score"][0] == 0
    assert r["no_interface"].to_list() == [True, False, False, False]
    assert r["score_reason"][0] == "no_peptide_contact"
    np.testing.assert_allclose(r["screen_score"][1:], expected[1:])
    np.testing.assert_allclose(r["term_confidence_energy"], .8 * _z(base))


@pytest.mark.parametrize("covered", [True, False])
def test_receptor_terms_use_the_correct_native_coverage_branch(covered):
    from tcren.native_reference import load_native_reference

    ref = load_native_reference()
    peptide = next(iter(ref["energy_by_peptide"])) if covered else "NOT_IN_THE_REFERENCE"
    t = _inputs().with_columns(
        peptide=pl.lit(peptide),
        plddt=pl.Series([80., None, 90., 95.]),
        contact_language=pl.Series([None] * 4, dtype=pl.Float64),
        pi_atom=pl.Series([4., 1., 3., 2.]),
        interlock_frac=pl.Series([.5, .4, .7, .8]),
        kept_edge_frac=pl.Series([.3, .2, .8, .5]),
        n_atom_contacts_tp=pl.Series([10, 165, 175, 200]),
    )
    p = t["iptm"].to_numpy()
    terms = [_z(np.log(p / (1 - p))), _z([80., 90., 90., 95.]), _z(t["pi_atom"])]
    if covered:
        anchor = ref["energy_by_peptide"][peptide]["centre"]
        terms += [_z(t["interlock_frac"]), -_z(np.abs(t["dPhi_peptide_sum"].to_numpy() - anchor))]
    else:
        terms += [_z(t["kept_edge_frac"])]
    terms += [_z(-np.abs(t["n_atom_contacts_tp"].to_numpy() - ref["size"]["centre"]))]
    r = rank_native_screen(t, task="tcr", ranking_unit="unit", min_cohort=1)
    assert r["native_reference_available"].to_list() == [covered] * 4
    assert r["language_source"].unique().to_list() == ["pi_atom"]
    np.testing.assert_allclose(r["workflow_score"], np.mean(terms, axis=0))


def test_incomplete_units_keep_rows_and_isolate_other_units():
    t = _inputs()
    bad = t.with_columns(unit=pl.lit("bad"), coverage=pl.Series([.3, None, .6, .1]))
    r = rank_native_screen(pl.concat([bad, t]), task="epitope", ranking_unit="unit", min_cohort=1)
    assert r.height == 8
    assert r["screen_score"].null_count() == 3
    assert r["n_missing_terms"].to_list() == [0, 1, 0, 0, 0, 0, 0, 0]
    assert r["workflow_score"].null_count() == 4
    assert r["screen_score"][0] == 0


def test_native_screen_cli_accepts_explicit_table_and_unit(tmp_path):
    features = tmp_path / "features.tsv"
    out = tmp_path / "scores.tsv"
    _inputs().write_csv(features, separator="\t")
    runner = CliRunner()
    args = ["screen", "--method", "native387-v1", "--features", str(features),
            "--task", "epitope", "--min-cohort", "1", "-o", str(out)]
    assert runner.invoke(app, args).exit_code == 2
    result = runner.invoke(app, [*args, "--ranking-unit", "unit"])
    assert result.exit_code == 0, result.output
    scored = pl.read_csv(out, separator="\t")
    assert scored["screen_score"][0] == 0
    assert scored.height == 4
