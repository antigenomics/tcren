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


def test_peptide_native_features_keep_contact_count_from_residue_weighted_pipeline(monkeypatch):
    from tcren import pipeline
    from tcren.native_ranking import native_epitope_features
    from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure
    from tcren.topology import surface

    def chain(id_, role, x):
        ca = Atom("CA", "C", np.asarray([x, 0., 0.]))
        cb = Atom("CB", "C", np.asarray([x, 1., 0.]))
        residues = [Residue(0, 0, "", "A", "ALA", (ca, cb))]
        return Chain(id_, residues, role, regions=[RegionMarkup(role, 0, 0, "A", residues)])

    structure = Structure("pose", [chain("A", "TRA", 0), chain("P", "PEPTIDE", 3)])
    structure.mhc_calls = []
    calls = []

    def run_with_count(s, **kwargs):
        calls.append(kwargs)
        return original_run(s, **kwargs)

    original_run = pipeline.run
    monkeypatch.setattr(pipeline, "run", run_with_count)
    monkeypatch.setattr(surface, "surface_map", lambda s, **kwargs: kwargs["side"])
    monkeypatch.setattr(surface, "surface_complementarity",
                        lambda *_: {"d_charge": .3, "coverage": .6})
    row = native_epitope_features([structure], threads=1).row(0, named=True)
    assert calls[0]["typed"] and calls[0]["reference_mode"] == "structural"
    assert row["complex.id"] == "pose"
    assert row["n_atom_contacts_tp"] == 4, row.get("native_energy_error")
    assert row["dPhi_peptide_sum"] == pytest.approx(0)
    assert row["d_charge"] == pytest.approx(.3)
    assert row["coverage"] == pytest.approx(.6)


def test_native_screen_from_peptide_structures_reads_confidence_and_types_once(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from tcren import cli, native_ranking, screen

    base = _inputs()
    structures = [SimpleNamespace(pdb_id=item) for item in base["complex.id"]]
    batches = []

    def typed(*args, **kwargs):
        batches.append(kwargs)
        return structures

    monkeypatch.setattr(cli, "_typed_batch", typed)
    monkeypatch.setattr(screen, "identities",
                        lambda batch: (list("PQRS"), [("A", "B")] * len(batch)))
    monkeypatch.setattr(native_ranking, "native_epitope_features",
                        lambda *_args, **_kwargs: base.drop("iptm", "unit"))
    confidence = tmp_path / "confidence.tsv"
    base.select("complex.id", "iptm").write_csv(confidence, separator="\t")
    out = tmp_path / "scores.tsv"
    result = CliRunner().invoke(app, [
        "screen", "--method", "native387-v1", "-s", "models", "--task", "epitope",
        "--confidence", str(confidence), "--min-cohort", "1", "-o", str(out), "-t", "1",
    ])
    assert result.exit_code == 0, result.output
    assert batches == [{"threads": 1}]
    expected = rank_native_screen(base, task="epitope", ranking_unit="unit", min_cohort=1)
    actual = pl.read_csv(out, separator="\t")
    np.testing.assert_allclose(actual["screen_score"], expected["screen_score"])
