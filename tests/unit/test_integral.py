"""Fixed-reference integration must not learn from the query batch or its labels."""
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from tcren.integral import build_reference, integral_scores


def table():
    return pl.DataFrame({"complex.id": ["a", "b", "c"], "mhc_class_bin": [0., 0., 0.],
        "Phi_tcr_pep": [2., 0., -2.], "n_contacts_tp": [2., 2., 2.],
        "sc_shape": [-1., 0., 1.], "g_cyclo_frac": [0., .5, 1.],
        "iptm": [0., .5, 1.]})


def test_fixed_reference_midrank_and_batch_invariance():
    data = table()
    ref = build_reference(data, provenance={"input_sha256": "a" * 64})
    out = integral_scores(data, ref)
    np.testing.assert_allclose(out["integral_score"], [.125, .5, .875])
    assert out["n_components"].to_list() == [4, 4, 4]
    one = integral_scores(data.slice(1, 1), ref)
    assert one["integral_score"][0] == out["integral_score"][1]
    assert integral_scores(data.with_columns(pl.lit(1).alias("y")), ref).equals(out)
    assert out["complex.id"].to_list() == data["complex.id"].to_list()


def test_missing_is_flagged_and_neutral_not_dropped():
    ref = build_reference(table(), provenance={"input_sha256": "a" * 64})
    data = table().head(1).with_columns(pl.lit(None, dtype=pl.Float64).alias("iptm"))
    out = integral_scores(data, ref)
    assert out["missing_iptm"][0] and out["n_components"][0] == 3
    assert out["integral_score"][0] == pytest.approx(.25)
    data = data.with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c)
                              for c in ("Phi_tcr_pep", "sc_shape", "g_cyclo_frac")])
    assert integral_scores(data, ref)["integral_score"][0] is None


@pytest.mark.parametrize("column,value", [("iptm", 1.1), ("sc_shape", -1.1),
    ("g_cyclo_frac", -.1), ("n_contacts_tp", -1.), ("n_contacts_tp", 1.5),
    ("Phi_tcr_pep", np.inf), ("mhc_class_bin", 1.)])
def test_invalid_domain_fails(column, value):
    ref = build_reference(table(), provenance={"input_sha256": "a" * 64})
    with pytest.raises(ValueError):
        integral_scores(table().with_columns(pl.lit(value).alias(column)), ref)


def test_zero_contacts_and_reference_validation():
    ref = build_reference(table(), provenance={"input_sha256": "a" * 64})
    data = table().head(1).with_columns(pl.lit(0.).alias("n_contacts_tp"))
    with pytest.raises(ValueError, match="zero contacts"):
        integral_scores(data, ref)
    data = data.with_columns(pl.lit(0.).alias("Phi_tcr_pep"))
    assert integral_scores(data, ref)["missing_energy"][0]
    ref["potential_sha256"] = "b" * 64
    with pytest.raises(ValueError, match="potential"):
        integral_scores(data, ref)


def test_reference_excludes_class_ii_and_rejects_duplicate_ids():
    data = table().with_columns(pl.Series("mhc_class_bin", [0., 0., 1.]))
    ref = build_reference(data, provenance={"input_sha256": "a" * 64})
    assert ref["n_reference"] == 2 and ref["n_excluded_class_ii"] == 1
    with pytest.raises(ValueError, match="unique"):
        build_reference(pl.concat([table(), table()]), provenance={"input_sha256": "a" * 64})


def test_cli_roundtrip_keyed_confidence_and_provenance(tmp_path):
    from typer.testing import CliRunner
    from tcren.cli import app
    runner = CliRunner()
    data, ref, conf, out = [tmp_path / p for p in ("features.tsv", "ref.json", "conf.tsv", "out.tsv")]
    table().drop("iptm").write_csv(data, separator="\t")
    table().select("complex.id", "iptm").reverse().write_csv(conf, separator="\t")
    made = runner.invoke(app, ["integral-reference", "--features", str(data), "-o", str(ref)])
    assert made.exit_code == 0, made.output
    result = runner.invoke(app, ["integrate", "--features", str(data), "--reference", str(ref),
                                "--confidence", str(conf), "-o", str(out)])
    assert result.exit_code == 0, result.output
    np.testing.assert_allclose(pl.read_csv(out, separator="\t")["integral_score"], [.125, .5, .875])
    assert Path(str(out) + ".provenance.json").exists()
