"""Confidence proxy fitting must not learn from a held-out epitope."""

import numpy as np
import polars as pl
import pytest
from tcren.contacts.confidence_proxy import allele_group, proxy_predictions


def test_proxy_holdout_and_allele_resolution():
    rows = pl.DataFrame(
        {
            "id": ["a", "b", "c", "d", "e", "f"],
            "dataset": ["x"] * 6,
            "stratum": ["x"] * 6,
            "epitope": ["a", "a", "b", "b", "c", "c"],
            "iptm_scored": [0.6, 0.7, 0.8, 0.9, 0.4, 0.5],
            "plddt": [91.0, 93.0, 95.0, 94.0, 90.0, 89.0],
            "ranking_confidence": [0.7, None, 0.8, 0.9, 0.5, 0.4],
            "betti_scored": [0.2] * 6,
            "enrichment_scored": [0.4] * 6,
            "mhc": ["HLA-A*02:01"] * 6,
        }
    )
    _, target, before, _ = proxy_predictions(rows)
    changed = rows.with_columns(
        pl.when(pl.col("id") == "b").then(80.0).otherwise(pl.col("plddt")).alias("plddt")
    )
    _, target2, after, _ = proxy_predictions(changed)
    assert target[0, 0] != target2[0, 0]
    # Row a's features stay fixed while its target cohort changes: no target leakage.
    for name in [
        "linear_confidence",
        "linear_structure",
        "iptm_decile",
        "mhc_field1",
        "mhc_field2",
    ]:
        np.testing.assert_allclose(before[name][0], after[name][0], atol=1e-12)
    assert allele_group(["HLA-A*02:01:01", "B*07:02", "unknown"], 1).tolist() == [
        "A*02",
        "B*07",
        None,
    ]
    assert allele_group(["HLA-A*02:01:01"], 2).tolist() == ["A*02:01"]
    with pytest.raises(ValueError):
        allele_group(["HLA-A*02:01"], 3)


def test_allele_proxy_fallback_and_frozen_cli(tmp_path):
    from tcren.contacts.confidence_proxy import allele_confidence, score_command

    context = allele_confidence(["HLA-A*02:01", "HLA-A*99:99"])
    assert context["confidence_proxy_fallback"].to_list() == [False, True]
    assert context["n_models"].to_list() == [2101, 3056]
    source, output = tmp_path / "features.tsv", tmp_path / "scores.tsv"
    pl.DataFrame(
        {
            "mhc": ["HLA-A*02:01", "HLA-A*99:99"],
            "iptm": [0.8, 0.7],
            "betti": [0.3, 0.4],
            "enrichment": [0.6, 0.5],
        }
    ).write_csv(source, separator="\t")
    score_command(source, output, 2, None)
    result = pl.read_csv(output, separator="\t")
    np.testing.assert_allclose(result["score"], [0.6272, 0.6086])
    assert result["confidence_context"].unique().to_list() == ["allele_proxy_full_reference"]


def test_proxy_invalid_reference_and_absent_components(tmp_path):
    import gzip
    from importlib.resources import files
    from tcren.contacts.confidence_proxy import allele_confidence, score_command

    reference = pl.read_csv(
        gzip.decompress(files("tcren.data").joinpath("confidence_alleles.tsv.gz").read_bytes()),
        separator="\t",
    )
    for invalid in [float("nan"), float("inf"), -1.0, 101.0]:
        with pytest.raises(ValueError):
            allele_confidence(
                ["A*02:01"], reference=reference.with_columns(pl.lit(invalid).alias("mean_plddt"))
            )
    source, output = tmp_path / "missing.tsv", tmp_path / "scores.tsv"
    pl.DataFrame(
        {"mhc": ["A*02:01"], "iptm": [None], "betti": [None], "enrichment": [None]}
    ).write_csv(source, separator="\t")
    score_command(source, output, 2, None)
    result = pl.read_csv(output, separator="\t")
    assert result["score"].null_count() == 1
    assert result.select("missing_iptm", "missing_betti", "missing_enrichment").row(0) == (
        True,
        True,
        True,
    )
