"""Production confidence-motif scoring and frozen resource contracts."""
import gzip
import hashlib
import json
from importlib.resources import files

import numpy as np
import polars as pl
import pytest

from tcren.contacts.motif import bundled_enrichment
from tcren.contacts.motif_optimize import apply_confidence_score
from tcren.screen import CohortKind, bundled_reference, confidence_weight, screen_scores


def table(plddt=(92.612607434279, 92.612607434279)):
    return pl.DataFrame({"complex.id": ["a", "b"], "iptm": [.8, .7],
                         "g_recognition_cyclo_frac": [.4, .3], "enrichment": [.6, .5],
                         "plddt": plddt, "dFT_tcr_pep_cdr": [0., 0.],
                         "sc_shape": [.5, .5], "n_contacts_tp": [10, 10]})


def test_selected_formula_and_inclusive_confidence_ties():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    out = screen_scores(table(), ref, min_cohort=1)
    assert out["gate_q"].to_list() == [.4, .4]
    with pytest.raises(ValueError, match="cohort pLDDT"):
        confidence_weight(ref)
    assert confidence_weight(ref, .3, native_template=False) == .3
    expected = .704 * np.array([.8, .7]) + .164 * np.array([.4, .3]) + .132 * np.array([.6, .5])
    np.testing.assert_allclose(out["screen_score"], expected, atol=1e-15)
    assert screen_scores(table(), ref, min_cohort=1, native_template=False)["screen_score"].equals(out["screen_score"])
    np.testing.assert_allclose(screen_scores(table(), ref, min_cohort=1, weight=1.)["screen_score"], [.8,.7])
    assert screen_scores(table(), ref)["screen_score"].null_count() == 2
    absent = table().with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c["column"])
                                   for c in ref["components"]])
    assert screen_scores(absent, ref, min_cohort=1)["screen_score"].null_count() == 2
    absent_active = table().with_columns(*[pl.lit(None).alias(c)
                            for c in ["iptm", "g_recognition_cyclo_frac", "enrichment"]])
    assert screen_scores(absent_active, ref, min_cohort=1)["screen_score"].null_count() == 2
    discarded = table().with_columns(pl.Series("n_contacts_tp", [0, 10]))
    assert screen_scores(discarded, ref, min_cohort=1)["screen_score"][0] is None


def test_cohort_context_missingness_and_invalid_confidence():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    # Mixed observations share one cohort mean, rather than per-model CDF ranks.
    assert screen_scores(table((90., 96.)), ref, min_cohort=1)["gate_q"].to_list() == [.6, .6]
    assert screen_scores(table((None, None)), ref, min_cohort=1)["gate_q"].to_list() == [.8,.8]
    with pytest.raises(ValueError, match="plddt"):
        screen_scores(table((-1., 100.)), ref, min_cohort=1)
    with pytest.raises(ValueError, match="plddt"):
        screen_scores(table().drop("plddt"), ref, min_cohort=1)
    with pytest.raises(ValueError, match="finite component"):
        screen_scores(table().with_columns(enrichment=pl.lit(2.)), ref, min_cohort=1)


def test_frozen_resources_no_identifiers_and_mock_count():
    data = files("tcren.data")
    ref = bundled_reference(CohortKind.TCR_RANKING)
    controls = data.joinpath("motif_mock_controls.tsv.gz").read_bytes()
    assert hashlib.sha256(controls).hexdigest() == ref["mock_controls_sha256"]
    rows = pl.read_csv(gzip.decompress(controls), separator="\t")
    assert rows.columns == ["fisher_statistic", "weight"]
    assert rows.height == 2056 and rows["weight"].sum() == 2056
    assert hashlib.sha256(data.joinpath("motif_residue_reference.json").read_bytes()).hexdigest() == ref["motif_residue_sha256"]
    model = json.loads(data.joinpath("motif_confidence_reference.json").read_bytes())
    assert model == ref["motif_confidence"]
    rows = table().rename({"g_recognition_cyclo_frac":"betti"}).with_columns(
        dataset=pl.lit("sample"), stratum=pl.lit("sample"), epitope=pl.lit("sample"))
    np.testing.assert_allclose(apply_confidence_score(rows,model)["score"],
                               screen_scores(table(),ref,min_cohort=1)["screen_score"])


def test_zero_contact_mock_enrichment_and_six_loop_requirement():
    rows=[]
    for chain in ["TRA", "TRB"]:
        for i, region in enumerate(["CDR1","CDR2","CDR3"]):
            rows.append(dict(id="model", chain_id=chain, residue_index=i, chain=chain,
                             region=region, aa="A", available=5, complete=True,
                             cycles=0, pairs=0, occupied=0))
    rows.append(dict(rows[0], chain_id="P", chain="PEPTIDE", region="PEPTIDE"))
    data=pl.DataFrame(rows)
    out=bundled_enrichment(data)
    assert out["enrichment"].to_list()==[0.]
    assert out["null_models"].to_list()==[2056]
    with pytest.raises(ValueError,match="six CDR"):
        bundled_enrichment(data.filter(pl.col("region")!="CDR1"))


def test_shared_api_typed_nulls_withholds_all_missing_components():
    model = bundled_reference(CohortKind.TCR_RANKING)["motif_confidence"]
    rows = pl.DataFrame({"iptm": [None], "betti": [None], "enrichment": [None],
                         "plddt": [None], "dataset": ["x"], "stratum": ["x"], "epitope": ["x"]})
    out = apply_confidence_score(rows, model)
    assert out["score"].to_list() == [None]
    assert out.select("missing_iptm", "missing_betti", "missing_enrichment").row(0) == (True, True, True)


def test_raw_structure_residue_pipeline_exact_physical_counts():
    from dataclasses import replace
    from tcren.cli import _screen_residue_chunk
    from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure

    def residue(i, z):
        atoms = tuple(Atom(n, n[0], np.array([j*.1, i*.2, z]))
                      for j, n in enumerate(["N", "CA", "C", "O", "CB"]))
        return Residue(i, i+1, "", "A", "ALA", atoms)

    chains = []
    for chain in ["TRA", "TRB"]:
        residues = [residue(i, 0.) for i in range(3)]
        regions = [RegionMarkup(f"CDR{i+1}", i, i, "A", [r]) for i, r in enumerate(residues)]
        chains.append(Chain(chain, residues, chain, "TRAB", regions=regions))
    chains.append(Chain("P", [residue(0, 2.)], "PEPTIDE", "PEPTIDE"))
    mhc = residue(0, 10.)
    chains.append(Chain("M", [mhc], "MHCa", "MHCI",
                        regions=[RegionMarkup("HELIX_A1", 0, 0, "A", [mhc])]))
    structure = Structure("a", chains)
    rows = _screen_residue_chunk([structure, replace(structure, pdb_id="b")])
    receptor = rows.filter(pl.col("chain") != "PEPTIDE")
    peptide = rows.filter(pl.col("chain") == "PEPTIDE")
    assert receptor["pairs"].to_list() == [25]*12
    assert receptor["cycles"].to_list() == [24]*12
    assert peptide["pairs"].to_list() == [150]*2
    assert peptide["cycles"].to_list() == [144]*2
    assert rows["complete"].all()
    result = bundled_enrichment(rows)
    assert result.height == 2 and result["enrichment"][0] == result["enrichment"][1]


def test_nearest_screen_context_preserves_floor_and_missingness():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    data = table((90., 96.))
    per_model = screen_scores(data, ref, min_cohort=1, confidence_context="nearest-plddt")
    assert per_model["gate_q"].to_list() == [.2, 1.]
    grouped = screen_scores(data, ref, min_cohort=1, confidence_context="nearest-mean-plddt")
    assert grouped["gate_q"].to_list() == [.6, .6]
    assert grouped["proxy_context_n"].to_list() == [2, 2]
    assert grouped["cohort_mean_plddt"].to_list() == [93., 93.]
    assert screen_scores(data, ref, confidence_context="nearest-plddt")["screen_score"].null_count() == 2
    missing = screen_scores(table((None, 93.)), ref, min_cohort=1, confidence_context="nearest-plddt")
    assert missing["screen_score"][0] is None
    assert "no observed pLDDT" in missing["score_unavailable"][0]
    assert missing["screen_score"][1] is not None
    with pytest.raises(ValueError, match="requires the TCR"):
        screen_scores(data, bundled_reference(CohortKind.EPITOPE_RANKING), confidence_context="nearest-plddt")
