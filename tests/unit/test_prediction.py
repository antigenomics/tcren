"""Native generator confidence joins preserve identities, directions and missing evidence."""

import json

import numpy as np
import polars as pl
import pytest

from tcren.prediction import confidence_on_contacts, read_prediction_confidence, inspect_predictions
from tcren.structure.model import Atom, Chain, Residue, Structure
from tcren.validation.md import read_xvg


def test_boltz_native_confidence_requires_producer_mapping_and_retains_direction(tmp_path):
    summary = tmp_path / "confidence_target_model_0.json"
    summary.write_text(
        json.dumps(
            {
                "iptm": 0.8,
                "confidence_score": 0.7,
                "pair_chains_iptm": {"0": {"0": 0.9, "1": 0.4}, "1": {"0": 0.6, "1": 0.8}},
            }
        )
    )
    chains = np.array(
        [("D", 0, 0, 0, 1), ("P", 0, 1, 1, 1)],
        dtype=[
            ("name", "U5"),
            ("mol_type", "i4"),
            ("asym_id", "i4"),
            ("res_idx", "i4"),
            ("res_num", "i4"),
        ],
    )
    residues = np.array(
        [("GLY", 0, True), ("GLY", 0, True)],
        dtype=[("name", "U5"), ("res_idx", "i4"), ("is_standard", "?")],
    )
    producer = tmp_path / "producer.npz"
    np.savez(producer, chains=chains, residues=residues, mask=[True, True])
    pae = tmp_path / "pae_target_model_0.npz"
    np.savez(pae, pae=[[0, 5], [9, 0]])
    np.savez(tmp_path / "plddt_target_model_0.npz", plddt=[0.75, 0.85])
    with pytest.raises(ValueError, match="producer_structure"):
        read_prediction_confidence(summary, generator="boltz2")
    confidence = read_prediction_confidence(
        summary, generator="boltz2", producer_structure_path=producer, detail_path=pae
    )
    joined = confidence_on_contacts(_contacts(), _structure(), confidence)
    assert joined["chain_pair_iptm"].to_list() == [0.4, 0.6]
    assert joined["token_pair_pae"].to_list() == [5, 9]
    assert joined["residue_plddt.from"].to_list() == [75, 85]
    assert confidence.global_values["confidence_score"] == 0.7
    assert "ranking_score" not in confidence.global_values
    assert joined["token_pair_contact_probability"].null_count() == 2
    np.savez(producer, chains=chains, residues=residues, mask=[True, False])
    with pytest.raises(ValueError, match="masked"):
        read_prediction_confidence(summary, generator="boltz2", producer_structure_path=producer)


def test_tfold_retained_confidence_keeps_native_units_and_missing_pae(tmp_path):
    p = tmp_path / "confidence.json"
    raw = {
        "generator": "tfold-tcr",
        "residue_indexing": "sequence_index",
        "plddt_units": "probability_0_1",
        "scalars": {"iptm": 0.8},
        "chain_sequences": {"D": "G", "P": "G"},
        "token_chain_ids": ["D", "P"],
        "token_res_ids": [0, 0],
        "residue_plddts": [0.75, 0.85],
    }
    p.write_text(json.dumps(raw))
    c = read_prediction_confidence(p, generator="tfold-tcr")
    rows = confidence_on_contacts(_contacts(), _structure(), c)
    assert rows["residue_plddt.from"].to_list() == [75, 85]
    assert rows["chain_pair_iptm"].null_count() == 2
    assert rows["token_pair_pae"].null_count() == 2
    raw["plddt_units"] = "unspecified"
    p.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="schema"):
        read_prediction_confidence(p, generator="tfold-tcr")


def test_generic_af2_bundle_retains_pae_and_plddt_without_pair_imputation(tmp_path):
    p = tmp_path / "af2.json"
    p.write_text(
        json.dumps(
            {
                "generator": "af2-multimer",
                "residue_indexing": "sequence_index",
                "plddt_units": "score_0_100",
                "scalars": {"iptm": 0.8},
                "chain_sequences": {"D": "G", "P": "G"},
                "token_chain_ids": ["D", "P"],
                "token_res_ids": [0, 0],
                "residue_plddts": [75, 85],
                "pae": [[0, 5], [9, 0]],
            }
        )
    )
    c = read_prediction_confidence(p, generator="af2-multimer")
    rows = confidence_on_contacts(_contacts(), _structure(), c)
    assert rows["token_pair_pae"].to_list() == [5, 9]
    assert rows["residue_plddt.from"].to_list() == [75, 85]
    assert rows["chain_pair_iptm"].null_count() == 2


def _files(tmp_path):
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "iptm": 0.8,
                "ptm": 0.7,
                "chain_ids": ["D", "P"],
                "chain_pair_iptm": [[0.9, 0.4], [0.6, 0.8]],
                "chain_pair_pae_min": [[1.0, 2.0], [3.0, 4.0]],
            }
        )
    )
    detail = tmp_path / "detail.json"
    detail.write_text(
        json.dumps(
            {
                "token_chain_ids": ["D", "P"],
                "token_res_ids": [203, 8],
                "pae": [[0.0, 5.0], [9.0, 0.0]],
                "contact_probs": [[1.0, 0.3], [0.3, 1.0]],
            }
        )
    )
    return summary, detail


def test_af2_renamed_coordinate_chains_use_exact_sequence_bijection(tmp_path):
    from dataclasses import replace

    structure = _structure()
    structure.chains[1].residues = [replace(structure.chains[1].residues[0], aa="A", resname="ALA")]
    path = tmp_path / "af2.json"
    raw = {
        "generator": "af2-multimer",
        "residue_indexing": "sequence_index",
        "plddt_units": "score_0_100",
        "scalars": {"iptm": 0.8},
        "chain_sequences": {"input_peptide": "A", "input_tcr": "G"},
        "token_chain_ids": ["input_peptide", "input_tcr"],
        "token_res_ids": [0, 0],
        "residue_plddts": [85, 75],
        "pae": [[0, 9], [5, 0]],
    }
    path.write_text(json.dumps(raw))
    confidence = read_prediction_confidence(path, generator="af2-multimer")
    joined = confidence_on_contacts(_contacts(), structure, confidence)
    assert joined["residue_plddt.from"].to_list() == [75, 85]
    assert joined["token_pair_pae"].to_list() == [5, 9]
    assert joined["chain_pair_iptm"].null_count() == 2
    assert confidence.token_sequences == raw["chain_sequences"]
    with pytest.raises(ValueError, match="not bijective"):
        confidence_on_contacts(_contacts().head(0), _structure(), confidence)


def test_af2_renamed_identical_chains_are_ambiguous(tmp_path):
    path = tmp_path / "af2.json"
    path.write_text(
        json.dumps(
            {
                "generator": "af2-multimer",
                "residue_indexing": "sequence_index",
                "plddt_units": "score_0_100",
                "scalars": {},
                "chain_sequences": {"A": "G", "B": "G"},
                "token_chain_ids": ["A", "B"],
                "token_res_ids": [0, 0],
                "residue_plddts": [75, 85],
            }
        )
    )
    confidence = read_prediction_confidence(path, generator="af2-multimer")
    with pytest.raises(ValueError, match="Absent/ambiguous"):
        confidence_on_contacts(_contacts(), _structure(), confidence)


def _structure():
    def residue(author):
        return Residue(0, author, "", "G", "GLY", (Atom("CA", "C", np.zeros(3)),))

    return Structure("sample", [Chain("D", [residue(203)]), Chain("P", [residue(8)])])


def _contacts():
    return pl.DataFrame(
        {
            "chain.id.from": ["D", "P"],
            "residue.index.from": [0, 0],
            "chain.id.to": ["P", "D"],
            "residue.index.to": [0, 0],
        }
    )


def test_directional_pair_and_token_confidence_uses_author_numbers(tmp_path):
    summary, detail = _files(tmp_path)
    confidence = read_prediction_confidence(summary, generator="af3", detail_path=detail)
    joined = confidence_on_contacts(_contacts(), _structure(), confidence)
    assert joined["chain_pair_iptm"].to_list() == [0.4, 0.6]
    assert joined["token_pair_pae"].to_list() == [5.0, 9.0]
    assert joined["token_pair_contact_probability"].to_list() == [0.3, 0.3]
    assert joined["iptm"].to_list() == [0.8, 0.8]
    assert confidence_on_contacts(_contacts().head(0), _structure(), confidence).height == 0


def test_tcrmodel2_requires_literal_rank_and_retains_zero(tmp_path):
    path = tmp_path / "statistics.json"
    path.write_text(
        json.dumps({"ranked_0": {"iptm": "0.0", "plddt": "93.0"}, "ranked_1": {"iptm": 0.2}})
    )
    with pytest.raises(ValueError, match="ranked key"):
        read_prediction_confidence(path, generator="tcrmodel2")
    c = read_prediction_confidence(path, generator="tcrmodel2", model_key="ranked_0")
    joined = confidence_on_contacts(_contacts(), _structure(), c)
    assert joined["iptm"].to_list() == [0.0, 0.0]
    assert joined["token_pair_pae"].null_count() == 2
    assert joined["chain_pair_iptm"].null_count() == 2


def test_malformed_dimensions_and_unknown_token_ids_fail(tmp_path):
    summary, detail = _files(tmp_path)
    c = read_prediction_confidence(summary, generator="af3", detail_path=detail)
    s = _structure()
    s.chains[0].chain_id = "renamed"
    with pytest.raises(KeyError):
        confidence_on_contacts(_contacts(), s, c)
    raw = json.loads(summary.read_text())
    raw["chain_pair_iptm"] = [[0.8]]
    summary.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="matrix"):
        read_prediction_confidence(summary, generator="af3")


def test_manifest_does_not_collapse_repeated_ranked_filenames():
    manifest = pl.DataFrame(
        {
            "complex.id": ["ranked_0", "ranked_0"],
            "structure": ["one/ranked_0.pdb", "two/ranked_0.pdb"],
        }
    )
    with pytest.raises(ValueError, match="unique"):
        inspect_predictions(manifest)


def test_native_confidence_is_absent_not_zero():
    joined = confidence_on_contacts(_contacts(), _structure(), None)
    assert joined["iptm"].null_count() == 2
    assert joined["confidence.generator"].null_count() == 2


def test_reuse_annotation_dispatch_never_calls_search(monkeypatch):
    import tcren.annotation.batch as batch
    import tcren.mhc as mhc
    import tcren.descriptors.table as table

    def fail(*args, **kwargs):
        raise AssertionError("annotation was repeated")

    monkeypatch.setattr(batch, "_batch_annotate", fail)
    monkeypatch.setattr(mhc, "annotate_mhc_batch", fail)
    monkeypatch.setattr(table, "_featurise_one", lambda args: {"complex.id": args[0]})
    assert table.recognition_table([("s", _structure())], annotate=False) == [{"complex.id": "s"}]


def test_xvg_retains_ps_and_rejects_nonfinite_or_repeated_time(tmp_path):
    path = tmp_path / "data.xvg"
    path.write_text('@ xaxis label "Time (ps)"\n# comment\n0 2\n0.1 3\n')
    assert read_xvg(path, columns=2)[:, 0].tolist() == [0.0, 0.1]
    path.write_text("0 nan\n1 2\n")
    with pytest.raises(ValueError, match="nonfinite"):
        read_xvg(path)
    path.write_text("0 2\n0 3\n")
    with pytest.raises(ValueError, match="increase"):
        read_xvg(path)


def test_no_complete_rows_never_evaluate_gaussian(monkeypatch):
    from tcren.score import holdout_model, score_table

    model = holdout_model()
    table = pl.DataFrame(
        {"complex.id": ["unsupported"], **{col: [float("nan")] for col in model.transform.names}}
    )

    def fail(*args, **kwargs):
        raise AssertionError("Gaussian evaluation on zero complete-case rows")

    monkeypatch.setattr(model.joint, "log_odds", fail)
    monkeypatch.setattr(model.joint, "anomaly_on", fail)
    with np.errstate(all="raise"):
        out = score_table(table, iptm=np.asarray([0.72]), model=model)
    assert np.isnan(out["binder_score"][0])
    assert np.isnan(out["confidence_residual"][0])


def test_retained_tcrmodel2_pae_remaps_renumbered_chains_by_sequence(tmp_path):
    summary = tmp_path / "statistics.json"
    summary.write_text(json.dumps({"ranked_0": {"ranking_confidence": 0.8, "iptm": 0.7}}))
    detail = tmp_path / "detail.json"
    detail.write_text(
        json.dumps(
            {
                "generator": "tcrmodel2",
                "residue_indexing": "sequence_index",
                "chain_sequences": {"rawTCR": "A", "rawPep": "G"},
                "token_chain_ids": ["rawTCR", "rawPep"],
                "token_res_ids": [0, 0],
                "pae": [[0.0, 5.0], [9.0, 0.0]],
                "residue_plddts": [80.0, 90.0],
                "scalars": {"ranking_confidence": 0.8},
            }
        )
    )
    s = _structure()
    s.chains[0].residues = [Residue(0, 203, "", "A", "ALA", (Atom("CA", "C", np.zeros(3)),))]
    c = read_prediction_confidence(
        summary, generator="tcrmodel2", model_key="ranked_0", detail_path=detail
    )
    joined = confidence_on_contacts(_contacts(), s, c)
    assert joined["token_pair_pae"].to_list() == [5.0, 9.0]
    assert joined["residue_plddt.from"].to_list() == [80.0, 90.0]
    assert joined["residue_plddt.to"].to_list() == [90.0, 80.0]
    assert joined["chain_pair_iptm"].null_count() == 2
    assert joined["token_pair_contact_probability"].null_count() == 2
    raw = json.loads(detail.read_text())
    raw["scalars"]["ranking_confidence"] = 0.9
    detail.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="requested rank"):
        read_prediction_confidence(
            summary, generator="tcrmodel2", model_key="ranked_0", detail_path=detail
        )


def test_generated_tcrmodel2_retention_and_unique_rank_selection(tmp_path):
    import pickle
    from tcren.validation.tcrmodel2 import export_generated_confidence, ranked_detail_file

    coordinates = tmp_path / "unrelaxed_model.pdb"
    coordinates.write_text(
        "ATOM      1  CA  ALA B   1       0.000   0.000   0.000  1.00 92.00           C\nEND\n"
    )
    result = tmp_path / "result_model.pkl"
    result.write_bytes(
        pickle.dumps(
            {
                "plddt": np.asarray([92.0]),
                "predicted_aligned_error": np.zeros((1, 1)),
                "ranking_confidence": 0.8,
                "iptm": 0.7,
                "ptm": 0.9,
            }
        )
    )
    out = tmp_path / "retained/model_1.json"
    payload = export_generated_confidence(result, coordinates, out)
    assert payload["chain_sequences"] == {"B": "A"}
    assert payload["token_res_ids"] == [0]
    assert (
        ranked_detail_file(out.parent, {"ranked_0": {"ranking_confidence": 0.8}}, "ranked_0") == out
    )
    (out.parent / "duplicate.json").write_text(out.read_text())
    with pytest.raises(ValueError, match="unique producer identity"):
        ranked_detail_file(out.parent, {"ranked_0": {"ranking_confidence": 0.8}}, "ranked_0")


@pytest.mark.parametrize("entrypoint", ["stdin", "ipykernel"])
def test_interactive_parallel_reuses_annotated_structures(monkeypatch, entrypoint):
    import sys
    from tcren.descriptors.table import recognition_table
    from tcren.runtime import available_cpus

    if available_cpus() < 2:
        pytest.skip("two allocated CPUs required")
    if entrypoint == "stdin":
        monkeypatch.setattr(sys.modules["__main__"], "__file__", "<stdin>", raising=False)
    else:
        monkeypatch.setitem(sys.modules, "ipykernel", object())
    rows = recognition_table([("one", _structure()), ("two", _structure())],
                             threads=2, annotate=False, include=[])
    assert rows == [{"complex.id": "one"}, {"complex.id": "two"}]


@pytest.mark.parametrize("historical", [False, True])
def test_inspection_default_never_loads_historical_assessment(tmp_path, monkeypatch, historical):
    from dataclasses import dataclass
    import tcren.annotation
    import tcren.annotation.batch
    import tcren.descriptors.table
    import tcren.docking
    import tcren.mhc
    import tcren.score
    import tcren.structure

    path = tmp_path / "model.pdb"
    path.write_text("END\n")
    structure = _structure()
    structure.pdb_id = "synthetic"
    residue = structure.chains[0].residues[0]
    structure.chains = [Chain(cid, [residue], chain_type=role) for cid, role in
                        [("A", "TRA"), ("B", "TRB"), ("C", "PEPTIDE"), ("D", "MHCa")]]
    monkeypatch.setattr(tcren.structure, "import_structure", lambda *args, **kwargs: structure)
    monkeypatch.setattr(tcren.annotation.batch, "annotate_batch", lambda *args, **kwargs: [None])
    monkeypatch.setattr(tcren.annotation, "classify_chains", lambda *args, **kwargs: None)
    monkeypatch.setattr(tcren.mhc, "annotate_mhc_batch", lambda *args, **kwargs: None)
    monkeypatch.setattr(tcren.descriptors.table, "recognition_table",
                        lambda *args, **kwargs: [{"complex.id": "synthetic"}])
    @dataclass
    class Geometry:
        value: float = 0
    monkeypatch.setattr(tcren.docking, "docking_angles", lambda *args: Geometry())
    monkeypatch.setattr(tcren.docking, "tcr_placement", lambda *args: Geometry())
    def fail(*args, **kwargs):
        raise AssertionError("default inspection must not load fitted readouts")
    monkeypatch.setattr(tcren.score, "score_table", fail)
    monkeypatch.setattr(tcren.score, "holdout_model", fail)
    if historical:
        def historical_scores(table, *, iptm=None):
            if iptm is not None:
                return pl.DataFrame({"complex.id": ["synthetic"], "confidence_residual": [float("nan")],
                                     "binder_iptm": [float("nan")], "confidence_residual_available": [False],
                                     "binder_iptm_available": [False]})
            return pl.DataFrame({"complex.id": ["synthetic"], "assessment_available": [True],
                                 "missing_descriptors": ["[]"]})
        monkeypatch.setattr(tcren.score, "score_table", historical_scores)
    tables = inspect_predictions(pl.DataFrame({"complex.id": ["synthetic"], "structure": [str(path)],
                                                "iptm": [0.7]}), historical_assessment=historical)
    if historical:
        assert tables["assessment"]["confidence_residual_available"].to_list() == [False]
        assert tables["assessment"]["binder_iptm_available"].to_list() == [False]
    else:
        assert tables["assessment"].is_empty()
        assert tables["models"]["assessment_status"].to_list() == ["historical_not_requested"]


@pytest.mark.parametrize("generator", ["protenix", "rf3"])
def test_native_summary_only_providers_keep_scalar_semantics(tmp_path, generator):
    path = tmp_path / "summary.json"
    path.write_text(json.dumps({"iptm": 0.0, "ptm": 0.8, "ranking_score": -99.84,
                                "chain_pair_iptm": [[0.1, 0.2], [0.3, 0.4]]}))
    confidence = read_prediction_confidence(path, generator=generator)
    assert confidence.generator == generator
    assert confidence.global_values["iptm"] == 0.0
    assert confidence.global_values["ranking_score"] == -99.84
    assert confidence.chain_pairs == {}
    assert confidence.token_index == {}
    assert confidence.token_mapping_status == "unavailable_unverified_token_mapping"
    with pytest.raises(ValueError, match="token mapping is unverified"):
        read_prediction_confidence(path, generator=generator, detail_path=path)
