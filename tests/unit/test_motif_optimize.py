"""Confidence gates preserve a convex three-component formula."""

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score
from tcren.contacts.motif_optimize import confidence_gate, gated_score, _gate_chunk


def test_gate_endpoints_and_cohort_means():
    x = np.array([[0.8, 0.2, 0.4], [0.3, 0.7, 0.5]])
    np.testing.assert_allclose(gated_score(x, [0, 0], [0.2, 0.5, 0.3]), x[:, 0])
    np.testing.assert_allclose(
        gated_score(x, [1, 1], [0.2, 0.5, 0.3]), x @ np.array([0.2, 0.5, 0.3])
    )
    rows = pl.DataFrame(
        dict(
            dataset=["d", "d"],
            stratum=["s", "s"],
            epitope=["p", "p"],
            iptm=[0.2, 0.8],
            plddt=[80.0, 100.0],
        )
    )
    np.testing.assert_allclose(confidence_gate(rows, "cohort-iptm"), [0.5, 0.5])
    np.testing.assert_allclose(confidence_gate(rows, "complex-plddt"), [0.8, 1])
    np.testing.assert_allclose(confidence_gate(rows, "cohort-plddt"), [0.9, 0.9])


def test_gate_candidate_ap_matches_scalar_formula():
    x = np.array([[0.8, 0.2, 0.4], [0.3, 0.7, 0.5], [0.6, 0.4, 0.8], [0.5, 0.3, 0.7]])
    y = np.array([1, 0, 1, 0])
    q = np.array([0.8, 0.3, 0.6, 0.5])
    w = np.array([[1, 0, 0], [0.2, 0.5, 0.3]])
    ap = _gate_chunk((x, q, y, np.array(["p"] * 4), w, False, False))
    for j, weights in enumerate(w):
        np.testing.assert_allclose(ap[j, 0], average_precision_score(y, gated_score(x, q, weights)))
    score = gated_score(x, q, w[1])
    assert ((score >= 0) & (score <= 1)).all()


def test_joint_residue_tail_retains_type_occurrence_information():
    from tcren.contacts.residue_surprise import fit_reference, score_residues

    rows = pl.DataFrame(
        dict(
            id=[str(i) for i in range(10)],
            chain=["TRB"] * 10,
            region=["CDR3"] * 10,
            aa=["A"] * 9 + ["R"],
            available=[5] * 10,
            complete=[True] * 10,
            occupied=[1] * 10,
            pairs=[3] * 10,
            cycles=[2] * 10,
        )
    )
    ref = fit_reference(rows)
    query = rows.tail(1)
    assert score_residues(query, ref, "cycles")["tail_probability"].item() == 1
    assert (
        score_residues(query, ref, "cycles", joint_residue=True)["tail_probability"].item()
        == 2 / 11
    )
    absent = query.with_columns(pl.lit("W").alias("aa"))
    assert (
        score_residues(absent, ref, "cycles", joint_residue=True)["tail_probability"].item()
        == 1 / 11
    )
    zero = query.with_columns(pl.lit(0).alias("cycles"))
    assert score_residues(zero, ref, "cycles", joint_residue=True)["tail_probability"].item() == 1


def test_residue_audit_marks_rs_flanks_and_peptide_gigi(tmp_path):
    from tcren.contacts.motif_optimize import residue_audit_command

    rows = []
    for chain, region, sequence in [("TRB", "CDR3", "ARSG"), ("PEPTIDE", "PEPTIDE", "GIGI")]:
        for i, aa in enumerate(sequence):
            rows.append(
                dict(
                    id="a",
                    chain=chain,
                    region=region,
                    residue_index=i + 1,
                    aa=aa,
                    pairs=1,
                    cycles=0,
                    bits=1.0,
                )
            )
    pl.DataFrame(rows).write_parquet(tmp_path / "res.parquet")
    pl.DataFrame(
        dict(id=["a"], dataset=["d"], stratum=["s"], epitope=["GILGFVFTL"], y=[1])
    ).write_parquet(tmp_path / "labels.parquet")
    residue_audit_command(
        residues=tmp_path / "res.parquet",
        labels=tmp_path / "labels.parquet",
        output=tmp_path / "out",
    )
    detail = pl.read_parquet(tmp_path / "out/residues.parquet")
    assert detail["GIL_RS_flanks"].sum() == 4
    assert detail["peptide_GIGI"].sum() == 4
    assert detail["fisher_share"].sum() == 1


def test_covalent_pair_patch_counts_shared_partners_once(tmp_path):
    from tcren.contacts.motif_optimize import pair_counts_command

    records = []
    for chain_id, chain, region, sequence in [
        ("T", "TRA", "CDR3", "AR"),
        ("P", "PEPTIDE", "PEPTIDE", "GI"),
    ]:
        for i, aa in enumerate(sequence):
            records.append(
                dict(
                    id="a",
                    chain_id=chain_id,
                    chain=chain,
                    region=region,
                    residue_index=i + 1,
                    aa=aa,
                    available=5,
                    complete=True,
                    pairs=2,
                    occupied=1,
                    partners=2,
                    cycles=0,
                    atom_scope="all",
                )
            )
    pl.DataFrame(records).write_parquet(tmp_path / "res.parquet")
    contacts = []
    for t in [1, 2]:
        for p in [1, 2]:
            contacts.append(
                {
                    "pdb.id": "a",
                    "chain.id.from": "T",
                    "residue.index.from": t,
                    "atom.from": "CA",
                    "chain.id.to": "P",
                    "residue.index.to": p,
                    "atom.to": "CA",
                }
            )
    pl.DataFrame(contacts).write_parquet(tmp_path / "contacts.parquet")
    pair_counts_command(
        residues=tmp_path / "res.parquet",
        contacts=[tmp_path / "contacts.parquet"],
        output=tmp_path / "pairs.parquet",
    )
    out = pl.read_parquet(tmp_path / "pairs.parquet")
    assert set(out["aa"]) == {"AR", "GI"}
    assert out["pairs"].to_list() == [4, 4] and out["partners"].to_list() == [2, 2]
    assert out["cycles"].to_list() == [2, 2] and out["occupied"].to_list() == [2, 2]


def test_odds_gate_is_logit_offset_and_score_is_bounded():
    from scipy.special import expit, logit
    from tcren.contacts.motif_optimize import odds_gate, odds_score

    q = np.array([0.2, 0.5, 0.9])
    a = 0.13
    np.testing.assert_allclose(odds_gate(q, a), expit(logit(q) + logit(a)))
    x = np.array([[0.2, 0.4, 0.8], [0.5, 0.7, 0.4], [0.9, 0.2, 0.5]])
    for p in [[0, 0], [1, 1], [0.13, 0.3]]:
        result = odds_score(x, q, p)
        assert ((result >= 0) & (result <= 1)).all()
    np.testing.assert_allclose(odds_score(x, q, [0, 0]), x[:, 0])
    np.testing.assert_allclose(odds_score(x, q, [1, 1]), x[:, 2])


def test_complementary_weights_reverse_betti_and_motif_with_confidence():
    from tcren.contacts.motif_optimize import complementary_score

    x = np.array([[0.4, 0.8, 0.2], [0.4, 0.8, 0.2], [0.4, 0.8, 0.2]])
    np.testing.assert_allclose(complementary_score(x, [0, 0.5, 1], [1, 1]), [0.2, 0.5, 0.8])
    for b, m in [(0, 0), (0.5, 0.3), (1, 1)]:
        q = np.array([0, 0.5, 1])
        wB = b * q
        wM = m * (1 - q)
        np.testing.assert_allclose(
            complementary_score(x, q, [b, m]), (1 - wB - wM) * x[:, 0] + wB * x[:, 1] + wM * x[:, 2]
        )


def test_frozen_confidence_score_needs_no_labels_and_survives_row_order(tmp_path):
    import json
    from tcren.contacts.motif_optimize import score_command

    rows = pl.DataFrame(
        dict(
            id=["a", "b", "c", "d"],
            dataset=["new"] * 4,
            stratum=["s"] * 4,
            epitope=["e1", "e1", "e2", "e2"],
            iptm=[0.7, 0.8, 0.6, 0.5],
            betti=[0.2] * 4,
            enrichment=[0.4] * 4,
            plddt=[79.0, 81.0, 89.0, 91.0],
        )
    )
    model = dict(
        schema="tcren-confidence-motif-1",
        parameters=dict(b=0.41, m=0.22),
        gate="cohort-plddt",
        rank_confidence=True,
        confidence_cdf=[0.8, 0.9],
        imputation=[0.7, 0.2, 0.4],
        plddt_imputation=85.0,
        formula="S=(1-b*q-m*(1-q))*ipTM+b*q*Betti+m*(1-q)*enrichment",
    )
    (tmp_path / "model.json").write_text(json.dumps(model))
    for tag, table in [("forward", rows), ("reverse", rows.reverse())]:
        table.write_parquet(tmp_path / (tag + ".parquet"))
        score_command(
            features=tmp_path / (tag + ".parquet"),
            resource=tmp_path / "model.json",
            output=tmp_path / (tag + ".tsv"),
        )
    a = pl.read_csv(tmp_path / "forward.tsv", separator="\t").sort("id")
    b = pl.read_csv(tmp_path / "reverse.tsv", separator="\t").sort("id")
    np.testing.assert_array_equal(a["score"], b["score"])
    np.testing.assert_allclose(a["gate_q"], [0.5, 0.5, 1, 1])
    np.testing.assert_allclose(a["w_I"] + a["w_B"] + a["w_M"], 1)


def test_v_segment_audit_uses_recorded_one_based_alignment_end(tmp_path):
    import json
    from tcren.contacts.motif_optimize import residue_audit_command

    rows = pl.DataFrame(
        dict(
            id=["a"] * 4,
            chain_id=["B"] * 4,
            chain=["TRB"] * 4,
            region=["CDR3"] * 4,
            residue_index=[0, 1, 2, 3],
            aa=list("ARSG"),
            pairs=[1] * 4,
            cycles=[0] * 4,
            bits=[1.0] * 4,
        )
    )
    rows.write_parquet(tmp_path / "res.parquet")
    pl.DataFrame(
        dict(id=["a"], dataset=["d"], stratum=["s"], epitope=["GLCTLVAML"], y=[1])
    ).write_parquet(tmp_path / "labels.parquet")
    (tmp_path / "annotation.json").write_text(
        json.dumps(
            {
                "a": {
                    "human": {
                        "B": dict(
                            locus="TRB", v_call="TRBV20-1*01", v_sequence_start=1, v_sequence_end=2
                        )
                    }
                }
            }
        )
    )
    residue_audit_command(
        residues=tmp_path / "res.parquet",
        labels=tmp_path / "labels.parquet",
        annotations=[tmp_path / "annotation.json"],
        output=tmp_path / "out",
    )
    detail = pl.read_parquet(tmp_path / "out/residues.parquet").sort("residue_index")
    assert detail["V_aligned"].to_list() == [True, True, False, False]
    assert detail["GLC_V_aligned_CDR3"].sum() == 2
