"""End-to-end pipeline: annotate → superimpose → contacts → per-interface scores."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")

from tcren.pipeline import run, score_row  # noqa: E402

_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def test_pipeline_no_superimpose_scores_three_interfaces():
    res = run(_FIXTURE, superimpose=False)
    assert set(res.scores) == {"tcr_peptide", "tcr_mhc", "peptide_mhc", "total"}
    assert res.oriented is None and res.rmsd is None
    assert res.markup.height > 0 and res.contacts.height > 0
    # total is the sum of the three interface energies
    assert res.scores["total"] == pytest.approx(
        res.scores["tcr_peptide"] + res.scores["tcr_mhc"] + res.scores["peptide_mhc"]
    )


def test_pipeline_superimpose_sets_canonical_frame():
    res = run(_FIXTURE, superimpose=True)
    assert res.oriented is not None and res.rmsd is not None
    assert {c.chain_id for c in res.oriented.chains} <= set("ABCDE")
    row = score_row(res)
    assert row["pdb.id"] == "1ao7" and row["mhc.class"] == "MHCI"
    assert row["F_tcr_pep"] == res.scores["tcr_peptide"]


def test_pipeline_reference_aa_adds_delta_f():
    """--delta path: ΔF per interface, ΔF_tcr_mhc ≡ 0, ΔF_total = ΔF_TP + ΔF_PM."""
    res = run(_FIXTURE, superimpose=False, reference_aa="A")
    assert res.scores["delta_tcr_mhc"] == 0.0          # the peptide is not in that interface
    assert res.scores["delta_total"] == pytest.approx(
        res.scores["delta_tcr_peptide"] + res.scores["delta_peptide_mhc"]
    )
    # ΔF = F(peptide) − F(poly-Ala) is a genuine difference, not a copy of F
    assert res.scores["delta_tcr_peptide"] != res.scores["tcr_peptide"]
    row = score_row(res)
    assert row["dF_total"] == res.scores["delta_total"]
    assert "dF_total" not in score_row(run(_FIXTURE, superimpose=False))
