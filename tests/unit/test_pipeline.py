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
    # Phi_TCRpMHC = c_TP Phi_TP + c_TM Phi_TM + c_PM Phi_PM. The three interfaces are scored with
    # different potentials, whose matrices are not on a common scale, so the total is the
    # Native2026-normalised sum rather than the raw one -- an unweighted sum would be 2.6x more
    # sensitive to a presentation contact than to a recognition one.
    from tcren.pipeline import _INTERFACE_POTENTIAL, _phi_scale, _resolve_potentials

    pots = _resolve_potentials(None)
    want = sum(res.scores[i] / _phi_scale(i, pots[i]) for i in _INTERFACE_POTENTIAL)
    assert res.scores["total"] == pytest.approx(want)
    assert res.scores["total"] != pytest.approx(
        res.scores["tcr_peptide"] + res.scores["tcr_mhc"] + res.scores["peptide_mhc"])


def test_pipeline_superimpose_sets_canonical_frame():
    res = run(_FIXTURE, superimpose=True)
    assert res.oriented is not None and res.rmsd is not None
    assert {c.chain_id for c in res.oriented.chains} <= set("ABCDE")
    row = score_row(res)
    assert row["pdb.id"] == "1ao7" and row["mhc.class"] == "MHCI"
    assert row["Phi_tcr_pep"] == res.scores["tcr_peptide"]


def test_pipeline_reference_aa_adds_delta_f():
    """--delta path: ΔΦ per interface, ΔΦ_tcr_mhc ≡ 0, ΔΦ_total the weighted TP + PM sum."""
    from tcren.pipeline import _phi_scale, _resolve_potentials

    res = run(_FIXTURE, superimpose=False, reference_aa="A")
    assert res.scores["delta_tcr_mhc"] == 0.0          # the peptide is not in that interface
    pots = _resolve_potentials(None)
    assert res.scores["delta_total"] == pytest.approx(
        res.scores["delta_tcr_peptide"] / _phi_scale("tcr_peptide", pots["tcr_peptide"])
        + res.scores["delta_peptide_mhc"] / _phi_scale("peptide_mhc", pots["peptide_mhc"])
    )
    # ΔΦ = Φ(peptide) − Φ(poly-Ala) is a genuine difference, not a copy of Φ
    assert res.scores["delta_tcr_peptide"] != res.scores["tcr_peptide"]
    row = score_row(res)
    assert row["dPhi_total"] == res.scores["delta_total"]
    assert "dPhi_total" not in score_row(run(_FIXTURE, superimpose=False))
