"""Configurable per-interface potentials (S1): defaults unchanged, overrides take effect."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")

from tcren.pipeline import _resolve_potentials, run, score_row  # noqa: E402
from tcren.potential import mj, tcren  # noqa: E402

_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def test_resolve_defaults_match_interface_potential():
    resolved = _resolve_potentials(None)
    assert resolved["tcr_peptide"].name == tcren().name
    assert resolved["tcr_mhc"].name == mj().name
    assert resolved["peptide_mhc"].name == mj().name


def test_default_equals_explicit_equal_mapping():
    # Default (None) must reproduce the explicit per-interface family mapping byte-for-byte.
    res_default = run(_FIXTURE, superimpose=False)
    res_explicit = run(
        _FIXTURE,
        superimpose=False,
        potentials={"tcr_peptide": "tcren", "tcr_mhc": "mj", "peptide_mhc": "mj"},
    )
    assert res_default.scores == res_explicit.scores


def test_swapping_tcr_mhc_to_tcren_changes_score():
    res_default = run(_FIXTURE, superimpose=False)
    res_swapped = run(_FIXTURE, superimpose=False, potentials={"tcr_mhc": "tcren"})
    # The TCR↔MHC interface now uses TCRen instead of MJ, so its energy must change.
    assert res_swapped.scores["tcr_mhc"] != res_default.scores["tcr_mhc"]
    # The other interfaces keep their default potential and are unchanged.
    assert res_swapped.scores["tcr_peptide"] == res_default.scores["tcr_peptide"]
    assert res_swapped.scores["peptide_mhc"] == res_default.scores["peptide_mhc"]


def test_intra_peptide_term_is_off_by_default():
    scores = run(_FIXTURE, superimpose=False).scores
    assert "peptide_internal" not in scores
    assert "F_pep_int" not in score_row(run(_FIXTURE, superimpose=False))


def test_intra_weight_reports_the_term_raw_and_folds_it_into_the_total():
    default = run(_FIXTURE, superimpose=False)
    weighted = run(_FIXTURE, superimpose=False, intra_weight=3.0)
    # The three interface energies are untouched; only the total absorbs the new term.
    for iface in ("tcr_peptide", "tcr_mhc", "peptide_mhc"):
        assert weighted.scores[iface] == default.scores[iface]
    term = weighted.scores["peptide_internal"]
    assert term != 0.0
    assert weighted.scores["total"] == pytest.approx(default.scores["total"] + 3.0 * term)
    assert score_row(weighted)["F_pep_int"] == term


def test_intra_peptide_potential_defaults_to_mj_and_is_overridable():
    assert _resolve_potentials(None)["peptide_internal"].name == mj().name
    default = run(_FIXTURE, superimpose=False, intra_weight=1.0)
    swapped = run(_FIXTURE, superimpose=False, intra_weight=1.0,
                  potentials={"peptide_internal": "keskin"})
    assert swapped.scores["peptide_internal"] != default.scores["peptide_internal"]
    assert swapped.scores["tcr_peptide"] == default.scores["tcr_peptide"]
