"""Configurable per-interface potentials (S1): defaults unchanged, overrides take effect."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")

from tcren.pipeline import _resolve_potentials, run  # noqa: E402
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
