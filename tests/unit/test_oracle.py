"""Unit tests for the S5 facade :func:`tcren.summarize_structure`.

The facade only orchestrates S1-S4, so the checks are structural (all five frames
with the expected columns) plus the byte-exact identity that its ``scores`` frame
reproduces :func:`tcren.pipeline.run`'s scores under default arguments.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")  # facade runs the full pipeline (arda / mmseqs annotation)

from tcren import summarize_structure  # noqa: E402
from tcren.pipeline import run  # noqa: E402
from tcren.structure.io import import_structure  # noqa: E402

_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"

_EXPECTED_COLUMNS = {
    "scores": ["pdb.id", "rmsd", "tcr_peptide", "tcr_mhc", "peptide_mhc", "total"],
    "rank": ["pdb.id", "peptide", "score", "rank_pct", "n_background"],
    "ddg": ["pos", "wt_aa", "ddG"],
    "markup": [
        "chain.id", "residue.index", "chain.type", "chain.supertype",
        "region.type", "region.start", "residue.aa",
    ],
}


def test_summarize_returns_five_frames_with_expected_columns():
    out = summarize_structure(_FIXTURE, superimpose=False, background=50)
    assert set(out) == {"scores", "rank", "ddg", "markup", "contacts"}
    for key, cols in _EXPECTED_COLUMNS.items():
        assert out[key].columns == cols, key
    # one summary row each for scores and rank
    assert out["scores"].height == 1
    assert out["rank"].height == 1
    # markup / contacts are non-empty per-residue / per-contact tables
    assert out["markup"].height > 0
    assert out["contacts"].height > 0
    assert "pdb.id" in out["contacts"].columns


def test_ddg_empty_without_alanine_flag():
    out = summarize_structure(_FIXTURE, superimpose=False, background=50)
    # default alanine=False -> empty ddg frame, but with the alanine-scan schema
    assert out["ddg"].height == 0
    assert out["ddg"].columns == ["pos", "wt_aa", "ddG"]


def test_alanine_flag_runs_per_position_scan():
    out = summarize_structure(_FIXTURE, superimpose=False, background=50, alanine=True)
    ddg = out["ddg"]
    # one row per peptide position
    native = out["rank"]["peptide"][0]
    assert ddg.height == len(native)
    assert ddg["wt_aa"].to_list() == list(native)


def test_scores_reproduce_run_byte_exact():
    # Facade scores must equal run()'s scores exactly under matching arguments.
    out = summarize_structure(_FIXTURE, superimpose=False, background=50)
    res = run(import_structure(_FIXTURE), superimpose=False)
    row = out["scores"].to_dicts()[0]
    assert row["pdb.id"] == res.pdb_id
    for key in ("tcr_peptide", "tcr_mhc", "peptide_mhc", "total"):
        assert row[key] == res.scores[key]


def test_rank_uses_native_peptide_and_is_deterministic():
    a = summarize_structure(_FIXTURE, superimpose=False, background=100, seed=0)
    b = summarize_structure(_FIXTURE, superimpose=False, background=100, seed=0)
    ra, rb = a["rank"].to_dicts()[0], b["rank"].to_dicts()[0]
    assert ra == rb  # same seed -> identical rank row
    assert ra["peptide"] == "LLFGYPVYV"  # the structure's native peptide
    assert ra["n_background"] == 100
    assert 0.0 <= ra["rank_pct"] <= 1.0
