"""The compiled fitting-alignment extension (``tcren._align``) + Bio fallback agreement."""

from __future__ import annotations

import pytest

_align = pytest.importorskip("tcren._align")


def test_fitting_score_identical_is_self_match():
    # A sequence fully placed against itself scores its BLOSUM62 self-similarity (no gaps).
    s = "ACDEFGHIKLMNPQRSTVWY"
    assert _align.fitting_score(s, s) > 100  # all positives, no gaps


def test_fitting_score_free_target_skips():
    # The short "placed" string threads through the longer "free" string at no skip cost,
    # so embedding it in flanking residues does not lower the score.
    placed = "WWWW"
    bare = _align.fitting_score(placed, "WWWW")
    padded = _align.fitting_score(placed, "AAAWWWWAAA")
    assert padded == bare


def test_best_hit_picks_exact_candidate():
    free = "AAAAKLMNPQRSTAAAA"
    candidates = ["KLMNPQRST", "WWWWWWWWW", "DDDDDDDDD"]
    idx, score = _align.best_hit(free, candidates)
    assert idx == 0
    assert score == _align.fitting_score(candidates[0], free)


def test_align_returns_matched_positions():
    placed = "KLMNP"
    free = "AAAKLMNPAAA"
    pairs = _align.align(placed, free)
    # every placed residue lands on its identical free residue (offset by the 3-residue pad)
    assert [p for p, _ in pairs] == [0, 1, 2, 3, 4]
    assert [c for _, c in pairs] == [3, 4, 5, 6, 7]


def test_matches_biopython_fitting_aligner():
    # The C++ score must equal Biopython's fitting configuration (the Python fallback).
    pytest.importorskip("Bio")
    from tcren.mhc.pseudo import _pseudo_aligner

    placed = "YFAMYGEKVAHTHVDTLYVRYHYYTWAVWAYTWY"
    free = "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGE"
    assert _align.fitting_score(placed, free) == _pseudo_aligner().score(placed, free)
