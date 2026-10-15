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


def test_traceback_attains_fitting_score_across_gap_transitions():
    """Both gap-state traceback transitions must preserve the DP optimum."""
    from itertools import product
    import random

    def check(placed, free):
        pairs = _align.align(placed, free)
        score = sum(11 if placed[i] == free[j] == 'W' else
                    4 if placed[i] == free[j] == 'A' else -3 for i, j in pairs)
        for (i, j), (next_i, next_j) in zip(pairs, pairs[1:]):
            assert next_i > i and next_j > j
            gap = next_i - i - 1
            if gap:
                score -= 11 + gap - 1
        assert score == _align.fitting_score(placed, free), (placed, free, pairs)

    check('WAWAAWAWAWAA', 'WWAAWWA')
    short = [''.join(s) for n in range(5) for s in product('WA', repeat=n)]
    for placed in short:
        for free in short:
            check(placed, free)
    rng = random.Random(3)
    for _ in range(500):
        check(''.join(rng.choices('WA', k=12)), ''.join(rng.choices('WA', k=7)))
