"""The bundled AAindex3 resource: parsing, the catalogue, identification, and the component split.

The load-bearing claim these guard is that ``mj()`` and ``keskin()`` ARE named AAindex entries.
That identification is what lets the manuscript cite them, so a silent parser regression would be
a citation defect, not just a test failure.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.potential import (aaindex, betancourt, catalogue, entry, identify,
                             keskin, mj, mj1996, parse_aaindex3)


def test_catalogue_shape_and_kinds():
    c = catalogue()
    assert c.height == 47
    kinds = dict(c.group_by("kind").len().iter_rows())
    assert kinds == {"energy": 42, "count": 2, "distance": 3}
    # the three environment-dependent entries whose row and column environments differ
    asym = set(c.filter(~pl.col("symmetric"))["accession"].to_list())
    assert asym == {"ZHAC000102", "ZHAC000103", "ZHAC000105"}


def test_bundled_mj_and_keskin_are_named_aaindex_entries():
    """400 of 400 cells, and unique: the runner-up must be far away, or it is not an identification."""
    for pot, acc, runner_up_at_least in ((mj(), "MIYS990106", 0.5),
                                         (keskin(), "KESO980101", 1.0),
                                         (betancourt(), "BETM990101", 0.5)):
        hits = identify(pot)
        assert hits[0][0] == acc and hits[0][1] == pytest.approx(0.0, abs=1e-9), hits[:3]
        assert hits[1][1] > runner_up_at_least, hits[:3]


def test_our_mj1996_transcription_differs_from_aaindex_in_five_cells():
    """Ours was hand-transcribed from the published table; AAindex is curated. They disagree.

    Five of 210 unique pairs differ by 0.04 to 0.28, correlating at 0.99978, and four of the five
    involve Met, Arg or His. Pinned rather than fixed: the file is left byte-for-byte untouched
    under the same rule as the MJ/Keskin file, and this test makes the discrepancy visible instead
    of letting a later reader rediscover it. ``aaindex("MIYS960101")`` is the curated alternative.
    """
    hits = identify(mj1996())
    assert hits[0][0] == "MIYS960101"
    assert hits[0][1] == pytest.approx(0.28, abs=0.005)
    assert hits[1][1] > 1.0                                # still unique: runner-up off by 1.58
    ours = np.asarray(mj1996().as_matrix()[0], float)
    idx = mj1996().as_matrix()[1]
    e = entry("MIYS960101")
    theirs = np.empty_like(ours)
    for i, a in enumerate(e.rows):
        for j, b in enumerate(e.cols):
            theirs[idx[a], idx[b]] = e.matrix[i, j]
    off = {tuple(sorted((a, b))) for a in idx for b in idx
           if abs(ours[idx[a], idx[b]] - theirs[idx[a], idx[b]]) > 1e-9}
    assert off == {("M", "V"), ("D", "M"), ("E", "M"), ("H", "R"), ("A", "P")}


def test_aaindex_refuses_non_energy_tables():
    aaindex("MOOG990101")                                  # an energy: fine
    for acc in ("MIYS960103", "BONM030104"):               # counts, distances
        with pytest.raises(ValueError, match="not contact energies"):
            aaindex(acc)
        assert entry(acc).kind in ("count", "distance")    # still reachable deliberately
    with pytest.raises(KeyError):
        aaindex("NOSUCH000101")


def test_components_reassemble_and_size_is_the_contact_count():
    for pot in (mj(), mj1996(), keskin(), betancourt()):
        parts = pot.components()
        M = np.asarray(pot.as_matrix()[0], float)
        S = sum(np.asarray(parts[t].as_matrix()[0], float) for t in ("size", "comp", "pair"))
        assert np.allclose(M, S), pot.name
        # "size" is constant, so its interface sum is proportional to the number of contacts
        size = np.asarray(parts["size"].as_matrix()[0], float)
        assert np.allclose(size, size.flat[0])
        # "pair" is double-centred: every row and column sums to zero
        pair = np.asarray(parts["pair"].as_matrix()[0], float)
        assert np.allclose(pair.sum(axis=0), 0.0) and np.allclose(pair.sum(axis=1), 0.0)


def test_parser_rejects_a_malformed_matrix_block():
    bad = ("H FAKE000101\nD test\nA nobody\nT none\nJ none\n"
           "M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV\n"
           "  0.1 0.2\n  0.3\n//\n")
    with pytest.raises(ValueError, match="lower triangle nor a full rectangle"):
        parse_aaindex3(bad)


def test_missing_cells_survive_as_nan_and_are_dropped_from_the_potential():
    e = entry("PARB960101")                                # one residue's row and column omitted
    assert e.n_missing == 39
    p = aaindex("PARB960101")
    assert p.matrix.height == 400 - 39
