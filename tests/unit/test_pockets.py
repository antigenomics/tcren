"""A–F pocket markers along the projected peptide track (pure polars)."""

from __future__ import annotations

import polars as pl

from tcren.project2d.pockets import pocket_markers


def _markup(uv):
    return pl.DataFrame({
        "complex_chain": ["peptide"] * len(uv) + ["mhca"],
        "aa_index": list(range(len(uv))) + [0],
        "u": [p[0] for p in uv] + [99.0],
        "v": [p[1] for p in uv] + [99.0],
    })


def test_pocket_markers_span_the_track():
    # A straight peptide track from (0,0) to (10,0); A..F evenly spaced, endpoints anchored.
    pk = pocket_markers(_markup([(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]))
    assert pk["pocket"].to_list() == ["A", "B", "C", "D", "E", "F"]
    assert pk.row(0, named=True)["u"] == 0.0          # A at the N-terminus
    assert pk.row(5, named=True)["u"] == 10.0         # F at the C-terminus
    assert pk["u"].to_list() == sorted(pk["u"].to_list())  # monotone along the track


def test_pocket_markers_empty_without_peptide():
    assert pocket_markers(_markup([(0.0, 0.0)])).height == 0   # <2 projected peptide residues
