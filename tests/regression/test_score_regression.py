"""End-to-end scoring parity against the committed example output."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from tcren import ContactMap, parse_structure, score_peptides
from tcren.annotation import classify_chains
from tcren.potential import tcren

pytest.importorskip("arda")

REPO = Path(__file__).resolve().parents[2]
EXAMPLE = REPO / "example"


def _read_candidates(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and line.strip().lower() != "peptide"
    ]


def test_example_scores_match_oracle():
    sdir = EXAMPLE / "input_structures"
    expected = pl.read_csv(EXAMPLE / "output_TCRen" / "candidate_epitopes_TCRen.csv")
    candidates = _read_candidates(EXAMPLE / "candidate_epitopes.txt")

    frames = []
    for fp in sorted(sdir.glob("*.pdb")):
        s = parse_structure(fp, pdb_id=fp.name)
        classify_chains(s, organism="human")
        cm = ContactMap.from_structure(s)
        frames.append(score_peptides(cm, candidates, tcren()))
    got = pl.concat(frames)

    assert got.height == expected.height
    joined = got.join(
        expected, on=["complex.id", "peptide"], how="full", coalesce=True, suffix="_exp"
    )
    assert joined.filter(
        pl.col("score").is_null() | pl.col("score_exp").is_null()
    ).height == 0
    max_abs = joined.select((pl.col("score") - pl.col("score_exp")).abs().max()).item()
    assert max_abs == pytest.approx(0.0, abs=1e-9), f"max abs score diff = {max_abs}"

    # Ranking (sorted ascending) must match exactly.
    got_order = got.sort("score")["peptide"].to_list()
    exp_order = expected.sort("score")["peptide"].to_list()
    assert got_order == exp_order
