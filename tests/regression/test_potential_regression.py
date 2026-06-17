"""Numerical-fidelity tests: derived potentials must match the committed oracles.

These run without any PDB parsing or arda — they validate the potential derivation in
isolation against ``TCRen_potential.csv`` (classic) and ``tcren_am/tcren.txt`` (am).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from tcren.potential import derive_tcren
from tcren.potential.model import Potential

REPO = Path(__file__).resolve().parents[2]
CONTACT_MAPS = REPO / "tests" / "assets" / "oracle" / "data" / "contact_maps_PDB.csv"
SUMMARY = REPO / "tests" / "assets" / "oracle" / "data" / "summary_PDB_structures.csv"
TCREN_CSV = REPO / "tests" / "assets" / "oracle" / "data" / "TCRen_potential.csv"
TCREN_AM_TXT = REPO / "tests" / "assets" / "oracle" / "tcren_am" / "tcren.txt"


def _nonred_ids() -> list[str]:
    summary = pl.read_csv(SUMMARY)
    return summary.filter(pl.col("nonred"))["pdb.id"].to_list()


def _join_on_pairs(got: pl.DataFrame, want: pl.DataFrame, val_got: str, val_want: str):
    keys = ["residue.aa.from", "residue.aa.to"]
    j = got.select(*keys, pl.col(val_got)).join(
        want.select(*keys, pl.col(val_want)), on=keys, how="full", coalesce=True
    )
    assert j.filter(pl.col(val_got).is_null() | pl.col(val_want).is_null()).height == 0, (
        "pair sets differ between derived and oracle potentials"
    )
    return j


def test_classic_tcren_matches_oracle():
    contacts = pl.read_csv(CONTACT_MAPS)
    pot = derive_tcren(contacts, include=_nonred_ids(), variant="classic")
    oracle = pl.read_csv(TCREN_CSV)  # columns: residue.aa.from, residue.aa.to, TCRen

    assert pot.matrix.height == 380  # 19 (no Cys on 'from') x 20

    j = _join_on_pairs(pot.matrix, oracle, "value", "TCRen")
    max_abs = j.select((pl.col("value") - pl.col("TCRen")).abs().max()).item()
    assert max_abs == pytest.approx(0.0, abs=1e-9), f"max abs diff = {max_abs}"


# Frozen snapshot of the contact data that tcren_am/tcren.txt was derived from
# (repo commit 8f5dfe9, 2025-08-15). The committed tcren.txt is stale relative to the
# current contact_maps_PDB.csv (regenerated 2026-05-18), so the am oracle is pinned
# against this historical input rather than the live data. See test below.
AM_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "am_oracle"


def test_am_tcren_matches_oracle():
    contacts = pl.read_csv(AM_FIXTURE / "contact_maps.csv")
    nonred = pl.read_csv(AM_FIXTURE / "summary.csv").filter(pl.col("nonred"))["pdb.id"].to_list()
    pot = derive_tcren(contacts, include=nonred, variant="am")
    oracle = pl.read_csv(TCREN_AM_TXT, separator="\t")  # from, to, TCRen, count

    assert pot.matrix.height == 21 * 21

    j = _join_on_pairs(pot.matrix, oracle, "value", "TCRen")
    max_abs = j.select((pl.col("value") - pl.col("TCRen")).abs().max()).item()
    assert max_abs == pytest.approx(0.0, abs=1e-9), f"max abs diff = {max_abs}"

    # The gap/gap cell is seeded with the contact count and must round-trip.
    cj = pot.matrix.join(
        oracle.select("residue.aa.from", "residue.aa.to", "count"),
        on=["residue.aa.from", "residue.aa.to"],
        how="inner",
    )
    assert (cj["count"] == cj["count_right"]).all()


def test_am_derivation_runs_on_current_data():
    """The am variant must also derive cleanly on the live contact map (no oracle)."""
    contacts = pl.read_csv(CONTACT_MAPS)
    pot = derive_tcren(contacts, include=_nonred_ids(), variant="am")
    assert pot.matrix.height == 21 * 21
    assert pot.matrix["value"].is_finite().all()


def test_potential_roundtrip_csv(tmp_path):
    contacts = pl.read_csv(CONTACT_MAPS)
    pot = derive_tcren(contacts, include=_nonred_ids(), variant="classic")
    out = tmp_path / "p.csv"
    pot.to_csv(out)
    again = Potential.from_csv(out, name="TCRen", value_col="value")
    assert again.matrix.sort("residue.aa.from", "residue.aa.to").equals(
        pot.matrix.sort("residue.aa.from", "residue.aa.to")
    )
