"""Unit tests for the leave-one-out n199 path in ``scripts/bench_harness.py``.

Scoring the real oracle structures needs arda/mmseqs and the manuscript PDBs, so these
tests exercise the LOO *mechanism* in isolation: the per-structure held-out potential
derivation (``_loo_potential_csv``) and the ``n199_r2(loo=True, ...)`` argument
validation. The full-oracle LOO refit is covered by the manuscript sweep, not here.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import pytest

# scripts/ is shipped in the sdist but is not on the import path by default; the harness
# drivers add it the same way.
_SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(_SCRIPTS))

import bench_harness as bh  # noqa: E402
from tcren.potential.model import Potential  # noqa: E402


def _contacts() -> pl.DataFrame:
    # Three structures, each contributing a distinct, recognisable aa-pair signature.
    return pl.DataFrame(
        [{"pdb.id": "s1", "residue.aa.from": "A", "residue.aa.to": "D"}] * 4
        + [{"pdb.id": "s2", "residue.aa.from": "L", "residue.aa.to": "A"}] * 3
        + [{"pdb.id": "s3", "residue.aa.from": "K", "residue.aa.to": "E"}] * 2
    )


def test_loo_potential_excludes_left_out(tmp_path):
    contacts = _contacts()
    ids = ["s1", "s2", "s3"]
    csv = bh._loo_potential_csv(contacts, ids, "s1", tmp_path)
    assert Path(csv).exists()

    # The held-out potential must be loadable and differ from the all-in derivation,
    # because dropping s1 removes all the A->D contacts from the statistics.
    held_out = Potential.from_csv(csv, name="TCRen")
    from tcren.potential import derive_tcren

    full = derive_tcren(contacts, include=ids)
    assert held_out.value("A", "D") != pytest.approx(full.value("A", "D"))
    # And it must equal the explicit include=ids\{s1} derivation byte-for-byte.
    expect = derive_tcren(contacts, include=["s2", "s3"])
    j = held_out.matrix.join(
        expect.matrix, on=["residue.aa.from", "residue.aa.to"], suffix="_e"
    )
    max_abs = j.select((pl.col("value") - pl.col("value_e")).abs().max()).item()
    assert max_abs == pytest.approx(0.0, abs=1e-12)


def test_loo_requires_contacts_and_ids():
    with pytest.raises(ValueError, match="contacts and derivation_ids"):
        bh.n199_r2(loo=True)


def test_non_loo_requires_candidate_csv():
    with pytest.raises(ValueError, match="candidate_csv"):
        bh.n199_r2(loo=False)


def test_loo_forwards_derive_kwargs(tmp_path):
    # derive_kwargs (e.g. weights) must reach derive_tcren through the LOO helper.
    contacts = _contacts()
    ids = ["s1", "s2", "s3"]
    plain_dir = tmp_path / "plain"
    weighted_dir = tmp_path / "weighted"
    plain_dir.mkdir()
    weighted_dir.mkdir()
    base = bh._loo_potential_csv(contacts, ids, "s3", plain_dir)
    weighted = bh._loo_potential_csv(
        contacts, ids, "s3", weighted_dir, weights={"s1": 0.1, "s2": 1.0}
    )
    p_base = Potential.from_csv(base, name="TCRen")
    p_w = Potential.from_csv(weighted, name="TCRen")
    # Down-weighting s1 (the only A->D source left after holding out s3) shifts (A,D).
    assert p_base.value("A", "D") != pytest.approx(p_w.value("A", "D"))
