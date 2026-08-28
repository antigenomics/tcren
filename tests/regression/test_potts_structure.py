"""The bridge between a real Structure and the Potts stack. 2026-08-28

Every test in ``tests/unit/test_potts.py`` starts from a hand-built site frame, so the one function
that actually reads coordinates -- :func:`tcren.potts.available_pairs` -- and the one-shot wrapper
over it were the only parts of the module with no coverage at all. They are also the parts a
structure-format change breaks first.
"""

from __future__ import annotations

from pathlib import Path

import pytest

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"
pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def structure():
    """Chain-typed, exactly as ``tcren potts`` does it: parse, then classify."""
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.structure import parse_structure

    s = parse_structure(ASSET, pdb_id="1ao7")
    classify_chains(s, organism="human", autodetect_species=True)
    return s


def test_available_pairs_reads_a_structure_into_the_site_schema(structure):
    from tcren.potts import available_pairs

    p = available_pairs(structure)
    assert p.height > 0
    assert {"pdb.id", "aa.rec", "chain.rec", "region.rec", "pos.rec", "aa.par", "pos.par",
            "role.par", "cls", "d_ca", "sigma"} <= set(p.columns), p.columns
    assert p["d_ca"].max() <= 15.0                       # the availability radius
    assert set(p["sigma"].unique()) <= {0.0, 1.0}
    assert 0.0 < p["sigma"].mean() < 1.0, "a real interface has both contacts and near-misses"
    assert p["pos.par"].min() == 0, "pos.par is a within-region offset, 0-based"


def test_available_pairs_rejects_an_unknown_partner(structure):
    from tcren.potts import available_pairs

    with pytest.raises(ValueError, match="partner"):
        available_pairs(structure, partner="nonsense")


def test_score_structure_is_the_one_shot_path(structure):
    """The convenience wrapper has to agree with enumerate-then-score, or it is a second answer."""
    import polars as pl

    from tcren.potts import PottsModel, available_pairs, score_sites, score_structure

    m = PottsModel.bundled()
    one = score_structure(structure, model=m, seed=0)
    ref = score_sites(available_pairs(structure), m, seed=0).to_dicts()[0]
    assert isinstance(one, dict)
    assert one["n_sites"] == ref["n_sites"] and one["n_contacts"] == ref["n_contacts"]
    assert one["energy"] == pytest.approx(ref["energy"])
    assert one["neg_energy"] == pytest.approx(-one["energy"])
    assert isinstance(pl.DataFrame([one]), pl.DataFrame)
