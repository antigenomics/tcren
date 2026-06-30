"""Unit tests for the percentile-rank engine (S3).

This is a new method, so there is no external oracle CSV. The checks are
determinism (same seed => identical result) plus analytic properties on a tiny
hand-built contact map / potential.
"""

from __future__ import annotations

import polars as pl
import pytest

from tcren.contactmap import ContactMap
from tcren.potential import Potential
from tcren.scoring_rank import background_peptides, percentile_rank

_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def _toy_potential() -> Potential:
    # 3-letter alphabet with distinct, hand-checkable values (mirrors test_scoring).
    vals = {("A", "A"): 1.0, ("A", "K"): -2.0, ("L", "A"): 0.5, ("L", "K"): 3.0,
            ("A", "G"): 0.1, ("L", "G"): 0.2}
    rows = [{"residue.aa.from": fr, "residue.aa.to": to, "value": v}
            for (fr, to), v in vals.items()]
    return Potential(name="toy", matrix=pl.DataFrame(rows), alphabet=("A", "L", "K", "G"))


def _toy_contact_map() -> ContactMap:
    # TCR 'A' contacts peptide pos 0; TCR 'L' contacts peptide pos 2.
    contacts = pl.DataFrame(
        {
            "chain.type.from": ["TRA", "TRB"],
            "chain.type.to": ["PEPTIDE", "PEPTIDE"],
            "residue.aa.from": ["A", "L"],
            "residue.aa.to": ["G", "G"],
            "region.type.from": ["CDR3", "CDR3"],
            "residue.index.from": [10, 20],
            "residue.index.to": [0, 2],
            "region.start.from": [8, 18],
            "region.start.to": [0, 0],
            "pdb.id": ["toy", "toy"],
        }
    )
    return ContactMap(pdb_id="toy", contacts=contacts, peptide_length=3)


# --- background_peptides --------------------------------------------------------------

def test_background_determinism():
    a = background_peptides(9, n=50, seed=42)
    b = background_peptides(9, n=50, seed=42)
    assert a == b
    assert background_peptides(9, n=50, seed=43) != a


def test_background_shape_and_alphabet():
    bg = background_peptides(9, n=100, seed=0)
    assert len(bg) == 100
    assert all(len(p) == 9 for p in bg)
    assert set("".join(bg)) <= set(_AMINO_ACIDS)


def test_background_from_source(tmp_path):
    fasta = tmp_path / "epi.fasta"
    fasta.write_text(">e1\nAGKAGKAGK\n>e2\nKKAAGGLLM\n>short\nAGK\n")
    bg = background_peptides(9, n=20, seed=1, source=str(fasta))
    assert len(bg) == 20
    assert all(len(p) == 9 for p in bg)
    # permutations of length-9 epitopes only; the length-3 entry is excluded.
    pool = {"AGKAGKAGK", "KKAAGGLLM"}
    assert all(sorted(p) in [sorted(s) for s in pool] for p in bg)


# --- percentile_rank ------------------------------------------------------------------

def test_rank_determinism():
    cm, pot = _toy_contact_map(), _toy_potential()
    r1 = percentile_rank(cm, "AGK", pot, n_background=200, seed=7)
    r2 = percentile_rank(cm, "AGK", pot, n_background=200, seed=7)
    assert r1 == r2


def test_rank_properties():
    cm, pot = _toy_contact_map(), _toy_potential()
    res = percentile_rank(cm, "AGK", pot, n_background=200, seed=0)
    assert set(res) == {"peptide", "score", "rank_pct", "n_background"}
    assert res["peptide"] == "AGK"
    assert 0.0 <= res["rank_pct"] <= 1.0
    assert res["n_background"] == 200
    # native "AGK": pos0 'A' x TCR 'A' = (A,A)=1.0 ; pos2 'K' x TCR 'L' = (L,K)=3.0 -> 4.0
    assert res["score"] == pytest.approx(4.0)


def test_rank_analytic_with_explicit_background():
    cm, pot = _toy_contact_map(), _toy_potential()
    # Native "AGK" -> 4.0. Build a hand-chosen background:
    #   "AGA": (A,A)=1.0 + (L,A)=0.5 = 1.5   (< 4.0  -> counts)
    #   "AGK": 4.0                            (== 4.0 -> tie, counts)
    #   "LGK": (A,L) missing -> dropped; (L,K)=3.0 -> 3.0  (< 4.0 -> counts)
    #          ('L' at pos0 pairs TCR 'A'; (A,L) absent from potential, contributes 0)
    # All three have score <= 4.0 => rank_pct == 1.0
    bg = ["AGA", "AGK", "LGK"]
    res = percentile_rank(cm, "AGK", pot, background=bg)
    assert res["n_background"] == 3
    assert res["rank_pct"] == pytest.approx(1.0)

    # A background strictly worse than native gives rank_pct 0.
    # "GGG": pos0 'G' x TCR 'A' = (A,G)=0.1 ; pos2 'G' x TCR 'L' = (L,G)=0.2 -> 0.3 (better!)
    # Use a peptide that scores worse: need score > 4.0. "AGK" axis only reaches 4.0 max here,
    # so instead rank the weakest binder against a strong-only background.
    res2 = percentile_rank(cm, "GGG", pot, background=["AGK"])
    # native "GGG" -> 0.3 ; background "AGK" -> 4.0 (> 0.3, not <=) => rank_pct 0
    assert res2["score"] == pytest.approx(0.3)
    assert res2["rank_pct"] == pytest.approx(0.0)


def test_rank_native_missing_raises():
    cm, pot = _toy_contact_map(), _toy_potential()
    with pytest.raises(ValueError):
        # wrong length -> dropped by score_peptides' length filter
        percentile_rank(cm, "AGKK", pot, background=["AGK"])
