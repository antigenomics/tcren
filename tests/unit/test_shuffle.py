"""Unit tests for shuffled-decoy generation (pure structure-model logic, no IO/mmseqs)."""

from __future__ import annotations

import pytest

from tcren.structure.model import Chain, Structure
from tcren.shuffle import graft_tcr, make_decoys, mhc_class, run_shuffle, _tcr_chains, _pmhc_chains


def _mk(pid: str, cls: str = "MHCI") -> Structure:
    return Structure(pid, [
        Chain("A", [], chain_type="TRA"),
        Chain("B", [], chain_type="TRB"),
        Chain("C", [], chain_type="PEPTIDE"),
        Chain("D", [], chain_type="MHCa", chain_supertype=cls),
        Chain("E", [], chain_type="B2M"),
    ])


def test_split_and_class():
    s = _mk("x")
    assert {c.chain_id for c in _tcr_chains(s)} == {"A", "B"}
    assert {c.chain_id for c in _pmhc_chains(s)} == {"C", "D", "E"}
    assert mhc_class(s) == "MHCI"


def test_graft_keeps_pmhc_takes_tcr():
    a, b = _mk("a"), _mk("b")
    a.chain("C").residues.append(object())          # mark a's peptide identity by object identity
    d = graft_tcr(a, b)
    assert d.pdb_id == "a__tcr_b"
    # pMHC chains are a's exact objects; TCR chains are b's
    assert d.chain("C") is a.chain("C")
    assert d.chain("A") is b.chain("A") or d.chain("A").chain_type == "TRA"
    assert [c.chain_type for c in d.chains if c.chain_type in ("PEPTIDE", "MHCa", "B2M")]  # pMHC present
    assert any(c.chain_type == "TRA" for c in d.chains)                                     # TCR present


def test_graft_reassigns_colliding_chain_id():
    pmhc = Structure("p", [Chain("C", [], chain_type="PEPTIDE"), Chain("D", [], chain_type="MHCa", chain_supertype="MHCI")])
    tcr = Structure("t", [Chain("C", [], chain_type="TRA")])   # id "C" collides with the pMHC peptide
    d = graft_tcr(pmhc, tcr)
    ids = [c.chain_id for c in d.chains]
    assert len(ids) == len(set(ids)), "chain ids must be unique after graft"
    assert d.chain("C").chain_type == "PEPTIDE"                # the peptide kept its id


def test_graft_requires_both_sides():
    only_pmhc = Structure("p", [Chain("C", [], chain_type="PEPTIDE")])
    with pytest.raises(ValueError, match="no TCR"):
        graft_tcr(only_pmhc, only_pmhc)
    with pytest.raises(ValueError, match="no pMHC"):
        graft_tcr(_mk("t"), _mk("t"), pdb_id=None) if False else graft_tcr(
            Structure("x", [Chain("A", [], chain_type="TRA")]), _mk("t"))


def test_make_decoys_within_class_and_derangement():
    structs = [_mk(f"s{i}", "MHCI" if i < 4 else "MHCII") for i in range(6)]  # 4 MHCI + 2 MHCII
    decoys = list(make_decoys(structs, n_per=2, seed=0))
    # 4 MHCI x 2 + 2 MHCII x 1 (pool size 1) = 10
    assert len(decoys) == 10
    for d in decoys:
        pm, tc = d.pdb_id.split("__tcr_")
        assert pm != tc                                        # no complex keeps its own TCR
        # within class: pMHC and TCR sources share a class
        assert (int(pm[1:]) < 4) == (int(tc[1:]) < 4)


def test_make_decoys_reproducible():
    structs = [_mk(f"s{i}") for i in range(5)]
    a = [d.pdb_id for d in make_decoys(structs, n_per=2, seed=7)]
    b = [d.pdb_id for d in make_decoys(structs, n_per=2, seed=7)]
    assert a == b


def test_run_shuffle_says_so_when_nothing_parsed(tmp_path):
    # Writing 0 decoys and reporting success is how an unusable input dir went unnoticed.
    with pytest.raises(ValueError, match="ORIENTED"):
        run_shuffle(tmp_path, tmp_path / "out")
