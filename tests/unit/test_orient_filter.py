"""Unit tests for the oriented-complex sanity filter (pure geometry, no arda)."""

from __future__ import annotations

import numpy as np

from tcren.orient import check_oriented_complex
from tcren.structure.model import Atom, Chain, Residue, Structure


def _chain(cid, centre, n=8):
    res = [Residue(i, i, "", "A", "ALA", (Atom("CA", "C", np.asarray(centre, float) + [i, 0, 0]),))
           for i in range(n)]
    return Chain(cid, res)


def _complex(pep_centre=(0, 0, 0), tcr_centre=(0, 0, 10), e_centre=(0, 0, -42), pep_n=9):
    # A/B TCR on top (+z), C peptide near origin, D MHC below, E b2m further below
    return Structure("t", [
        _chain("A", tcr_centre), _chain("B", tcr_centre),
        _chain("C", pep_centre, n=pep_n),
        _chain("D", (0, 0, -25)), _chain("E", e_centre),
    ])


def test_sane_complex_passes():
    ok, reason = check_oriented_complex(_complex())
    assert ok and reason == "ok"


def test_peptide_off_center_rejected():
    ok, reason = check_oriented_complex(_complex(pep_centre=(40, 0, 0)))
    assert not ok and reason == "peptide_off_center"


def test_peptide_too_long_rejected():
    ok, reason = check_oriented_complex(_complex(pep_n=28))
    assert not ok and reason == "peptide_too_long"


def test_orphan_chain_rejected():
    ok, reason = check_oriented_complex(_complex(e_centre=(120, 0, 0)))
    assert not ok and reason.startswith("orphan_chain")


def test_tcr_not_engaged_rejected():
    ok, reason = check_oriented_complex(_complex(tcr_centre=(0, 0, 60)))
    assert not ok and reason == "tcr_not_engaged"
