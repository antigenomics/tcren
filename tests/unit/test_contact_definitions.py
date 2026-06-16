"""Unit tests for the flexible d1/d2/d3 contact definition."""

from __future__ import annotations

import numpy as np

from tcren.contacts import multi_contacts, representative_atom_contacts
from tcren.structure.model import Atom, Chain, Residue, Structure


def _res(seq, aa, ca, cb=None):
    atoms = [Atom("CA", "C", np.asarray(ca, dtype=float))]
    if cb is not None:
        atoms.append(Atom("CB", "C", np.asarray(cb, dtype=float)))
    resname = "GLY" if aa == "G" else "ALA"
    return Residue(seq, seq, "", aa, resname, tuple(atoms))


def _structure():
    # A:0 Gly (CA only), A:1 Ala (CA+CB); B:0 Ala (CA+CB) placed so CB-CB is close but
    # CA-CA is farther — exercises the Cβ layer and the glycine Cα fallback.
    a = Chain("A", [_res(0, "G", ca=(0, 0, 0)), _res(1, "A", ca=(0, 10, 0), cb=(1, 10, 0))])
    b = Chain("B", [_res(0, "A", ca=(7, 0, 0), cb=(3, 0, 0))])
    return Structure("toy", [a, b])


def test_cb_or_ca_glycine_fallback():
    s = _structure()
    gly, ala = s.chains[0].residues
    assert gly.cb is None and gly.cb_or_ca is gly.ca  # glycine: no CB, falls back to CA
    assert ala.cb is not None and np.allclose(ala.cb_or_ca, ala.cb)


def test_cb_layer_uses_cb_with_glycine_ca_fallback():
    s = _structure()
    # Cβ contacts within 8 A: A:0(Gly→CA at x=0) vs B:0(CB at x=3) = 3 A -> contact;
    # A:1(CB at (1,10)) vs B:0(CB at (3,0)) ~ 10.2 A -> no contact.
    cb = representative_atom_contacts(s, kind="cb", cutoff=8.0)
    pairs = {(r["chain.id.from"], r["residue.index.from"], r["chain.id.to"], r["residue.index.to"])
             for r in cb.iter_rows(named=True)}
    assert ("A", 0, "B", 0) in pairs
    assert cb["atom.from"].to_list() and set(cb["atom.from"].to_list()) == {"CB"}


def test_ca_layer_distance():
    s = _structure()
    ca = representative_atom_contacts(s, kind="ca", cutoff=12.0)
    hit = ca.filter((ca["chain.id.from"] == "A") & (ca["residue.index.from"] == 0)
                    & (ca["chain.id.to"] == "B") & (ca["residue.index.to"] == 0))
    assert hit.height == 1 and abs(hit["dist"][0] - 7.0) < 1e-9  # CA-CA = 7 A


def test_multi_contacts_has_three_layers():
    s = _structure()
    layers = set(multi_contacts(s)["layer"].to_list())
    assert layers <= {"d1", "d2", "d3"} and "d3" in layers
