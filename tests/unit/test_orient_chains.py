"""Canonical chain renaming + multi-copy detection (synthetic structures, no arda)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.orient.chains import _has_multiple_copies, rename_chains
from tcren.structure.model import Atom, Chain, Residue, Structure


def _chain(chain_id, chain_type, n=2):
    res = [Residue(i, i + 1, "", "A", "ALA", (_a("CA", i),)) for i in range(n)]
    return Chain(chain_id, res, chain_type=chain_type)


def _a(name, i):
    return Atom(name, name[0], np.array([float(i), 0.0, 0.0]))


def _complex(extra=()):
    chains = [_chain("D", "MHCa"), _chain("B2", "B2M"), _chain("P", "PEPTIDE"),
              _chain("A1", "TRA"), _chain("B1", "TRB")]
    chains += list(extra)
    return Structure(pdb_id="syn", chains=chains)


def test_rename_chains_maps_roles_drops_untyped():
    s = _complex(extra=[_chain("X", None)])  # an untyped tag chain
    out, chain_map = rename_chains(s)
    assert [c.chain_id for c in out.chains] == ["A", "B", "C", "D", "E"]  # role order
    assert chain_map == {"A1": "A", "B1": "B", "P": "C", "D": "D", "B2": "E"}
    assert all(c.chain_type is not None for c in out.chains)  # untyped 'X' dropped


def test_has_multiple_copies_and_collision():
    single = _complex()
    assert _has_multiple_copies(single) is False
    two = _complex(extra=[_chain("P2", "PEPTIDE")])  # a second peptide copy
    assert _has_multiple_copies(two) is True
    with pytest.raises(ValueError, match="collision"):
        rename_chains(two)  # two -> role C without select_primary_complex first
