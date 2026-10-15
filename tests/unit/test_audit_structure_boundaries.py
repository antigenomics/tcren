"""Synthetic regressions for the remaining structure-boundary audit findings."""
import sys
from types import ModuleType

import numpy as np
import pytest

from tcren.structure.model import Atom, Chain, Residue, Structure


def _chain(cid, role, points):
    return Chain(cid, [Residue(i, i + 1, "", "A", "ALA",
                              (Atom("CA", "C", np.asarray(p, dtype=float)),))
                       for i, p in enumerate(points)], chain_type=role)


def test_capped_annotation_propagates_runtime_failure(monkeypatch):
    from tcren.annotation.batch import _annotate

    mapper = ModuleType("arda.annotate.mapper")
    def fail(*args, **kwargs):
        raise RuntimeError("search failed")
    mapper.annotate_records = fail
    monkeypatch.setitem(sys.modules, "arda.annotate.mapper", mapper)
    class Backend:
        def annotate_sequences(self, *args, **kwargs):
            pytest.fail("runtime failure must not retry without the thread cap")
    with pytest.raises(RuntimeError, match="search failed"):
        _annotate(Backend(), [("a", "AAA")], "human", 2)


def test_primary_complex_omits_disconnected_copy():
    from tcren.docking.chains import select_primary_complex

    s = Structure("synthetic", [_chain("P", "PEPTIDE", [[0, 0, 0]]),
                               _chain("D", "MHCa", [[1, 0, 0]]),
                               _chain("A", "TRA", [[2, 0, 0]]),
                               _chain("B1", "TRB", [[100, 0, 0]]),
                               _chain("B2", "TRB", [[200, 0, 0]])])
    assert {c.chain_id for c in select_primary_complex(s).chains} == {"P", "D", "A"}
    assert len(s.chains) == 5


def test_primary_complex_keeps_receptor_connected_through_mhc():
    from tcren.docking.chains import select_primary_complex

    s = Structure("synthetic", [_chain("P", "PEPTIDE", [[0, 0, 0]]),
                               _chain("D", "MHCa", [[6, 0, 0]]),
                               _chain("A", "TRA", [[12, 0, 0]]),
                               _chain("A2", "TRA", [[100, 0, 0]])])
    assert {c.chain_id for c in select_primary_complex(s).chains} == {"P", "D", "A"}


@pytest.mark.parametrize("points", [[[0, 0, 0]] * 3,
                                    [[0, 0, 0], [1, 0, 0], [2, 0, 0]]])
def test_canonical_basis_rejects_degenerate_cloud(points):
    from tcren.docking.frame import _canonical_basis

    with pytest.raises(ValueError, match="degenerate"):
        _canonical_basis(Structure("synthetic", [_chain("A", "TRA", points)]))


def test_canonical_frame_does_not_hide_backend_failure(monkeypatch):
    from tcren.docking import frame

    def fail(*args, **kwargs):
        raise RuntimeError("alignment backend failed")
    monkeypatch.setattr(frame, "align_to_native", fail)
    with pytest.raises(RuntimeError, match="alignment backend failed"):
        frame.canonical_frame(Structure("synthetic", []))


def test_empty_contacts_validate_options():
    from tcren.contacts.geometry import all_atom_contacts, representative_atom_contacts

    s = Structure("empty", [])
    with pytest.raises(ValueError, match="scope"):
        all_atom_contacts(s, scope="invalid")
    with pytest.raises(ValueError, match="kind"):
        representative_atom_contacts(s, kind="invalid")
    for cutoff in (float("nan"), -1):
        with pytest.raises(ValueError, match="cutoff"):
            all_atom_contacts(s, cutoff=cutoff)
        with pytest.raises(ValueError, match="cutoff"):
            representative_atom_contacts(s, cutoff=cutoff)
