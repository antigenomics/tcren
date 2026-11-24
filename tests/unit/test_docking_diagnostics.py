"""Pose diagnostics on explicit contacts and an MHC-fixed groove plane."""

import numpy as np
import pytest

from tcren.docking.diagnostics import diagnose_docking
from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure


def _pose(a, b, *, floor=True, rotate=False):
    def chain(chain_id, role, points):
        coords = np.asarray(points, dtype=float)
        if rotate:
            coords = coords[:, [2, 0, 1]] + [7, -3, 11]
        residues = [
            Residue(i, i + 1, "", "A", "ALA", (Atom("CA", "C", point),))
            for i, point in enumerate(coords)
        ]
        return Chain(chain_id, residues, role)

    mhc = chain("M", "MHCa", [(-8, -9, 0), (-8, 9, 0), (8, -9, 0), (8, 9, 0), (8, 0, 0)])
    if floor:
        mhc.regions = [RegionMarkup("GROOVE_FLOOR", 0, 3, "AAAA", mhc.residues[:4])]
    return Structure("pose", [
        mhc, chain("P", "PEPTIDE", [(-3, 0, 4), (0, 0, 4), (3, 0, 4)]),
        chain("A", "TRA", [a]), chain("B", "TRB", [b]),
    ])


def test_docking_flags_are_independent_and_rigid_transform_invariant():
    above = diagnose_docking(_pose((-2, 0, 7), (2, 0, 7)))
    assert above.reason_codes == ()
    assert above.peptide_contacts == {"A": 2, "B": 2}
    assert above.side_angle_deg == pytest.approx(0)

    one_chain = diagnose_docking(_pose((-2, 0, 7), (8, 0, 3)))
    assert one_chain.peptide_contacts == {"A": 2, "B": 0}
    assert one_chain.pmhc_contacts == {"A": 2, "B": 1}
    assert one_chain.reason_codes == ("chain_no_peptide_contact",)

    one_chain_detached = diagnose_docking(_pose((-2, 0, 7), (12, 0, 15)))
    assert one_chain_detached.pmhc_contacts == {"A": 2, "B": 0}
    assert one_chain_detached.reason_codes == (
        "chain_no_peptide_contact", "chain_no_pmhc_contact",
    )

    for rotate in (False, True):
        lateral = diagnose_docking(_pose((-7, 9, 1), (7, 9, 1), rotate=rotate))
        assert lateral.peptide_contacts == {"A": 0, "B": 0}
        assert lateral.pmhc_contacts == {"A": 1, "B": 1}
        assert lateral.reason_codes == ("no_peptide_contact", "side_docking")
        assert lateral.side_angle_deg == pytest.approx(np.degrees(np.arctan2(9, -3)))

    detached = diagnose_docking(_pose((-2, 0, 15), (2, 0, 15)))
    assert detached.reason_codes == ("no_peptide_contact", "no_pmhc_contact")


def test_missing_floor_preserves_contact_diagnostics():
    result = diagnose_docking(_pose((-2, 0, 7), (2, 0, 7), floor=False))
    assert result.peptide_contacts == {"A": 2, "B": 2}
    assert result.side_angle_deg is None
    assert "groove_geometry_unavailable" in result.reason_codes


def test_missing_receptor_chain_is_an_invalid_input():
    pose = _pose((-2, 0, 7), (2, 0, 7))
    pose.chains.pop()
    with pytest.raises(ValueError, match="receptor"):
        diagnose_docking(pose)
