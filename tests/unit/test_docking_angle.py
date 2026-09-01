"""Unit tests for the TCR docking-angle geometry (pure trig, no structure/mmseqs)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.docking import crossing_incident_from_vector


@pytest.mark.parametrize(
    "v,exp_crossing,exp_signed,exp_incident",
    [
        ([0.0, 10.0, 0.0], 0.0, 0.0, 0.0),     # Vα→Vβ along the groove long axis
        ([10.0, 0.0, 0.0], 90.0, 90.0, 0.0),   # perpendicular, in plane
        ([-10.0, 0.0, 0.0], 90.0, -90.0, 0.0),  # perpendicular, other handedness
        ([10.0, 10.0, 0.0], 45.0, 45.0, 0.0),  # 45° crossing
        ([0.0, 10.0, 10.0], 0.0, 0.0, 45.0),   # along groove but tilted +45° out of plane
        ([0.0, 10.0, -10.0], 0.0, 0.0, -45.0),  # tilted down
    ],
)
def test_crossing_incident_from_vector(v, exp_crossing, exp_signed, exp_incident):
    crossing, signed, incident = crossing_incident_from_vector(np.asarray(v))
    assert crossing == pytest.approx(exp_crossing, abs=1e-6)
    assert signed == pytest.approx(exp_signed, abs=1e-6)
    assert incident == pytest.approx(exp_incident, abs=1e-6)


def test_crossing_undefined_when_normal_to_plane():
    with pytest.raises(ValueError, match="normal to the groove plane"):
        crossing_incident_from_vector(np.asarray([0.0, 0.0, 10.0]))


# --- TCR placement (translations, not rotations) ---------------------------------------------

from tcren.docking import tcr_placement  # noqa: E402
from tcren.docking.angles import _groove_frame  # noqa: E402
from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure  # noqa: E402


def _res(i, xyz):
    return Residue(i, i + 1, "", "G", "GLY", (Atom("CA", "C", np.asarray(xyz, float)),))


def _placement_complex(cdr_xyz):
    """Peptide along +y at the origin; a TCR whose CDR loop sits at ``cdr_xyz``."""
    pep = Chain("C", [_res(i, [0.0, float(i) - 1.0, 0.0]) for i in range(3)],
                chain_type="PEPTIDE")
    # Framework residues place the TCR centroid straight above the peptide; the CDR loop is a
    # separate, offsettable landmark. Many framework residues, because the groove normal is fit
    # from the TCR centroid itself -- an off-axis loop tilts the frame in proportion to its weight
    # in that centroid, and with only a handful of residues the tilt is several degrees.
    n_fr = 39
    fr = [_res(i, [0.0, 0.0, 20.0 + 0.1 * i]) for i in range(n_fr)]
    cdr = [_res(n_fr, list(cdr_xyz))]
    tcr = Chain("B", fr + cdr, chain_type="TRB")
    tcr.regions = [RegionMarkup("CDR3", n_fr, n_fr, "G", cdr)]
    return Structure("synth", [pep, tcr])


def test_placement_height_is_a_distance_above_the_groove_plane():
    s = _placement_complex([0.0, 0.0, 12.0])
    p = tcr_placement(s)
    assert p.height == pytest.approx(12.0, abs=1e-6)
    assert p.offset == pytest.approx(0.0, abs=1e-6)
    assert p.n_cdr == 1


def test_placement_shift_u_follows_the_peptide_long_axis():
    # peptide runs N->C along +y, so a loop displaced along +y is shifted toward the C-terminus
    p = tcr_placement(_placement_complex([0.0, 4.0, 12.0]))
    assert p.shift_u == pytest.approx(4.0, abs=1e-6)
    assert p.shift_w == pytest.approx(0.0, abs=1e-6)
    assert p.offset == pytest.approx(4.0, abs=1e-6)
    # height is measured to the loop, not to the whole domain sitting at z ~ 20
    assert p.height == pytest.approx(12.0, abs=0.05)


def test_placement_shift_w_is_the_across_groove_displacement():
    p = tcr_placement(_placement_complex([3.0, 0.0, 12.0]))
    assert abs(p.shift_w) == pytest.approx(3.0, abs=0.05)
    assert p.shift_u == pytest.approx(0.0, abs=1e-6)


def test_whole_tcr_centroid_would_have_no_lateral_component_by_construction():
    # This is why the footprint centroid is the reference point: the groove frame builds w
    # perpendicular to the peptide->TCR-centroid vector, so that vector's w component is 0 exactly.
    s = _placement_complex([3.0, 4.0, 12.0])
    u, w, n = _groove_frame(s)
    from tcren.docking.angles import _TCR_TYPES, _chain_ca
    pep, tcr = _chain_ca(s, ("PEPTIDE",)), _chain_ca(s, _TCR_TYPES)
    d = tcr.mean(axis=0) - pep.mean(axis=0)
    assert float(np.dot(d, w)) == pytest.approx(0.0, abs=1e-9)
    # ...while the CDR footprint does carry one
    assert abs(tcr_placement(s).shift_w) > 1.0


def test_placement_needs_region_markup():
    s = _placement_complex([0.0, 0.0, 12.0])
    s.chains[1].regions = []
    with pytest.raises(ValueError, match="region markup"):
        tcr_placement(s)
