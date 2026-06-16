"""Unit tests for the TCR docking-angle geometry (pure trig, no structure/mmseqs)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.orient import crossing_incident_from_vector


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
