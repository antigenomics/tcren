"""Fast tests for the interface-sanity (assay-noise) filter."""

from __future__ import annotations

import math

from tcren.binder.noise import is_real_interface


def test_good_interface_passes():
    assert is_real_interface(25, 45.0, 5.0) is True


def test_no_contacts_fails():
    assert is_real_interface(0, 45.0, 5.0) is False


def test_nan_angle_fails():
    assert is_real_interface(25, math.nan, 5.0) is False


def test_none_input_fails():
    assert is_real_interface(25, None, 5.0) is False


def test_out_of_range_scanning_fails():
    assert is_real_interface(25, 90.0, 5.0) is False
