"""Interface-sanity filter: flag trivial non-interfaces as assay noise before scoring.

A TCR-pMHC "complex" can fail to be a real interface for two reasons that need no energy
model to spot: the dock collapsed (few or no TCR:peptide contacts) or the docking geometry
is out of the physiological range (a failed pose lands the TCR at an implausible crossing or
tilt). :func:`is_real_interface` rejects both from three cheap descriptors so the downstream
energy score is only ever asked about interfaces that could plausibly be binders.

The thresholds are the p01-p99 range of the 309 real (binder) interfaces from the TCRvdb
AlphaFold set, rounded *inward* (lower bounds up, upper bounds down): at least
:data:`N_CONTACTS_MIN` TCR:peptide residue contacts, a crossing (scanning) angle in
:data:`SCANNING_RANGE` and an incident (pitch) angle in :data:`PITCH_RANGE`. Any missing
descriptor (``NaN``/``None``, i.e. an undocked or un-oriented complex) is treated as noise.

The derivation is ``bench/fit_models.py::envelope`` in the benchmark repo, which reproduces all
three constants exactly and is regression-tested there. Two caveats it records: the pitch *lower*
bound of 0 deg is the domain floor of an unsigned angle rather than a percentile (p01 = 0.04), and
the pitch axis is derived from a cached ``pitch_angle`` column carrying AlphaFold-confidence
leakage -- the scanning and contact-count axes are clean.
"""

from __future__ import annotations

import math

N_CONTACTS_MIN = 15
SCANNING_RANGE = (24.0, 70.0)  # crossing_angle (== scanning_angle), degrees
PITCH_RANGE = (0.0, 14.0)  # incident_angle (== pitch_angle), degrees


def _is_missing(x: float | None) -> bool:
    """True if ``x`` is ``None`` or a NaN float (an undocked / un-oriented descriptor)."""
    if x is None:
        return True
    try:
        return math.isnan(x)
    except (TypeError, ValueError):
        return True


def is_real_interface(
    n_contacts: float | None,
    scanning_angle: float | None,
    pitch_angle: float | None,
) -> bool:
    """False if the TCR:pMHC interface is assay noise / a failed dock.

    ``NaN`` (or ``None``) in any input => ``False`` (undocked). Thresholds are the ~p01-p99
    range of real (binder) interfaces from the TCRvdb AF set.

    Args:
        n_contacts: number of TCR:peptide residue-pair contacts at the interface.
        scanning_angle: TCR crossing (scanning) angle in degrees
            (``DockingAngles.crossing_angle``).
        pitch_angle: TCR incident (pitch) angle in degrees
            (``DockingAngles.incident_angle``).

    Returns:
        ``True`` only if ``n_contacts >= N_CONTACTS_MIN`` and ``scanning_angle`` lies in
        :data:`SCANNING_RANGE` and ``pitch_angle`` lies in :data:`PITCH_RANGE`; ``False``
        otherwise, including when any input is missing.

    Example:
        >>> is_real_interface(25, 45.0, 5.0)
        True
        >>> is_real_interface(0, 45.0, 5.0)
        False
    """
    if _is_missing(n_contacts) or _is_missing(scanning_angle) or _is_missing(pitch_angle):
        return False
    return (
        n_contacts >= N_CONTACTS_MIN
        and SCANNING_RANGE[0] <= scanning_angle <= SCANNING_RANGE[1]
        and PITCH_RANGE[0] <= pitch_angle <= PITCH_RANGE[1]
    )
