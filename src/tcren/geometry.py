"""Closed-form CDR3 Ω-loop reachability — the feasibility filter.

The CDR3 loop is a chain of ``n = L + 1`` virtual Cα–Cα bonds of length :data:`BOND_CA`, pinned at
its neck separation ``R`` (the ``cdr3{a,b}_ext`` frame descriptor). Pure geometry then bounds how far
the loop's apex can reach from the neck axis, independent of sequence or energy:

* :func:`reach_max` — the maximum apex reach for a loop of length ``L`` pinned at ``R``.
* :func:`reachability_floor` — the shortest loop that can span a target at distance ``d``.
* :func:`span_saturation` — the fraction of maximal reach a structure actually uses,
  ``reach / reach_max``; ``> 1`` is geometrically impossible and real repertoires sit well below 1
  (~0.3–0.5 at CDR3 length 10–20), i.e. with large conformational slack.

These are **feasibility / mode** descriptors, not binder classifiers — backbone Cα geometry does not
discriminate binders (every apparent signal was a provenance / epitope confound). Bill them as
``structure -> geometry`` featurisation only.

Ported from ``model/src/Feasibility.jl`` of the 2026-tcren-algem monograph.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = ["BOND_CA", "reach_max", "reachability_floor", "span_saturation",
           "cdr3_internal_coords", "LoopInternalCoords"]

#: Virtual Cα–Cα bond length (Å) used for the loop chain.
BOND_CA = 3.80


def reach_max(length: int, neck: float, bond: float = BOND_CA) -> float:
    """Maximum apex reach ``sqrt(((L+1)b/2)^2 - (R/2)^2)`` of a length-``length`` loop pinned at ``neck``.

    Args:
        length: CDR3 loop length ``L`` (number of residues).
        neck: neck separation ``R`` in Å (the ``ext`` frame descriptor).
        bond: virtual Cα–Cα bond length (default :data:`BOND_CA`).

    Returns:
        Maximum reach in Å, or ``0.0`` if the loop is too short to span the neck.
    """
    half = (length + 1) * bond / 2.0
    return math.sqrt(half * half - (neck / 2.0) ** 2) if half > neck / 2.0 else 0.0


def reachability_floor(distance: float, neck: float, bond: float = BOND_CA) -> int:
    """Shortest loop length that can reach a target at ``distance`` from a neck of separation ``neck``.

    Inverts :func:`reach_max`: the minimum ``L`` with ``reach_max(L, neck) >= distance``.
    """
    return int(math.ceil((2.0 / bond) * math.sqrt(distance * distance + (neck / 2.0) ** 2) - 1.0))


def span_saturation(structure) -> dict[str, float]:
    """Per-loop ``reach / reach_max`` for ``cdr3a`` and ``cdr3b`` of a chain-typed structure.

    Composes the existing CDR3 frame descriptors: ``reach`` (loop-centroid distance from the groove
    origin) and ``ext`` (neck separation ``R``), with ``L`` the resolved CDR3 length. ``> 1`` flags an
    infeasible pose (the loop cannot reach that far); real structures sit ~0.3–0.5. Returns ``NaN`` for
    a loop whose frame or CDR3 span is undefined.

    Note:
        The structure must be chain-typed (``classify_chains``) so its CDR3 regions are populated.
    """
    from .recognition import _cdr3_frame_features

    frame = _cdr3_frame_features(structure)
    lengths = _cdr3_lengths(structure)
    out: dict[str, float] = {}
    for loop in ("cdr3a", "cdr3b"):
        reach, neck, length = frame.get(f"{loop}_reach"), frame.get(f"{loop}_ext"), lengths.get(loop)
        rmax = reach_max(length, neck) if (length and neck is not None and math.isfinite(neck)) else 0.0
        out[loop] = reach / rmax if (rmax > 0 and reach is not None and math.isfinite(reach)) else math.nan
    return out


def _cdr3_lengths(structure) -> dict[str, int]:
    """Resolved CDR3 Cα count per loop (``cdr3a``/``cdr3b``) from the region markup."""
    out: dict[str, int] = {}
    for loop, ctype in (("cdr3a", "TRA"), ("cdr3b", "TRB")):
        for c in structure.chains:
            if getattr(c, "chain_type", None) != ctype:
                continue
            for reg in getattr(c, "regions", []) or []:
                if reg.region_type == "CDR3":
                    out[loop] = sum(1 for r in reg.residues if r.ca is not None)
                    break
    return out


@dataclass(slots=True)
class LoopInternalCoords:
    """Virtual-bond internal coordinates of a CDR3 Cα trace (the input to the Julia loop model).

    A Cα chain of ``N`` points has ``N-1`` bonds, ``N-2`` bond angles (interior vertices) and ``N-3``
    pseudo-torsions — the inverse of a discrete-Frenet reconstruction. ``neck`` is the end-to-end
    span ``|p_last - p_first|`` (the ``ext`` frame descriptor). Angles and torsions are in **radians**.
    """

    bonds: object            #: (N-1,) virtual Cα–Cα bond lengths, Å
    angles: object           #: (N-2,) Cα–Cα–Cα pseudo bond angles θ, radians
    torsions: object         #: (N-3,) Cα pseudo-dihedrals τ, radians in (−π, π]
    neck: float              #: end-to-end span R = |p_last − p_first|, Å
    n_ca: int


def _bond_angle(p1, p2, p3) -> float:
    """Cα–Cα–Cα pseudo bond angle at ``p2`` (radians)."""
    a = p1 - p2
    b = p3 - p2
    c = float(a.dot(b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))
    return math.acos(max(-1.0, min(1.0, c)))


def cdr3_internal_coords(structure, loop: str = "cdr3b") -> LoopInternalCoords | None:
    """Virtual-bond internal coordinates (θ, τ, neck R) of a CDR3 loop's Cα trace.

    A thin ``structure -> coordinates`` extractor: the loop's *model* (torsion sampling, closure,
    ensemble entropy) lives in the algem monograph's Julia code and does not migrate — tcren only
    hands off the intrinsic internal coordinates of the observed loop. Reuses the CDR3 region markup
    (``classify_chains``) for the span and the shared ``_dihedral`` primitive for torsions.

    Args:
        structure: a chain-typed structure.
        loop: ``"cdr3a"`` (TRA) or ``"cdr3b"`` (TRB).

    Returns:
        A :class:`LoopInternalCoords`, or ``None`` if the loop has fewer than 4 resolved Cα (a
        torsion needs 4 points).
    """
    from .contactmap import _cdr3_ca
    from .orient.tcrdock_geometry import _dihedral

    ctype = {"cdr3a": "TRA", "cdr3b": "TRB"}.get(loop)
    if ctype is None:
        raise ValueError(f"loop must be 'cdr3a' or 'cdr3b', got {loop!r}")
    ca = _cdr3_ca(structure, ctype)
    if ca is None or len(ca) < 4:
        return None
    bonds = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    angles = np.array([_bond_angle(ca[k], ca[k + 1], ca[k + 2]) for k in range(len(ca) - 2)])
    torsions = np.array([_dihedral(ca[k], ca[k + 1], ca[k + 2], ca[k + 3])
                         for k in range(len(ca) - 3)])
    neck = float(np.linalg.norm(ca[-1] - ca[0]))
    return LoopInternalCoords(bonds=bonds, angles=angles, torsions=torsions, neck=neck, n_ca=len(ca))
