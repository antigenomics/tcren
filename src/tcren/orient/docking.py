"""TCR docking geometry: crossing angle and incident (tilt) angle.

The TCR "docking angle" (Rudolph, Stanfield & Wilson 2006; Garcia et al.) describes how the
αβ (or γδ) TCR sits on top of the peptide-MHC groove. It is computed here directly from the
canonical frame (:func:`tcren.orient.frame.canonical_frame`), so no external package
(TCRdock / STCRpy) is required:

* the **crossing angle** is the angle between the Vα→Vβ pseudo-axis projected into the MHC
  groove plane and the groove long axis (peptide N→C, canonical ``+y`` — collinear with the
  MHC α1 helix to within a few degrees). Reported on ``[0, 180)``; canonical αβ TCRs cluster
  around ~20–70°. A signed variant (``[-180, 180)``) carries the handedness of the docking.
* the **incident (tilt) angle** is the elevation of the same Vα→Vβ vector out of the groove
  plane (canonical ``z`` is the MHC→TCR normal): positive when Vβ rides higher above the groove
  than Vα.

The Vα/Vβ landmarks are the centroids of the variable-domain Cα atoms of the two receptor
chains (TRA/TRB for αβ, TRG/TRD for γδ). The frame is fit from the query itself (PCA), so the
calculation needs neither the native database nor mmseqs once the structure is chain-typed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..structure.model import Structure

# Receptor chain pairing: (V-alpha-like, V-beta-like) for each cell type.
_AB_PAIR = ("TRA", "TRB")
_GD_PAIR = ("TRG", "TRD")
_TCR_TYPES = _AB_PAIR + _GD_PAIR


@dataclass(slots=True)
class DockingAngles:
    """TCR docking geometry relative to the MHC groove (all angles in degrees)."""

    crossing_angle: float  # [0, 180): Vα→Vβ vs groove long axis, in the groove plane
    crossing_angle_signed: float  # [-180, 180): same, carrying docking handedness
    incident_angle: float  # [-90, 90]: elevation of Vα→Vβ out of the groove plane
    cell_type: str  # "ab" | "gd"
    n_va: int  # Vα(-like) Cα atoms used
    n_vb: int  # Vβ(-like) Cα atoms used


def _domain_centroid(structure: Structure, chain_type: str) -> np.ndarray | None:
    """Centroid of the variable-domain Cα atoms of the first chain of ``chain_type``."""
    pts = [r.ca for c in structure.chains if c.chain_type == chain_type
           for r in c.residues if r.ca is not None]
    return np.asarray(pts).mean(axis=0) if pts else None


def _count_ca(structure: Structure, chain_type: str) -> int:
    return sum(1 for c in structure.chains if c.chain_type == chain_type
               for r in c.residues if r.ca is not None)


def _chain_ca(structure: Structure, types) -> np.ndarray:
    pts = [r.ca for c in structure.chains if c.chain_type in types
           for r in c.residues if r.ca is not None]
    return np.asarray(pts) if pts else np.empty((0, 3))


def _groove_frame(structure: Structure) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Orthonormal groove frame ``(u, w, n)`` built from the peptide + MHC→TCR direction.

    Unlike the whole-complex PCA basis, this is immune to the TCR's own Vα–Vβ spread:

    * ``u`` = groove long axis — principal axis of the peptide Cα cloud, signed N→C;
    * ``n`` = groove normal (MHC→TCR), the peptide→TCR direction orthogonalised against ``u``;
    * ``w`` = ``n × u`` — the in-plane groove-width axis (right-handed ``(u, w, n)``).
    """
    pep = _chain_ca(structure, ("PEPTIDE",))
    if len(pep) < 2:
        raise ValueError("peptide chain with ≥2 Cα required for the groove frame")
    _, _, vt = np.linalg.svd(pep - pep.mean(axis=0), full_matrices=False)
    u = vt[0]
    if np.dot(pep[-1] - pep[0], u) < 0:  # sign N→C
        u = -u
    tcr = _chain_ca(structure, _TCR_TYPES)
    if len(tcr) == 0:
        raise ValueError("no TCR Cα to orient the groove normal")
    prov = tcr.mean(axis=0) - pep.mean(axis=0)  # peptide → TCR (MHC below, TCR on top)
    n = prov - np.dot(prov, u) * u
    norm = np.linalg.norm(n)
    if norm < 1e-9:
        raise ValueError("TCR centroid is collinear with the groove axis; normal undefined")
    n = n / norm
    w = np.cross(n, u)
    return u, w, n


def crossing_incident_from_vector(v_canon: np.ndarray) -> tuple[float, float, float]:
    """``(crossing, crossing_signed, incident)`` degrees from a Vα→Vβ vector in canonical axes.

    ``v_canon`` is ``[vx, vy, vz]`` along canonical x (groove width), y (groove long axis,
    peptide N→C) and z (MHC→TCR normal). The crossing angle is measured in the groove plane
    (xy) from the long axis; the incident angle is the elevation out of that plane.
    """
    vx, vy, vz = float(v_canon[0]), float(v_canon[1]), float(v_canon[2])
    in_plane = float(np.hypot(vx, vy))
    if in_plane < 1e-9:
        raise ValueError("Vα→Vβ axis is normal to the groove plane; crossing angle undefined")
    crossing = float(np.degrees(np.arccos(np.clip(vy / in_plane, -1.0, 1.0))))  # [0, 180]
    crossing_signed = float(np.degrees(np.arctan2(vx, vy)))                     # [-180, 180]
    incident = float(np.degrees(np.arctan2(vz, in_plane)))                      # [-90, 90]
    return crossing, crossing_signed, incident


def docking_angles(structure: Structure) -> DockingAngles:
    """Crossing + incident angle of a chain-typed TCR-pMHC complex.

    The structure must already be chain-typed (``classify_chains``) and MHC-annotated
    (``annotate_mhc``) so the canonical frame can be fit. The frame is taken from the query's
    own Cα cloud (PCA), so the result needs no native database.

    Args:
        structure: a chain-typed, MHC-annotated TCR-pMHC structure.

    Returns:
        A :class:`DockingAngles` with the crossing and incident angles.

    Raises:
        ValueError: if a receptor chain pair (TRA/TRB or TRG/TRD) is missing, or the canonical
            frame is degenerate.
    """
    u, w, n = _groove_frame(structure)

    for cell_type, (va_t, vb_t) in (("ab", _AB_PAIR), ("gd", _GD_PAIR)):
        a, b = _domain_centroid(structure, va_t), _domain_centroid(structure, vb_t)
        if a is not None and b is not None:
            n_va, n_vb = _count_ca(structure, va_t), _count_ca(structure, vb_t)
            break
    else:
        raise ValueError("no complete receptor chain pair (TRA/TRB or TRG/TRD) for docking angle")

    # Vα→Vβ vector expressed in the groove frame: [width (w), long axis (u), normal (n)].
    d = b - a
    v = np.array([np.dot(d, w), np.dot(d, u), np.dot(d, n)])
    crossing, crossing_signed, incident = crossing_incident_from_vector(v)
    return DockingAngles(crossing, crossing_signed, incident, cell_type, n_va, n_vb)
