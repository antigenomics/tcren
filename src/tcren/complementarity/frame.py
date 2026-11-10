"""Explicit right-handed frames for side (x,z) and top (x,y) contact views.

Transforms use row vectors. A frame is constructed once on a reference complex;
mobile structures inherit it through matched fixed-component landmarks.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np


def _points(value, name: str, minimum: int = 3) -> np.ndarray:
    points = np.asarray(value, dtype=float)
    if (points.ndim != 2 or points.shape[1] != 3 or len(points) < minimum
            or not np.isfinite(points).all()):
        raise ValueError(f"{name} requires >= {minimum} finite 3D landmarks")
    return points


def _unit(vector: np.ndarray, name: str) -> np.ndarray:
    length = np.linalg.norm(vector)
    if length < 1e-8:
        raise ValueError(f"degenerate {name}")
    return vector / length


@dataclass(frozen=True)
class CanonicalFrame:
    """Rigid transform and reference identity; distances remain in Angstroms."""

    rotation: np.ndarray
    translation: np.ndarray
    reference_id: str
    origin: str
    rmsd: float = 0.0
    anchor_count: int = 0

    def transform(self, points) -> np.ndarray:
        """Apply this frame without modifying the source coordinates."""
        array = np.asarray(points, dtype=float)
        if array.shape[-1:] != (3,) or not np.isfinite(array).all():
            raise ValueError("coordinates must be finite with final dimension 3")
        return array @ self.rotation + self.translation


def canonical_frame(
    peptide_ca,
    groove_ca,
    *,
    reference_id: str,
    origin: Literal["apex", "center-of-mass"] = "center-of-mass",
    peptide_atoms=None,
    peptide_masses=None,
) -> CanonicalFrame:
    """Construct the reference frame from ordered peptide and groove landmarks.

    +x follows the resolved peptide termini N→C. +z is the groove-plane normal
    toward the peptide, orthogonalized against x. +y completes the right-handed
    basis. The apex is the peptide C-alpha of greatest height over that plane
    (first in sequence on ties). Center of mass uses explicitly supplied peptide
    atom coordinates and positive masses; it never substitutes a CA centroid.
    Groove landmarks must span a plane; a peptide lying in that plane is ambiguous.
    """
    pep = _points(peptide_ca, "peptide", 2)
    groove = _points(groove_ca, "groove")
    center = groove.mean(axis=0)
    _, singular, axes = np.linalg.svd(groove - center, full_matrices=False)
    if singular[1] < 1e-8:
        raise ValueError("groove landmarks are collinear")
    normal = axes[-1]
    toward = np.dot(pep.mean(axis=0) - center, normal)
    if abs(toward) < 1e-8:
        raise ValueError("peptide does not define the upper groove face")
    normal *= np.sign(toward)
    x = _unit(pep[-1] - pep[0], "peptide termini")
    z = _unit(normal - np.dot(normal, x) * x, "groove normal")
    y = _unit(np.cross(z, x), "transverse axis")
    rotation = np.column_stack((x, y, z))
    if origin == "apex":
        position = pep[np.argmax((pep - center) @ z)]
    elif origin == "center-of-mass":
        if peptide_atoms is None or peptide_masses is None:
            raise ValueError("center-of-mass requires peptide atoms and masses")
        atoms = _points(peptide_atoms, "peptide atoms", 1)
        masses = np.asarray(peptide_masses, dtype=float)
        if (masses.shape != (len(atoms),) or not np.isfinite(masses).all()
                or np.any(masses <= 0)):
            raise ValueError("one finite positive mass is required per peptide atom")
        position = np.average(atoms, axis=0, weights=masses)
    else:
        raise ValueError("origin must be apex or center-of-mass")
    return CanonicalFrame(rotation, -position @ rotation, reference_id, origin)


def fit_frame(mobile_anchors, reference_anchors, reference: CanonicalFrame) -> CanonicalFrame:
    """Fit matched fixed-component landmarks and inherit one reference frame.

    Callers own sequence correspondence: use MHC groove anchors for TCR ranking,
    TCR framework anchors for epitope ranking. This function never recenters a
    mobile peptide, changes a frame sign or permits reflection.
    """
    mobile = _points(mobile_anchors, "mobile anchors")
    fixed = _points(reference_anchors, "reference anchors")
    if mobile.shape != fixed.shape:
        raise ValueError("anchor arrays must have identical shapes")
    mc, fc = mobile.mean(axis=0), fixed.mean(axis=0)
    a, b = mobile - mc, fixed - fc
    if min(np.linalg.matrix_rank(a), np.linalg.matrix_rank(b)) < 2:
        raise ValueError("matched anchors are collinear")
    u, _, vt = np.linalg.svd(a.T @ b)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    translation = fc - mc @ rotation
    rmsd = float(np.sqrt(np.mean(np.sum((mobile @ rotation + translation - fixed)**2, axis=1))))
    return CanonicalFrame(
        rotation @ reference.rotation,
        translation @ reference.rotation + reference.translation,
        reference.reference_id, reference.origin, rmsd, len(mobile),
    )
