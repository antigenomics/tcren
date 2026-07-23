"""Shared rigid-transform primitives for the orientation / RMSD code.

Both live here so the four call sites that used to inline the SVDSuperimposer boilerplate
(and the several that inlined ``points @ rot + tran``) share one implementation.
"""

from __future__ import annotations

import numpy as np


def kabsch(mob_pts: np.ndarray, ref_pts: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Least-squares rigid transform mapping ``mob_pts`` onto the fixed ``ref_pts``.

    Args:
        mob_pts: ``(N, 3)`` mobile points (transformed onto the reference).
        ref_pts: ``(N, 3)`` fixed reference points.

    Returns:
        ``(rotation, translation, rmsd)`` such that ``mob_pts @ rotation + translation``
        best matches ``ref_pts``.
    """
    from Bio.SVDSuperimposer import SVDSuperimposer

    sup = SVDSuperimposer()
    sup.set(ref_pts, mob_pts)  # reference is fixed; map mobile onto it
    sup.run()
    rot, tran = sup.get_rotran()
    return rot, tran, float(sup.get_rms())


def apply_rigid(points: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """Apply a rigid transform to an ``(N, 3)`` (or ``(3,)``) array: ``points @ rotation + translation``."""
    return np.asarray(points) @ rotation + translation
