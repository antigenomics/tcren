"""Canonical TCR-pMHC frame by PCA: z ≈ PC1 (MHC→TCR), y ≈ PC2 (peptide), x ≈ PC3.

Every structure is first superposed onto a per-class native reference by its MHC groove Cα
(:func:`tcren.native.align.align_to_native`); a fixed per-class rotation ``R_canon`` then maps
that reference frame into the canonical axes. ``R_canon`` is obtained by centring the reference
complex's Cα cloud at its centre of mass and taking its principal axes (PCA):

* ``z`` = PC1 (largest variance, the MHC→TCR long axis), signed ``+z`` toward the TCR so the
  MHC sits at ``−z``;
* ``y`` = PC2 (the groove/peptide axis), signed ``+y`` toward the peptide C-terminus;
* ``x`` = PC3 (the thin axis), signed for a right-handed frame.

``R_canon`` + the variance fractions are cached in the bundled ``tcren/data/canonical_frame.json``
so orientation is reproducible and inspectable. When no native database is available the same
PCA axes are fit directly from the query (the PCA fallback).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources
from typing import Literal

import numpy as np

from ..native.align import DEFAULT_REFERENCE, align_to_native, _reference_structure
from ..native.database import NativeDatabase
from ..structure.model import Structure

_GROOVE_FLOOR = "GROOVE_FLOOR"
_TCR_TYPES = ("TRA", "TRB", "TRD", "TRG")
_VJ_TYPES = ("TRA", "TRG")
_VDJ_TYPES = ("TRB", "TRD")


@dataclass(slots=True)
class CanonResult:
    """Composed rigid transform that maps a structure into the canonical frame."""

    rotation: np.ndarray  # (3, 3); canonical = coord @ rotation + translation
    translation: np.ndarray  # (3,)
    rmsd: float
    n_anchor_atoms: int
    reference_id: str | None
    frame: Literal["native", "pca"]
    reversed_dock: bool | None = None
    chain_map: dict[str, str] = field(default_factory=dict)


def _chain_ca(structure: Structure, types) -> np.ndarray:
    pts = [r.ca for c in structure.chains if c.chain_type in types
           for r in c.residues if r.ca is not None]
    return np.asarray(pts) if pts else np.empty((0, 3))


def _peptide_termini(structure: Structure) -> tuple[np.ndarray, np.ndarray] | None:
    """(N-terminal Cα, C-terminal Cα) of the peptide chain, or ``None`` if unavailable."""
    for c in structure.chains:
        if c.chain_type == "PEPTIDE":
            cas = [r.ca for r in c.residues if r.ca is not None]
            if len(cas) >= 2:
                return cas[0], cas[-1]
    return None


_FRAME_TYPES = ("MHCa", "MHCb", "PEPTIDE", "TRA", "TRB", "TRD", "TRG")


def _canonical_basis(structure: Structure) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Right-handed canonical basis (rows x,y,z), centre of mass, and variance fractions.

    The complex Cα cloud is centred at its centre of mass and decomposed by PCA: ``z`` = PC1
    (largest variance, the MHC→TCR long axis), ``y`` = PC2 (the groove/peptide axis), ``x`` =
    PC3 (the thin axis). Signs are fixed from biology — ``+z`` toward the TCR (MHC at ``−z``),
    ``+y`` toward the peptide C-terminus — and ``x`` is signed for a right-handed frame.
    Raises ``ValueError`` on degenerate geometry.
    """
    ca = np.asarray([r.ca for c in structure.chains if c.chain_type in _FRAME_TYPES
                     for r in c.residues if r.ca is not None])
    if len(ca) < 3:
        raise ValueError("too few Cα atoms to define the canonical frame")
    com = ca.mean(axis=0)
    _, s, vt = np.linalg.svd(ca - com, full_matrices=False)
    z, y, x = vt[0].copy(), vt[1].copy(), vt[2].copy()  # PC1, PC2, PC3
    var = (s ** 2) / float(np.sum(s ** 2))
    tcr = _chain_ca(structure, _TCR_TYPES)
    if len(tcr) == 0:
        raise ValueError("no TCR Cα to orient the long axis")
    if np.dot(tcr.mean(axis=0) - com, z) < 0:           # +z toward the TCR
        z = -z
    term = _peptide_termini(structure)
    if term is not None and abs(np.dot(term[1] - term[0], y)) > 1e-6:
        if np.dot(term[1] - term[0], y) < 0:            # +y toward the peptide C-terminus
            y = -y
    basis = np.array([x, y, z])
    if np.linalg.det(basis) < 0:                        # enforce right-handedness (x = ±PC3)
        x = -x
        basis = np.array([x, y, z])
    return basis, com, var[:3]


def _transform_from_basis(basis: np.ndarray, centroid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation/translation mapping a coord into the canonical axes: ``(p - c) @ basis.T``."""
    rot = basis.T
    return rot, -centroid @ rot


@lru_cache(maxsize=4)
def _reference_canon(reference_id: str, root: str | None) -> tuple[np.ndarray, np.ndarray]:
    """``R_canon`` (rotation, translation) for a class reference, from the bundled artifact
    if present, else computed from the reference structure's own groove/peptide geometry."""
    bundled = _load_canonical_frame()
    for entry in bundled.values():
        if entry.get("reference_id") == reference_id:
            return np.asarray(entry["rotation"]), np.asarray(entry["translation"])
    ref = _reference_structure(reference_id, root)
    basis, centroid, _ = _canonical_basis(ref)
    return _transform_from_basis(basis, centroid)


def _load_canonical_frame() -> dict:
    try:
        return json.loads(resources.files("tcren.data").joinpath("canonical_frame.json").read_text())
    except (FileNotFoundError, ModuleNotFoundError):
        return {}


def canonical_frame(
    structure: Structure,
    db: NativeDatabase | None = None,
    reference_id: str | None = None,
    force_pca: bool = False,
) -> CanonResult:
    """Compose the MHC superposition with the per-class ``R_canon`` (native), or fit the
    canonical axes directly from the query (PCA fallback when no DB / too few anchors)."""
    if not force_pca:
        try:
            db = db or NativeDatabase()
            align = align_to_native(structure, db=db, reference_id=reference_id)
            r_canon, t_canon = _reference_canon(align.reference_id, str(db.root))
            rot = align.rotation @ r_canon
            tran = align.translation @ r_canon + t_canon
            return CanonResult(rot, tran, align.rmsd, align.n_anchor_atoms,
                               align.reference_id, "native")
        except Exception:
            if reference_id is not None:
                raise
    basis, centroid, _ = _canonical_basis(structure)
    rot, tran = _transform_from_basis(basis, centroid)
    return CanonResult(rot, tran, float("nan"), 0, None, "pca")


def build_canonical_frame(db: NativeDatabase | None = None) -> dict:
    """(Re)compute ``R_canon`` for each class reference and return the artifact dict.

    Writes nothing; the caller serialises to ``tcren/data/canonical_frame.json``.
    """
    db = db or NativeDatabase()
    out: dict = {}
    for mhc_class, reference_id in DEFAULT_REFERENCE.items():
        try:
            ref = _reference_structure(reference_id, str(db.root))
            basis, centroid, var = _canonical_basis(ref)
            rot, tran = _transform_from_basis(basis, centroid)
            out[mhc_class] = {
                "reference_id": reference_id,
                "rotation": rot.tolist(),
                "translation": tran.tolist(),
                "variance_explained": {"PC1_z": float(var[0]), "PC2_y": float(var[1]),
                                       "PC3_x": float(var[2])},
            }
        except Exception as exc:  # noqa: BLE001
            out[mhc_class] = {"reference_id": reference_id, "error": str(exc)[:120]}
    return out
