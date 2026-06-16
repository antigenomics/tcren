"""Project an oriented TCR-peptide-MHC interface onto a single 2D plane.

The "optimal plane" is the MHC groove plane. Two routes produce it:

* **native** (default) — orient the structure onto a canonical reference via
  :func:`tcren.native.align.align_to_native`; the canonical groove plane is xy and its
  normal is z, so projecting is just dropping z. The transform is applied only to the
  extracted Cα coordinates (``apply_transform`` would clear region annotations).
* **pca** (fallback) — fit a plane to the groove-floor Cα by SVD when no native database
  is available. The normal is the lowest-variance axis; z is sign-oriented toward the
  peptide.

The residues projected are the CDR1-3 loops (TRA/TRB), the peptide, and the MHC groove
helices/floor — the residues that line the interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..structure.model import Structure

# Region types projected onto the map, by which chain they sit on.
_TCR_REGIONS = {"CDR1", "CDR2", "CDR3"}
_MHC_REGIONS = {"HELIX_A1", "HELIX_A2", "HELIX_B1", "GROOVE_FLOOR"}
_TCR_TYPES = ("TRA", "TRB", "TRD", "TRG", "IGH", "IGK", "IGL")


@dataclass(slots=True)
class ProjectionResult:
    """2D projection of selected interface residues onto the groove plane."""

    keys: list[tuple[str, int]]  # (chain_id, seq_index) in row order
    coords3d: np.ndarray  # (N, 3) coordinates in the groove frame
    frame: Literal["native", "pca"]
    reference_id: str | None = None
    rmsd: float | None = None

    @property
    def uv(self) -> np.ndarray:
        """In-plane (u, v) coordinates."""
        return self.coords3d[:, :2]

    @property
    def height(self) -> np.ndarray:
        """Signed distance above the groove plane (z)."""
        return self.coords3d[:, 2]


def _interface_residues(structure: Structure):
    """Yield ``(chain, residue)`` for CDR loops, peptide, and MHC groove residues."""
    for chain in structure.chains:
        if chain.chain_type == "PEPTIDE":
            for res in chain.residues:
                yield chain, res
        elif chain.chain_type in _TCR_TYPES:
            wanted = {r.seq_index for reg in chain.regions if reg.region_type in _TCR_REGIONS
                      for r in reg.residues}
            for res in chain.residues:
                if res.seq_index in wanted:
                    yield chain, res
        elif chain.chain_type in ("MHCa", "MHCb"):
            wanted = {r.seq_index for reg in chain.regions if reg.region_type in _MHC_REGIONS
                      for r in reg.residues}
            for res in chain.residues:
                if res.seq_index in wanted:
                    yield chain, res


def _collect_ca(structure: Structure):
    """Return ``(keys, ca_coords)`` for the interface residues that have a Cα."""
    keys, coords = [], []
    for chain, res in _interface_residues(structure):
        ca = res.ca
        if ca is not None:
            keys.append((chain.chain_id, res.seq_index))
            coords.append(ca)
    return keys, (np.asarray(coords, dtype=np.float64) if coords else np.empty((0, 3)))


def _groove_floor_ca(structure: Structure) -> np.ndarray:
    pts = []
    for chain in structure.chains:
        if chain.chain_type in ("MHCa", "MHCb"):
            for reg in chain.regions:
                if reg.region_type == "GROOVE_FLOOR":
                    pts.extend(r.ca for r in reg.residues if r.ca is not None)
    return np.asarray(pts, dtype=np.float64) if pts else np.empty((0, 3))


def _peptide_centroid(structure: Structure) -> np.ndarray | None:
    pts = [r.ca for c in structure.chains if c.chain_type == "PEPTIDE"
           for r in c.residues if r.ca is not None]
    return np.mean(pts, axis=0) if pts else None


def _pca_frame(structure: Structure, keys, ca: np.ndarray) -> ProjectionResult:
    floor = _groove_floor_ca(structure)
    if len(floor) < 3:
        raise ValueError("need >=3 groove-floor Cα for the PCA projection plane")
    centroid = floor.mean(axis=0)
    _u, _s, vt = np.linalg.svd(floor - centroid, full_matrices=True)
    axis_u, axis_v, normal = vt[0], vt[1], vt[2]
    # Orient the normal toward the peptide ("looking down into the groove").
    pep = _peptide_centroid(structure)
    if pep is not None and np.dot(pep - centroid, normal) < 0:
        normal = -normal
        axis_v = -axis_v  # keep a right-handed (u, v, normal) basis
    basis = np.vstack([axis_u, axis_v, normal])  # (3, 3) rows = axes
    coords3d = (ca - centroid) @ basis.T
    return ProjectionResult(keys=keys, coords3d=coords3d, frame="pca")


def _tcr_up(structure: Structure, keys, coords3d: np.ndarray) -> np.ndarray:
    """Flip z so the TCR (CDR loops) sits above the peptide ("TCR on top, pMHC below").

    Keeps a right-handed frame by also negating the u axis. No-op if there is no peptide
    or no CDR residue to compare.
    """
    row = {k: i for i, k in enumerate(keys)}
    cdr_idx, pep_idx = [], []
    for chain in structure.chains:
        if chain.chain_type == "PEPTIDE":
            pep_idx += [row[(chain.chain_id, r.seq_index)] for r in chain.residues
                        if (chain.chain_id, r.seq_index) in row]
        elif chain.chain_type in _TCR_TYPES:
            for reg in chain.regions:
                if reg.region_type in _TCR_REGIONS:
                    cdr_idx += [row[(chain.chain_id, r.seq_index)] for r in reg.residues
                                if (chain.chain_id, r.seq_index) in row]
    if cdr_idx and pep_idx and coords3d[cdr_idx, 2].mean() < coords3d[pep_idx, 2].mean():
        coords3d = coords3d.copy()
        coords3d[:, 2] *= -1.0
        coords3d[:, 0] *= -1.0
    return coords3d


def project_structure(
    structure: Structure,
    db=None,
    reference_id: str | None = None,
    force_pca: bool = False,
) -> ProjectionResult:
    """Project the interface residues of an annotated structure onto the groove plane.

    Args:
        structure: A chain-typed, MHC-annotated structure.
        db: Native database for the canonical reference (native frame). If ``None`` the
            default database is used; falls back to PCA if it is unavailable.
        reference_id: Canonical reference complex id (defaults per MHC class).
        force_pca: Skip the native frame and fit the plane by PCA.

    Returns:
        A :class:`ProjectionResult`.
    """
    keys, ca = _collect_ca(structure)
    if len(ca) == 0:
        raise ValueError("no interface residues with Cα to project")

    if not force_pca:
        try:
            from ..native.align import align_to_native

            result = align_to_native(structure, db=db, reference_id=reference_id)
            coords3d = _tcr_up(structure, keys, ca @ result.rotation + result.translation)
            return ProjectionResult(
                keys=keys, coords3d=coords3d, frame="native",
                reference_id=result.reference_id, rmsd=result.rmsd,
            )
        except Exception:
            pass  # fall back to PCA (no DB, too few anchors, etc.)
    pca = _pca_frame(structure, keys, ca)
    pca.coords3d = _tcr_up(structure, keys, pca.coords3d)
    return pca
