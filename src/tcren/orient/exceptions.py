"""Detect reverse-docked TCR-pMHC complexes (a biological exception, flagged not flipped).

The canonical frame is fixed by peptide polarity (``+y`` = peptide C-terminus). The conserved
diagonal docking then places the VDJ chain (TRB/TRD, Vβ) on the peptide-C side (``+y``) and the
VJ chain (TRA/TRG, Vα) on the peptide-N side (``−y``) — consistent with the CDR footprint
CDR1α·CDR2α·CDR3α·CDR3β·CDR2β·CDR1β laid out N→C. A genuinely reverse-docked TCR lands with the
α/β sides mirrored. We report it; we never force-flip, because the orientation is meaningful.
"""

from __future__ import annotations

import numpy as np

from ..structure.model import Structure

_VJ_TYPES = ("TRA", "TRG")
_VDJ_TYPES = ("TRB", "TRD")


def _mean_y(structure: Structure, types, rotation, translation) -> float | None:
    pts = [r.ca for c in structure.chains if c.chain_type in types
           for r in c.residues if r.ca is not None]
    if not pts:
        return None
    transformed = np.asarray(pts) @ rotation + translation
    return float(transformed[:, 1].mean())


def detect_reverse_dock(
    structure: Structure, rotation: np.ndarray, translation: np.ndarray, margin: float = 2.0
) -> bool | None:
    """Apply the canonical transform and check the TCR α/β handedness.

    Canonical: VDJ (TRB/TRD) at ``+y`` (peptide-C side) and VJ (TRA/TRG) at ``−y``. Returns
    ``True`` when the VJ chain is on the ``+y`` side of the VDJ chain by more than ``margin`` Å
    (reverse dock), ``False`` for a canonical dock, and ``None`` when a TCR side is missing.
    """
    vj = _mean_y(structure, _VJ_TYPES, rotation, translation)
    vdj = _mean_y(structure, _VDJ_TYPES, rotation, translation)
    if vj is None or vdj is None:
        return None
    return (vj - vdj) > margin
