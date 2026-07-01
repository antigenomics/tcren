"""Peptide RMSD between two poses of the same complex.

The benchmark question is "how close is a re-modelled peptide to its native crystal pose?" We answer it
the way the rest of the codebase measures geometry: superpose the two complexes on the conserved MHC
groove Cα (:func:`tcren.orient.align._matched_anchors` + Biopython's ``SVDSuperimposer``), then compute
the peptide backbone RMSD in that common MHC frame. Superposing on the MHC — not on the peptide —
means the peptide RMSD reflects how well the pose sits in the groove, which is the quantity that
matters for downstream TCR-contact scoring.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..structure.model import PEPTIDE_TYPE, Structure

_BACKBONE = ("N", "CA", "C", "O")


def _peptide_atoms(structure: Structure, names: tuple[str, ...]) -> dict[tuple[int, str], np.ndarray]:
    """Map (residue seq_index, atom name) -> coord for the peptide chain's selected atoms."""
    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in {structure.pdb_id!r}")
    out: dict[tuple[int, str], np.ndarray] = {}
    for res in pep.residues:
        for a in res.atoms:
            if a.name in names:
                out[(res.seq_index, a.name)] = a.coord
    return out


@dataclass(frozen=True, slots=True)
class PeptideRMSD:
    """Peptide RMSD of a model against a reference, after MHC-groove superposition."""

    backbone_rmsd: float  # over N, CA, C, O present in both
    ca_rmsd: float  # over Cα only
    anchor_ca_rmsd: float  # over Cα at the anchor positions only (nan if no anchors given)
    n_backbone: int
    n_ca: int
    groove_rmsd: float  # the MHC-groove superposition residual (frame quality)


def peptide_rmsd(
    model: Structure,
    reference: Structure,
    anchors: tuple[int, ...] = (),
) -> PeptideRMSD:
    """Peptide backbone / Cα / anchor-Cα RMSD of ``model`` vs ``reference`` in the MHC-groove frame.

    Both structures must be chain-typed and MHC-annotated. ``model`` is superposed onto ``reference``
    by their shared groove Cα; the peptide RMSDs are then computed over atoms (and residues) present
    in both. ``anchors`` are 0-based peptide residue indices for the anchor-Cα RMSD.
    """
    from Bio.SVDSuperimposer import SVDSuperimposer

    from ..orient.align import _matched_anchors

    mob_pts, ref_pts = _matched_anchors(model, reference)
    if len(mob_pts) < 3:
        raise ValueError(f"too few matched groove Cα to superpose {model.pdb_id!r}")
    sup = SVDSuperimposer()
    sup.set(ref_pts, mob_pts)
    sup.run()
    rot, tran = sup.get_rotran()
    groove_rmsd = float(sup.get_rms())

    mod_atoms = _peptide_atoms(model, _BACKBONE)
    ref_atoms = _peptide_atoms(reference, _BACKBONE)
    shared = sorted(mod_atoms.keys() & ref_atoms.keys())
    if not shared:
        raise ValueError("no shared peptide backbone atoms between model and reference")

    bb_m = np.array([np.dot(mod_atoms[k], rot) + tran for k in shared])
    bb_r = np.array([ref_atoms[k] for k in shared])
    bb_rmsd = float(np.sqrt(((bb_m - bb_r) ** 2).sum(1).mean()))

    ca_keys = [k for k in shared if k[1] == "CA"]
    ca_m = np.array([np.dot(mod_atoms[k], rot) + tran for k in ca_keys])
    ca_r = np.array([ref_atoms[k] for k in ca_keys])
    ca_rmsd = float(np.sqrt(((ca_m - ca_r) ** 2).sum(1).mean()))

    anchor_set = set(anchors)
    anc_keys = [k for k in ca_keys if k[0] in anchor_set]
    if anc_keys:
        a_m = np.array([np.dot(mod_atoms[k], rot) + tran for k in anc_keys])
        a_r = np.array([ref_atoms[k] for k in anc_keys])
        anchor_rmsd = float(np.sqrt(((a_m - a_r) ** 2).sum(1).mean()))
    else:
        anchor_rmsd = float("nan")

    return PeptideRMSD(bb_rmsd, ca_rmsd, anchor_rmsd, len(shared), len(ca_keys), groove_rmsd)
