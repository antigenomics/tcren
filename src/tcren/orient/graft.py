"""Build chimeric TCR:pMHC complexes by grafting one complex's TCR onto another's pMHC.

Given a **host** and a **donor** TCR:pMHC complex, :func:`substitute_tcr` produces a new complex that
keeps the **host peptide + MHC** and the **donor TCR**, rigidly positioned by one of two anchors:

* ``by="mhc"`` — superpose the donor MHC-groove Cα onto the host MHC-groove Cα, then carry the donor
  TCR into the host frame. The donor TCR keeps its **native docking geometry** relative to MHC (the
  pose it adopts on its own pMHC, transferred onto the host groove).
* ``by="tcr"`` — superpose the donor TCR Cα onto the host TCR Cα, then drop the host TCR. The donor
  TCR inherits the **host's docking pose** (it lands where the host TCR sat).

Both yield *host pMHC + donor TCR*; only the superposition anchor differs. Residue correspondence
between the two MHCs (or the two TCRs) is by sequence alignment, so different alleles / V-genes are
handled. The MHC path needs both inputs MHC-annotated (:func:`tcren.mhc.annotate_mhc`); both paths
need chain typing (:func:`tcren.annotation.classify_chains`).
"""

from __future__ import annotations

import copy

import numpy as np

from ..mhc.regions import _aligner
from ..structure.model import PEPTIDE_TYPE, TCR_TYPES, Structure
from .align import OrientationResult, _matched_anchors, apply_transform

# Fresh chain-id pool for relabelling grafted chains that would collide with the host.
_ID_POOL = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def _tcr_ca_matched(mobile: Structure, reference: Structure):
    """Matched (mobile, reference) TCR Cα arrays over shared TCR chain types, by sequence alignment."""
    def by_type(s: Structure) -> dict[str, tuple[str, dict[int, np.ndarray]]]:
        out: dict[str, tuple[str, dict[int, np.ndarray]]] = {}
        for c in s.chains:
            if c.chain_type in TCR_TYPES:
                ca = {i: r.ca for i, r in enumerate(c.residues) if r.ca is not None}
                out[c.chain_type] = (c.sequence(), ca)
        return out

    mob, ref = by_type(mobile), by_type(reference)
    mob_pts, ref_pts = [], []
    for role in mob.keys() & ref.keys():
        mseq, mca = mob[role]
        rseq, rca = ref[role]
        alignment = _aligner().align(mseq, rseq)[0]
        for (qs, qe), (ts, _te) in zip(*alignment.aligned):
            for off in range(qe - qs):
                qp, tp = qs + off, ts + off
                if qp in mca and tp in rca:
                    mob_pts.append(mca[qp])
                    ref_pts.append(rca[tp])
    return np.asarray(mob_pts), np.asarray(ref_pts)


def _superpose(mob_pts: np.ndarray, ref_pts: np.ndarray, what: str) -> OrientationResult:
    """Kabsch transform mapping ``mob_pts`` onto ``ref_pts`` (needs ≥ 3 matched anchors)."""
    from ._transform import kabsch

    if len(mob_pts) < 3:
        raise ValueError(f"too few matched {what} Cα anchors ({len(mob_pts)}) to superpose")
    rot, tran, rmsd = kabsch(mob_pts, ref_pts)
    return OrientationResult(rotation=rot, translation=tran, rmsd=rmsd,
                             n_anchor_atoms=len(mob_pts), reference_id="")


def _relabel(chains, used: set[str]):
    """Return copies of ``chains`` with ids that don't collide with ``used`` (updated in place)."""
    out = []
    for ch in chains:
        new = copy.copy(ch)
        if new.chain_id in used:
            new.chain_id = next(c for c in _ID_POOL if c not in used)
        used.add(new.chain_id)
        out.append(new)
    return out


def substitute_tcr(host: Structure, donor: Structure, *, by: str = "mhc") -> Structure:
    """Graft the ``donor`` TCR onto the ``host`` pMHC → a chimeric TCR:pMHC :class:`Structure`.

    The result keeps the host peptide + MHC chains and the donor TCR chains (relabelled to avoid id
    collisions), with the donor TCR rigidly placed by the ``by`` anchor:

    * ``"mhc"`` — superpose the donor MHC groove onto the host MHC groove (both inputs must be
      :func:`tcren.mhc.annotate_mhc`-annotated); the donor TCR keeps its native docking geometry.
    * ``"tcr"`` — superpose the donor TCR onto the host TCR; the donor TCR inherits the host pose.

    Both inputs must be chain-typed (:func:`tcren.annotation.classify_chains`). Raises ``ValueError``
    if ``by`` is invalid, the host lacks a peptide or MHC chain, the donor lacks a TCR chain, or too
    few matched Cα anchors are found to superpose.
    """
    if by not in ("mhc", "tcr"):
        raise ValueError(f"by must be 'mhc' or 'tcr', got {by!r}")

    # The host pMHC is everything that is not the TCR (peptide + MHC/β2m). Defining it by exclusion
    # keeps ``by="tcr"`` working without MHC annotation (generic-typed MHC chains still count).
    host_pmhc = [c for c in host.chains if c.chain_type not in TCR_TYPES]
    donor_tcr = donor.by_type(*TCR_TYPES)
    if not any(c.chain_type == PEPTIDE_TYPE for c in host_pmhc):
        raise ValueError(f"host {host.pdb_id!r} has no peptide chain")
    if not any(c.chain_type != PEPTIDE_TYPE for c in host_pmhc):
        raise ValueError(f"host {host.pdb_id!r} has no MHC chain")
    if not donor_tcr:
        raise ValueError(f"donor {donor.pdb_id!r} has no TCR chain")

    if by == "mhc":
        mob, ref = _matched_anchors(donor, host)  # MHC groove Cα (needs annotate_mhc on both)
        result = _superpose(mob, ref, "MHC groove")
    else:
        mob, ref = _tcr_ca_matched(donor, host)  # TCR Cα
        result = _superpose(mob, ref, "TCR")

    moved_tcr = apply_transform(donor, result).by_type(*TCR_TYPES)  # donor TCR in the host frame

    used: set[str] = set()
    host_chains = _relabel(host_pmhc, used)  # copies; ids preserved (no collision yet)
    tcr_chains = _relabel(moved_tcr, used)   # copies; collisions relabelled from _ID_POOL
    return Structure(
        pdb_id=f"{host.pdb_id}_{donor.pdb_id}_{by}",
        chains=host_chains + tcr_chains,
        complex_species=host.complex_species,
        cell_type=donor.cell_type,
    )
