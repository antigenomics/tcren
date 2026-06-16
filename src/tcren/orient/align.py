"""Bring a structure into a canonical reference frame by MHC superposition.

A query complex is oriented onto a native reference by superposing the conserved MHC
groove Cα atoms (the helix/floor residues from :mod:`tcren.mhc.regions`). Because every
structure is aligned to the same reference, all oriented complexes share one frame —
the basis for overlaying structures and for 2D interface projection. Correspondence
between query and reference groove residues is established by sequence alignment, so
different alleles/numbering are handled.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from ..annotation import classify_chains
from ..mhc import annotate_mhc
from ..mhc.regions import _aligner
from ..paths import reference_structure_path
from ..structure import parse_structure
from ..structure.model import Atom, Residue, Structure

# Default canonical reference complex per MHC class.
DEFAULT_REFERENCE = {"MHCI": "1ao7", "MHCII": "1fyt"}
_GROOVE = {"HELIX_A1", "HELIX_A2", "HELIX_B1", "GROOVE_FLOOR"}


@dataclass(slots=True)
class OrientationResult:
    """Rigid transform that maps a structure onto the canonical reference frame."""

    rotation: np.ndarray  # (3, 3)
    translation: np.ndarray  # (3,)
    rmsd: float
    n_anchor_atoms: int
    reference_id: str


def _groove_ca_by_seq(structure: Structure) -> dict[str, list[tuple[int, np.ndarray]]]:
    """Map each MHC role -> [(position-in-chain-sequence, Cα coord)] over groove residues.

    Position is the residue's index within its chain (for cross-structure alignment).
    """
    out: dict[str, list[tuple[int, np.ndarray]]] = {}
    for chain in structure.chains:
        if chain.chain_type not in ("MHCa", "MHCb"):
            continue
        index_of = {r.seq_index: i for i, r in enumerate(chain.residues)}
        anchors = []
        for region in chain.regions:
            if region.region_type not in _GROOVE:
                continue
            for res in region.residues:
                if res.ca is not None:
                    anchors.append((index_of[res.seq_index], res.ca, chain.sequence()))
        if anchors:
            out[chain.chain_type] = (chain.sequence(), [(p, c) for p, c, _ in anchors])
    return out


def _matched_anchors(mobile: Structure, reference: Structure):
    """Return matched (mobile, reference) groove Cα arrays across shared MHC roles."""
    mob = _groove_ca_by_seq(mobile)
    ref = _groove_ca_by_seq(reference)
    mob_pts, ref_pts = [], []
    for role in mob.keys() & ref.keys():
        mob_seq, mob_anchors = mob[role]
        ref_seq, ref_anchors = ref[role]
        mob_ca = dict(mob_anchors)
        ref_ca = dict(ref_anchors)
        # Align the two chain sequences; keep ungapped columns present in both anchor sets.
        alignment = _aligner().align(mob_seq, ref_seq)[0]
        for (qs, qe), (ts, te) in zip(*alignment.aligned):
            for off in range(qe - qs):
                qp, tp = qs + off, ts + off
                if qp in mob_ca and tp in ref_ca:
                    mob_pts.append(mob_ca[qp])
                    ref_pts.append(ref_ca[tp])
    return np.asarray(mob_pts), np.asarray(ref_pts)


@lru_cache(maxsize=4)
def _reference_structure(reference_id: str) -> Structure:
    """Load + annotate a canonical reference complex from the Native2026 dataset."""
    s = parse_structure(reference_structure_path(reference_id), pdb_id=reference_id)
    classify_chains(s, organism="human")
    annotate_mhc(s)
    return s


def align_to_native(
    structure: Structure,
    reference_id: str | None = None,
) -> OrientationResult:
    """Compute the transform orienting ``structure`` onto a native reference by MHC.

    ``structure`` must already be chain-typed and MHC-annotated (see
    :func:`tcren.mhc.annotate_mhc`). The reference (default a canonical complex for the
    structure's MHC class) is loaded from the ``Native2026`` dataset (``tcren.paths``).
    """
    from Bio.SVDSuperimposer import SVDSuperimposer

    mhc_class = next(
        (c.chain_supertype for c in structure.chains if c.chain_type in ("MHCa", "MHCb")),
        "MHCI",
    )
    reference_id = reference_id or DEFAULT_REFERENCE.get(mhc_class, "1ao7")
    reference = _reference_structure(reference_id)

    mob_pts, ref_pts = _matched_anchors(structure, reference)
    if len(mob_pts) < 3:
        raise ValueError(
            f"too few matched groove Cα anchors ({len(mob_pts)}) to orient {structure.pdb_id}"
        )
    sup = SVDSuperimposer()
    sup.set(ref_pts, mob_pts)  # reference is fixed; map mobile onto it
    sup.run()
    rot, tran = sup.get_rotran()
    return OrientationResult(
        rotation=rot,
        translation=tran,
        rmsd=float(sup.get_rms()),
        n_anchor_atoms=len(mob_pts),
        reference_id=reference_id,
    )


def apply_transform(structure: Structure, result: OrientationResult) -> Structure:
    """Return a copy of ``structure`` with the orientation transform applied to all atoms."""
    rot, tran = result.rotation, result.translation
    new_chains = []
    for chain in structure.chains:
        new_residues = [
            Residue(
                seq_index=r.seq_index,
                pdb_index=r.pdb_index,
                insertion_code=r.insertion_code,
                aa=r.aa,
                resname=r.resname,
                atoms=tuple(
                    Atom(a.name, a.element, np.dot(a.coord, rot) + tran) for a in r.atoms
                ),
            )
            for r in chain.residues
        ]
        new_chain = copy.copy(chain)
        new_chain.residues = new_residues
        new_chain.regions = []  # region residue references are stale after copy
        new_chains.append(new_chain)
    return Structure(pdb_id=structure.pdb_id, chains=new_chains, complex_species=structure.complex_species)
