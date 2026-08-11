"""Ring-stacking geometry between residue side chains.

A contact potential scores a pair of residues by their identities and nothing else, so it
treats two rings lying face to face at 3.5 Å exactly like the same two residues brushing
past edge-on. Stacking is a directional interaction and that difference is the whole of it.
This module measures it from coordinates instead: how far apart two ring centroids are, how
nearly parallel the ring planes are, and how far the rings are displaced sideways.

Proline is included among the rings although it is not aromatic. Its pyrrolidine ring packs
face-on against aromatic side chains through CH--pi contacts, and leaving it out would miss
exactly the interaction this module exists to measure.

The readout is deliberately geometric and carries no energy. Nothing here says a stack is
worth some number of kT; it says the rings are or are not arranged the way a stack is.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl

from .structure.model import Chain, Residue, Structure

#: Ring atoms per residue. Six-membered rings for the aromatics, the imidazole for His, the
#: pyrrolidine for Pro. Trp is represented by its six-membered ring, which is the face that
#: stacks.
RING_ATOMS: dict[str, tuple[str, ...]] = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "PRO": ("N", "CA", "CB", "CG", "CD"),
}


@dataclass(frozen=True)
class Ring:
    """One side-chain ring: where it is and which way it faces."""

    chain_id: str
    seq_index: int
    resname: str
    centroid: np.ndarray
    normal: np.ndarray


@dataclass(frozen=True)
class RingPair:
    """Geometry of two rings, in the terms that distinguish a stack from a brush-past.

    Attributes:
        centroid_distance: Between ring centres (Å).
        interplanar_angle: Between the ring planes (degrees, 0--90). Near 0 is face-to-face,
            near 90 is edge-to-face.
        vertical: Centroid separation along the first ring's normal (Å) --- the gap between
            the planes.
        lateral: Centroid separation within that plane (Å) --- how far the rings slide past
            each other. A parallel-displaced stack has a small vertical and a lateral of a
            couple of Å; a perfectly stacked pair has both small.
    """

    a: Ring
    b: Ring
    centroid_distance: float
    interplanar_angle: float
    vertical: float
    lateral: float


def ring_of(residue: Residue, chain_id: str) -> Ring | None:
    """The ring of one residue, or ``None`` if it has none or is missing ring atoms."""
    names = RING_ATOMS.get(residue.resname)
    if names is None:
        return None
    coords = np.array([a.coord for a in residue.atoms if a.name in names], dtype=float)
    if len(coords) < len(names):        # an incompletely modelled side chain has no plane
        return None
    centroid = coords.mean(axis=0)
    # The plane normal is the least-varying direction of the ring atoms.
    _, _, vt = np.linalg.svd(coords - centroid)
    return Ring(chain_id, residue.seq_index, residue.resname, centroid, vt[2])


def rings(source: Structure | Chain) -> list[Ring]:
    """Every ring in a structure or a single chain, in residue order."""
    chains = source.chains if isinstance(source, Structure) else [source]
    out = []
    for chain in chains:
        for residue in chain.residues:
            found = ring_of(residue, chain.chain_id)
            if found is not None:
                out.append(found)
    return out


def ring_pair(a: Ring, b: Ring) -> RingPair:
    """Geometry of one ring pair."""
    separation = b.centroid - a.centroid
    distance = float(np.linalg.norm(separation))
    angle = float(np.degrees(np.arccos(np.clip(abs(float(a.normal @ b.normal)), 0.0, 1.0))))
    vertical = abs(float(separation @ a.normal))
    lateral = float(np.sqrt(max(distance**2 - vertical**2, 0.0)))
    return RingPair(a, b, distance, angle, vertical, lateral)


def ring_stacking(
    source: Structure | Chain,
    cutoff: float = 7.5,
    min_seq_sep: int = 1,
) -> pl.DataFrame:
    """All ring pairs whose centroids fall within ``cutoff``.

    Args:
        source: A parsed structure, or one chain of it.
        cutoff: Maximum centroid separation (Å). The default is generous: a stack sits near
            5 Å, and pairs beyond that are worth seeing in order to say they are not stacks.
        min_seq_sep: Minimum ``|i - j|`` for two rings on the same chain, so that sequence
            neighbours held together by the backbone are not reported as stacks.

    Returns:
        One row per pair, sorted by centroid distance: chain and residue identifiers on both
        sides, then ``centroid_distance``, ``interplanar_angle``, ``vertical``, ``lateral``.
    """
    found = rings(source)
    rows = []
    for i, a in enumerate(found):
        for b in found[i + 1:]:
            if a.chain_id == b.chain_id and abs(b.seq_index - a.seq_index) < min_seq_sep:
                continue
            pair = ring_pair(a, b)
            if pair.centroid_distance > cutoff:
                continue
            rows.append({
                "chain.id.from": a.chain_id, "residue.index.from": a.seq_index,
                "resname.from": a.resname,
                "chain.id.to": b.chain_id, "residue.index.to": b.seq_index,
                "resname.to": b.resname,
                "centroid_distance": pair.centroid_distance,
                "interplanar_angle": pair.interplanar_angle,
                "vertical": pair.vertical, "lateral": pair.lateral,
            })
    schema = {
        "chain.id.from": pl.Utf8, "residue.index.from": pl.Int64, "resname.from": pl.Utf8,
        "chain.id.to": pl.Utf8, "residue.index.to": pl.Int64, "resname.to": pl.Utf8,
        "centroid_distance": pl.Float64, "interplanar_angle": pl.Float64,
        "vertical": pl.Float64, "lateral": pl.Float64,
    }
    if not rows:
        return pl.DataFrame(schema=schema)
    return pl.DataFrame(rows, schema=schema).sort("centroid_distance")
