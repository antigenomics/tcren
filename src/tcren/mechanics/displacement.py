"""Rigid interfacial displacement, spring response and contact-graph turnover."""
from dataclasses import dataclass, replace

import numpy as np

from ..contactmap import ContactMap
from ..contacts.observations import select_conformer
from ..structure.model import TCR_TYPES, Structure
from ..topology.graph import _contact_graph_counts, _recognition_contact_edges
from .springs import _stiffness_matrix, interface_springs


@dataclass(frozen=True, slots=True)
class DisplacementResponse:
    direction: str
    interface: str
    edges: int
    beta0: int
    beta1: int
    retained: int
    gained: int
    elastic_work: float
    resisting_force: float


def interface_displacement(
    structure: Structure, distance: float = 1.0
) -> list[DisplacementResponse]:
    """Translate the TCR along one outward and four shear directions.

    Distance is in Å. The two shear axes diagonalize the in-plane spring
    stiffness. Work and resisting force use the original contact counts as
    spring weights, in reduced units. Recompute graphs at each displacement
    to measure edge turnover independently of the original springs.
    """
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError("displacement must be finite and positive")
    pose = select_conformer(structure)
    springs = interface_springs(pose, cutoff=5.0, weight="count")
    if not len(springs) or not np.isclose(np.linalg.norm(springs.axis), 1.):
        return []
    _, _, axes = np.linalg.svd(springs.axis[None, :])
    plane = axes[1:].T
    _, rotation = np.linalg.eigh(plane.T @ _stiffness_matrix(springs) @ plane)
    first, second = (plane @ rotation).T
    directions = {
        "rest": np.zeros(3),
        "pull": -springs.axis,
        "shear+1": first,
        "shear-1": -first,
        "shear+2": second,
        "shear-2": -second,
    }
    original = _recognition_contact_edges(ContactMap.from_structure(pose))
    rows = []
    for label, direction in directions.items():
        offset = distance * direction
        if label == "rest":
            graphs = original
        else:
            chains = [
                replace(chain, residues=[
                    replace(residue, atoms=tuple(
                        replace(atom, coord=atom.coord + offset) for atom in residue.atoms
                    )) for residue in chain.residues
                ]) if chain.chain_type in TCR_TYPES else chain
                for chain in pose.chains
            ]
            graphs = _recognition_contact_edges(
                ContactMap.from_structure(replace(pose, chains=chains))
            )
        vectors = springs.a + offset - springs.b
        length = np.linalg.norm(vectors, axis=1)
        extension = length - springs.rest
        work = float(0.5 * np.dot(springs.k, extension**2))
        projection = np.divide(vectors @ direction, length, out=np.zeros_like(length), where=length > 0)
        force = float(np.dot(springs.k * extension, projection)) if np.all(length > 0) else float("nan")
        for interface in ("tcr_peptide", "tcr_mhc"):
            graph = graphs[interface]
            edges, _, beta0, beta1 = _contact_graph_counts(graph)
            rows.append(DisplacementResponse(
                label, interface, edges, beta0, beta1, len(graph & original[interface]),
                len(graph - original[interface]), work, force,
            ))
    return rows
