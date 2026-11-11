"""Select structural landmarks and match one fixed component across a group."""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, replace

from .frame import CanonicalFrame, canonical_frame, fit_frame
from .topology import conformer_rank

_MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06,
         "P": 30.974, "SE": 78.971}
_GROOVE = {"GROOVE_FLOOR", "HELIX_A1", "HELIX_A2", "HELIX_B1"}
_FRAMEWORK = {"FR1", "FR2", "FR3", "FR4"}


def _mass_points(residues):
    points, masses = [], []
    for res in residues:
        seen = set()
        for atom in sorted(res.atoms, key=lambda a: conformer_rank(a.altloc)):
            element = atom.element.upper()
            if element in {"H", "D"} or atom.name in seen:
                continue
            seen.add(atom.name)
            if element not in _MASS:
                raise ValueError(f"unsupported mass element {element}")
            points.append(atom.coord)
            masses.append(_MASS[element])
    if not points:
        raise ValueError("centering requires heavy atoms")
    return points, masses


def reference_frame(structure, *, origin="center-of-mass", task="tcr"):
    """Build a frame from one peptide, annotated groove floor and heavy-atom mass.

    Center of mass explicitly excludes hydrogen/deuterium and unrecognized
    elements fail instead of receiving a guessed mass.
    """
    peptides = structure.by_type("PEPTIDE")
    if len(peptides) != 1:
        raise ValueError("reference requires exactly one peptide chain")
    peptide = peptides[0]
    ca = [r.ca for r in peptide.residues if r.ca is not None]
    groove = [r.ca for c in structure.by_type("MHCa", "MHCb") for reg in c.regions
              if reg.region_type == "GROOVE_FLOOR" for r in reg.residues if r.ca is not None]
    if task not in {"tcr", "epitope"} or origin not in {"apex", "center-of-mass"}:
        raise ValueError("task must be tcr/epitope and origin apex/center-of-mass")
    initial_origin = origin if task == "tcr" else "apex"
    atoms, masses = _mass_points(peptide.residues) if initial_origin == "center-of-mass" else (None, None)
    frame = canonical_frame(ca, groove, reference_id=structure.pdb_id, origin=initial_origin,
                            peptide_atoms=atoms, peptide_masses=masses)
    if task == "tcr":
        return frame
    tcr = _fixed_chains(structure, task)
    selected = [res for c in tcr.values() for reg in c.regions
                if reg.region_type in _FRAMEWORK | {"CDR1", "CDR2", "CDR3"} for res in reg.residues]
    unique = {id(res): res for res in selected}.values()
    if origin == "center-of-mass":
        points, weights = _mass_points(unique)
        position = np.average(points, axis=0, weights=weights)
    else:
        tips = [r.ca for c in tcr.values() for reg in c.regions
                if reg.region_type in {"CDR1", "CDR2", "CDR3"} for r in reg.residues if r.ca is not None]
        if not tips:
            raise ValueError("TCR apex requires CDR C-alpha landmarks")
        position = tips[np.argmin(frame.transform(tips)[:, 2])]
    return CanonicalFrame(frame.rotation, -position @ frame.rotation, frame.reference_id,
                          origin, origin_component="tcr-variable-domain")


def _fixed_chains(structure, task):
    types = ("MHCa", "MHCb") if task == "tcr" else ("TRA", "TRB")
    result = {}
    for chain in structure.by_type(*types):
        if chain.chain_type in result:
            raise ValueError(f"multiple {chain.chain_type} chains: explicit construct selection required")
        result[chain.chain_type] = chain
    if not result or (task == "epitope" and set(result) != {"TRA", "TRB"}):
        raise ValueError(f"missing fixed component for {task} ranking")
    return result


@dataclass(frozen=True)
class FixedComponent:
    """Execution-owned reference landmarks, without retaining a whole complex."""

    task: str
    sequences: dict
    landmarks: dict


def fixed_component(structure, task):
    """Keep only sequence identity and required reference CA landmarks."""
    if task not in {"tcr", "epitope"}:
        raise ValueError("task must be tcr or epitope")
    chains = _fixed_chains(structure, task)
    regions = _GROOVE if task == "tcr" else _FRAMEWORK
    sequences = {role: c.sequence() for role, c in chains.items()}
    if task == "tcr":
        peptides = structure.by_type("PEPTIDE")
        if len(peptides) != 1:
            raise ValueError("receptor-ranking group requires exactly one peptide")
        sequences["PEPTIDE"] = peptides[0].sequence()
    landmarks = {}
    for role, chain in chains.items():
        allowed = {r.seq_index for reg in chain.regions if reg.region_type in regions for r in reg.residues}
        for ordinal, residue in enumerate(chain.residues):
            if residue.seq_index in allowed and residue.ca is not None:
                landmarks[role, ordinal] = np.array(residue.ca, dtype=float, copy=True)
    return FixedComponent(task, sequences, landmarks)


def ranking_frame(structure, reference, frame, *, task="tcr"):
    """Sequence-match annotated fixed anchors, fit, and inherit reference origin.

    Receptor ranking matches identical peptide and MHC chain sequences before
    alignment. Epitope ranking requires identical alpha/beta sequences. All
    included CA anchors are paired by sequence alignment, never author numbering.
    Different constructs are rejected rather than mixed into a misleading overlay.
    """
    if task not in {"tcr", "epitope"}:
        raise ValueError("task must be tcr or epitope")
    mobile = fixed_component(structure, task)
    fixed = reference if isinstance(reference, FixedComponent) else fixed_component(reference, task)
    if fixed.task != task or mobile.sequences != fixed.sequences:
        raise ValueError("fixed-component sequences differ; use a separate reference group")
    keys = sorted(mobile.landmarks.keys() & fixed.landmarks.keys())
    source = np.array([mobile.landmarks[key] for key in keys], dtype=float).reshape((-1, 3))
    target = np.array([fixed.landmarks[key] for key in keys], dtype=float).reshape((-1, 3))
    return replace(fit_frame(source, target, frame), anchor_role_positions=tuple(keys))
