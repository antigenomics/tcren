"""Pose diagnostics for a chain-typed TCR–pMHC complex."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl

from ..complementarity.frame import canonical_frame
from ..contacts.geometry import _chain_ca, all_atom_contacts
from ..contacts.table import directed_views
from ..project2d.frame import _groove_floor_ca
from ..structure.model import MHC_TYPES, PEPTIDE_TYPE, TCR_TYPES, Structure


@dataclass(slots=True)
class DockingDiagnostic:
    """Per-chain residue-pair counts, groove measurements and stable reason codes."""

    peptide_contacts: dict[str, int]
    pmhc_contacts: dict[str, int]
    height_A: float | None
    lateral_offset_A: float | None
    side_angle_deg: float | None
    cutoff_A: float
    side_threshold_deg: float
    reason_codes: tuple[str, ...]


def diagnose_docking(
    structure: Structure, *, cutoff: float = 5.0, side_threshold: float = 45.0,
) -> DockingDiagnostic:
    """Report contact failures and sideways placement of a single typed complex.

    Count a residue pair once if any heavy-atom pair is within the inclusive
    cutoff (Å). pMHC counts include peptide and all MHC chains, including B2M.
    Count all receptor residues on the observed chains.

    Fit the annotated MHC floor plane and orient its normal toward the peptide.
    Measure the receptor Cα centroid relative to the peptide Cα centroid.
    The side angle is measured from the positive normal in the transverse
    plane, ignoring displacement along the groove. At the default 45°, the
    lateral offset equals the height. A negative height also flags a pose below
    the peptide plane. This is a geometric flag for review, with a configurable
    threshold. Unavailable groove geometry is reported separately from contacts.

    Requires one αβ or γδ receptor pair, one peptide and an MHC alpha chain.
    Annotate MHC groove floor regions to obtain the side angle.
    """
    if not np.isfinite(side_threshold) or not 0 < side_threshold < 90:
        raise ValueError("side_threshold must be between 0 and 90 degrees")
    receptor = structure.by_type(*TCR_TYPES)
    peptide = structure.by_type(PEPTIDE_TYPE)
    mhc = structure.by_type(*MHC_TYPES)
    if len(receptor) != 2 or {c.chain_type for c in receptor} not in (
        {"TRA", "TRB"}, {"TRG", "TRD"},
    ):
        raise ValueError("require one complete αβ or γδ receptor pair")
    if len(peptide) != 1 or len(structure.by_type("MHCa")) != 1:
        raise ValueError("require one peptide and one MHC alpha chain")
    selected = [*receptor, *peptide, *mhc]
    if len({c.chain_id for c in selected}) != len(selected):
        raise ValueError("require distinct receptor, peptide and MHC chain IDs")
    if not all(any(a.element != "H" for r in c.residues for a in r.atoms) for c in selected):
        raise ValueError("each receptor, peptide and MHC chain needs observed heavy atoms")

    ids = [c.chain_id for c in receptor]
    peptide_id = peptide[0].chain_id
    pmhc_ids = [peptide_id, *(c.chain_id for c in mhc)]
    physical = all_atom_contacts(Structure(structure.pdb_id, selected), cutoff=cutoff)
    counts = (
        directed_views(physical)
        .filter(pl.col("chain.id.from").is_in(ids), pl.col("chain.id.to").is_in(pmhc_ids))
        .group_by("chain.id.from")
        .agg((pl.col("chain.id.to") == peptide_id).sum().alias("peptide"), pl.len().alias("pmhc"))
    )
    peptide_contacts = dict.fromkeys(ids, 0)
    pmhc_contacts = dict.fromkeys(ids, 0)
    for chain_id, n_peptide, n_pmhc in counts.iter_rows():
        peptide_contacts[chain_id] = n_peptide
        pmhc_contacts[chain_id] = n_pmhc

    reasons = []
    for partner, values in (("peptide", peptide_contacts), ("pmhc", pmhc_contacts)):
        if not any(values.values()):
            reasons.append(f"no_{partner}_contact")
        elif not all(values.values()):
            reasons.append(f"chain_no_{partner}_contact")

    pep_ca = _chain_ca(structure, (PEPTIDE_TYPE,))
    tcr_ca = [_chain_ca(structure, (chain.chain_type,)) for chain in receptor]
    height = lateral = angle = None
    try:
        frame = canonical_frame(
            pep_ca, _groove_floor_ca(structure), reference_id=structure.pdb_id, origin="apex",
        )
        if any(not len(points) or not np.isfinite(points).all() for points in tcr_ca):
            raise ValueError("each receptor chain needs finite Cα coordinates")
        displacement = np.concatenate(tcr_ca).mean(axis=0) - pep_ca.mean(axis=0)
        _, width, z = displacement @ frame.rotation
        height, lateral = float(z), float(abs(width))
        if height == 0 and lateral == 0:
            raise ValueError("receptor lies on the groove axis")
        angle = float(np.degrees(np.arctan2(lateral, height)))
        if angle >= side_threshold:
            reasons.append("side_docking")
    except ValueError:
        reasons.append("groove_geometry_unavailable")

    return DockingDiagnostic(peptide_contacts, pmhc_contacts, height, lateral, angle,
                             cutoff, side_threshold, tuple(reasons))
