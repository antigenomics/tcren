"""Class-II core correspondence after MHC-only superposition, without sequence scoring."""
from __future__ import annotations

import numpy as np
import polars as pl

from ..docking._transform import kabsch
from ..docking.align import _matched_anchors


def select_peptide(structure, chain_id: str, first: int, last: int) -> None:
    """Select a source-verified antigen interval using inclusive author residue numbers.

    Preserve coordinates and original residue indices. Auxiliary peptide chains become
    unassigned; the selected interval is removed from its parent chain and regions.
    Insertion-code boundaries and annotated receptor regions are rejected.
    """
    from ..annotation.chains import _tag_peptide
    from ..structure.model import Chain

    chains = [c for c in structure.chains if c.chain_id == chain_id]
    if len(chains) != 1 or first > last:
        raise ValueError("invalid peptide chain or author-number interval")
    parent = chains[0]
    selected = [r for r in parent.residues if first <= r.pdb_index <= last]
    if (not selected or selected[0].pdb_index != first or selected[-1].pdb_index != last
            or any(r.insertion_code for r in selected)):
        raise ValueError("peptide endpoints absent or insertion codes require explicit selection")
    if parent.chain_type in {"TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"}:
        protected = {r.seq_index for region in parent.regions for r in region.residues}
        if not protected or any(r.seq_index in protected for r in selected):
            raise ValueError("peptide selection overlaps receptor markup or markup is absent")
    for chain in structure.chains:
        if chain.chain_type == "PEPTIDE":
            chain.chain_type = chain.chain_supertype = None
            chain.regions = []
    if len(selected) == len(parent.residues):
        _tag_peptide(parent)
        return
    used = {c.chain_id for c in structure.chains}
    new_id = next((x for x in "pqrstuvwxyz0123456789" if x not in used), None)
    if new_id is None:
        raise ValueError("no free peptide chain identifier")
    indices = {r.seq_index for r in selected}
    parent.residues = [r for r in parent.residues if r.seq_index not in indices]
    for region in parent.regions:
        region.residues = [r for r in region.residues if r.seq_index not in indices]
        region.sequence = "".join(r.aa for r in region.residues)
    peptide = Chain(new_id, selected)
    _tag_peptide(peptide)
    structure.chains.append(peptide)


def structural_register(structure, reference, reference_p1: int) -> pl.DataFrame:
    """Enumerate forward/reverse nine-residue cores in an annotated class-II complex.

    ``reference_p1`` is the verified, one-based observed peptide position of pocket P1.
    Fit both MHC groove chains only; rank candidate peptide Cα RMSDs in that fixed frame.
    Return every candidate and P1/P4/P6/P9 physical identities. No threshold converts
    a minimum into structural truth: residuals, runner-up margin and individual core
    deviations are exported for review. Incomplete nine-residue backbones fail explicitly.
    """
    for s in (structure, reference):
        roles = {c.chain_type for c in s.chains if c.chain_supertype == "MHCII"}
        if not {"MHCa", "MHCb"} <= roles or len(s.by_type("PEPTIDE")) != 1:
            raise ValueError("structural register requires class-II alpha/beta groove and one peptide")
    ref = reference.by_type("PEPTIDE")[0]
    if reference_p1 < 1 or reference_p1 + 8 > len(ref.residues):
        raise ValueError("reference P1 must start an observed nine-residue core")
    target = [r.ca for r in ref.residues[reference_p1 - 1:reference_p1 + 8]]
    if any(x is None for x in target):
        raise ValueError("reference core has missing Cα atoms")
    target = np.asarray(target)
    if not np.isfinite(target).all():
        raise ValueError("nonfinite reference peptide Cα coordinates")
    mob, fixed = _matched_anchors(structure, reference)
    if len(mob) < 3 or not np.isfinite(mob).all() or not np.isfinite(fixed).all():
        raise ValueError("insufficient finite MHC groove correspondences")
    rot, tran, mhc_rmsd = kabsch(mob, fixed)
    chain = structure.by_type("PEPTIDE")[0]
    candidates = []
    for start in range(len(chain.residues) - 8):
        core = chain.residues[start:start + 9]
        if any(r.ca is None for r in core):
            continue
        coords = np.asarray([r.ca for r in core]) @ rot + tran
        if not np.isfinite(coords).all():
            raise ValueError("nonfinite peptide Cα coordinates")
        for direction in (1, -1):
            ordered = core[::direction]
            deviations = np.linalg.norm(coords[::direction] - target, axis=1)
            row = {"pdb.id": structure.pdb_id, "chain.id": chain.chain_id,
                   "reference.id": reference.pdb_id, "reference_p1": reference_p1,
                   "core_start": start + 1 if direction == 1 else start + 9,
                   "direction": direction, "core_sequence": "".join(r.aa for r in ordered),
                   "core_rmsd_A": float(np.sqrt(np.mean(deviations ** 2))),
                   "core_max_deviation_A": float(max(deviations)),
                   "mhc_rmsd_A": mhc_rmsd, "n_mhc_ca": len(mob)}
            row.update({f"rotation.{i}{j}": float(rot[i, j]) for i in range(3) for j in range(3)})
            row.update({f"translation.{i}": float(tran[i]) for i in range(3)})
            for pos in (1, 4, 6, 9):
                residue = ordered[pos - 1]
                row.update({f"P{pos}.peptide_pos": start + pos if direction == 1 else start + 10 - pos,
                            f"P{pos}.residue_index": residue.seq_index,
                            f"P{pos}.pdb_index": residue.pdb_index,
                            f"P{pos}.insertion_code": residue.insertion_code,
                            f"P{pos}.deviation_A": float(deviations[pos - 1])})
            candidates.append(row)
    if not candidates:
        raise ValueError("no complete observed nine-residue peptide backbone")
    candidates.sort(key=lambda r: (r["core_rmsd_A"], r["core_start"], -r["direction"]))
    margin = candidates[1]["core_rmsd_A"] - candidates[0]["core_rmsd_A"]
    return pl.DataFrame([{**r, "rank": i + 1, "runner_up_margin_A": margin}
                         for i, r in enumerate(candidates)], infer_schema_length=None)
