r"""Per-residue attribution for a whole-structure score. 2026-09-04

:func:`tcren.score.peptide_score` decomposes per residue by construction: it is a sum over contacts,
and :func:`tcren.energetics.scoring.position_profile` already returns that sum split by peptide
position. The other four read-outs do not. Each is a function of ~149 *whole-structure* scalars
through a full covariance per class, so no residue owns a share of one and there is nothing to split.

What is defined for all five, and for every channel, is the leave-one-out difference

.. math:: \Delta_i \;=\; L(x) \;-\; L(x_{\setminus i})

the score of the complex minus the score of the same complex with residue *i*'s atoms removed.
Positive :math:`\Delta_i` means the residue carries the score. It answers what a per-residue
confidence colouring answers -- how much does this residue carry -- so one implementation colours a
figure by any read-out rather than one colouring per score.

**The cost, and the rule it obeys.** Chain typing and the MHC call are mmseqs searches. They happen
**once**, here, for the whole structure, and every masked copy is featurised against that one
annotation with ``annotate=False`` (CLAUDE.md 0-mmseqs). Only residues that actually touch an
interface are masked, because a residue making no contact moves no descriptor and its
:math:`\Delta_i` is zero by construction rather than by measurement.

**A masked row that goes non-finite is reported NaN, never imputed.** Removing a residue can leave a
descriptor undefined -- coverage entropy over an empty contact set is the case that occurs -- and
:meth:`tcren.score.ScoreModel.coordinates` drops such a row rather than filling it. That is the
right behaviour and it is surfaced: the residue gets a null delta, not a fabricated one.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..structure.model import Chain, RegionMarkup, Structure

#: The interfaces whose contacts define "touches an interface".
INTERFACES = ("tcr_peptide", "tcr_mhc", "peptide_mhc")


def _drop_residue(s: Structure, chain_id: str, seq_index: int) -> Structure:
    """A shallow copy of ``s`` with one residue's atoms gone, region markup included.

    ``seq_index`` values of the surviving residues are left as they are, so they carry gaps. That is
    deliberate: every consumer here addresses residues by ``seq_index``, and renumbering would move
    the region boundaries that CDR-loop descriptors are cut on.
    """
    chains = []
    for c in s.chains:
        if c.chain_id != chain_id:
            chains.append(c)
            continue
        keep = [r for r in c.residues if r.seq_index != seq_index]
        regions = [RegionMarkup(g.region_type, g.start_seq_index, g.end_seq_index, g.sequence,
                                [r for r in g.residues if r.seq_index != seq_index])
                   for g in c.regions]
        chains.append(Chain(c.chain_id, keep, c.chain_type, c.chain_supertype,
                            c.allele_info, regions))
    # mhc_calls is carried over: the chains are the same chains, so the search result still holds.
    return Structure(s.pdb_id, chains, s.complex_species, s.cell_type, s.mhc_calls)


def _read(name: str, table, model, receptor: bool) -> np.ndarray:
    """One named read-out over a descriptor table. Channels resolve by their own name."""
    from . import binder_score, channel_scores, peptide_score, pose_score
    if name == "binder":
        return binder_score(table, receptor=receptor, model=model)
    if name == "pose":
        return pose_score(table, model=model)
    if name == "peptide":
        return peptide_score(table)
    channels = channel_scores(table, receptor=receptor, model=model)
    if name in channels:
        return channels[name]
    raise ValueError(f"unknown score {name!r}; expected binder, pose, peptide "
                     f"or one of {sorted(channels)}")


def _interface_residues(cm, chain_types, interfaces) -> list[tuple[str, int]]:
    """Every (chain id, seq index) on either side of the requested interfaces, in structure order."""
    seen: dict[tuple[str, int], None] = {}
    for which in interfaces:
        sel = cm.interface(which) if which != "peptide_mhc" else cm.interface(which)
        for side in ("from", "to"):
            rows = sel.select(pl.col(f"chain.id.{side}"), pl.col(f"residue.index.{side}"),
                              pl.col(f"chain.type.{side}")).unique().iter_rows()
            for cid, idx, ctype in rows:
                if chain_types is None or ctype in chain_types:
                    seen[(cid, int(idx))] = None
    return sorted(seen)


def residue_deltas(source, *, score: str = "binder", model=None, organism: str = "human",
                   chain_types: tuple[str, ...] | None = None,
                   interfaces: tuple[str, ...] = INTERFACES,
                   radii: tuple[float, ...] = (7.0, 8.0)) -> pl.DataFrame:
    """Leave-one-residue-out contribution to a whole-structure score.

    Args:
        source: a path to a structure, or an already-parsed :class:`~tcren.structure.Structure`.
        score: ``"binder"``, ``"pose"``, ``"peptide"``, or a channel name (``"placement"``,
            ``"interface"``, ``"shape"``, ``"energetics"``, ``"mechanics"``).
        model: overrides :func:`tcren.score.holdout_model`.
        organism: passed to chain typing when the structure is not already annotated.
        chain_types: restrict to residues of these chain types, e.g. ``("PEPTIDE",)`` to colour a
            peptide only. ``None`` takes every residue on either side of ``interfaces``.
        interfaces: which interfaces define the residue set.
        radii: footprint radii, as in :func:`tcren.descriptors.table.recognition_table`.

    Returns:
        One row per masked residue: ``chain.id``, ``residue.index`` (seq index), ``residue.pdb``,
        ``residue.aa``, ``score.full``, ``score.without`` and ``delta`` = full − without. ``delta``
        is null where masking the residue left a descriptor undefined.
    """
    from ..annotation import classify_chains
    from ..contactmap import ContactMap
    from ..descriptors.table import _featurise_families
    from ..mhc import annotate_mhc
    from ..structure import import_structure

    s = source if isinstance(source, Structure) else import_structure(source)
    if all(c.chain_type is None for c in s.chains):        # the one mmseqs pass, for every mask
        classify_chains(s, organism=organism, autodetect_species=True)
    if s.mhc_calls is None:
        annotate_mhc(s)

    families = ["placement", "interface", "topology", "energetics", "potts", "kinetics"]
    targets = _interface_residues(ContactMap.from_structure(s), chain_types, interfaces)
    if not targets:
        raise ValueError(f"{s.pdb_id}: no residue touches any of {interfaces}")

    rows = [_featurise_families(s.pdb_id, s, organism, families, radii)]
    for cid, idx in targets:
        rows.append(_featurise_families(f"{s.pdb_id}:{cid}{idx}",
                                        _drop_residue(s, cid, idx), organism, families, radii))
    # One scoring call over baseline + masks, so the coordinate set is chosen once and every row is
    # read on identical columns. Scoring row by row would let one mask's missing column silently
    # change which descriptors the comparison is made over.
    vals = np.asarray(_read(score, pl.DataFrame(rows, infer_schema_length=None), model,
                            receptor=chain_types != ("PEPTIDE",)), float)
    full, without = float(vals[0]), vals[1:]
    lookup = {(c.chain_id, r.seq_index): r for c in s.chains for r in c.residues}
    return pl.DataFrame({
        "chain.id": [c for c, _ in targets],
        "residue.index": [i for _, i in targets],
        "residue.pdb": [lookup[t].pdb_index for t in targets],
        "residue.aa": [lookup[t].aa for t in targets],
        "score.full": full,
        "score.without": without,
        "delta": full - without,
    })
