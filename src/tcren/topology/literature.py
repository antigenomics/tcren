"""Published interface descriptors, computed on the objects this package already builds.

Three families the protein-protein interface literature defines and this catalogue did not reach,
each measured against the whole 19,213-structure benchmark corpus before being added here.

**Surface complementarity and the gap.** :mod:`tcren.topology.surface` has shipped
:func:`~tcren.topology.surface.surface_complementarity` since before this module existed and not
one of its twelve outputs was catalogued, so three channels an audit listed as unreachable were in
fact two lines away::

    surface_complementarity(surface_map(s, side="pmhc"), surface_map(s, side="tcr"))

Both faces are rasterised as height fields on **one shared grid** in a groove frame refit from the
structure, so the gap is a subtraction, not new geometry: ``gap = h_tcr - h_pmhc`` per cell, in
Angstrom. The module's own calibration over 60 Native2026 crystals records the fact that makes the
sign readable -- the median gap is **-1.7 A** and 71 % of cells are interdigitated, the receptor's
lowest point in a cell lying below the groove's highest point in the same cell. The two faces
interlock rather than stack, so a gap that grows is a receptor riding on a few high points.

The mean cancels the two signs, so the gap is also **integrated over the contact plane** with
them kept apart: ``sc_gap_vol`` is the void volume, ``sc_interlock`` the interdigitated volume,
both in Angstrom^3, and ``sc_gap_index`` is the void divided by the retained contact area, which is
the intensive form. On a real interface ``sc_interlock`` is the larger of the two.

Four more read the same field **by sign** rather than pooling it, because a cell the receptor
stands off from and one it reaches into are different events: ``sc_interlock_frac`` (the share of
cells that interlock, the per-structure form of the corpus 71 %), ``sc_gap_depth`` (how far in,
over the interlocked cells alone), ``sc_gap_height`` (how far off, over the void cells alone) and
``sc_gap_asym`` (the balance of the two volumes, in [-1, 1]). ``sc_gap_mean`` is
``interlock_frac`` weighting ``gap_height`` against ``gap_depth``, so these are the three numbers
it collapses into one, not a reparameterisation of it.

Measured over 4,907 labelled structures, ``gap_mean`` and ``gap_sd`` carry the largest
binder/non-binder contrasts of any published descriptor tested (Cohen's *d* -0.651 and -0.681
out-of-panel) at an R^2 on all 141 incumbent descriptors of **0.131 and 0.255** -- they are a
channel this catalogue could not otherwise reach. ``sc_shape`` is the more familiar quantity
(Lawrence & Colman's Sc is the same idea on a dot surface) and the *less* novel of the two: R^2
0.445, closest incumbent ``m_erank_tm`` at rho 0.414.

This is the gap **channel** Jones & Thornton 1996 named, not their gap volume index, which is a
gap volume from a Voronoi construction divided by interface ASA. Ours is a raster height-field
gap and is documented as such rather than cited to their formula.

**Contact order.** Plaxco et al. (*J Mol Biol* 1998, 277:985) normalise mean sequence separation
by chain length, and the whole lesson of that paper is that the normalisation is what makes the
quantity useful. The catalogue had no sequence-separation descriptor at all. On the MHC helix it
is genuinely new (R^2 0.275); on the peptide it is not (R^2 0.519), because our peptides are
overwhelmingly 9-mers and the length normalisation has almost no range to work with -- 619
distinct values over 12,662 structures.

**Participation coefficient.** Di Paola et al. (*Front Bioeng Biotechnol* 2015, 3:170) read PPI
interfaces as contact networks and found ``P_i = 1 - sum_s (k_si / k_i)^2`` their most
discriminative descriptor. On this corpus it is **redundant**: R^2 0.730 (TCR side) and 0.720
(pMHC side) on the incumbents, the pMHC side's nearest neighbour being ``g_loop_overlap`` from the
same bipartite object. It is emitted anyway, at no extra cost over the contact map already built,
and :data:`tcren.recognition.STATUS` says what it duplicates.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from ..contactmap import ContactMap
from ..structure.model import Structure

__all__ = [
    "LITERATURE_FEATURES",
    "SURFACE_FEATURES",
    "contact_order",
    "literature_features",
    "participation_coefficient",
    "surface_features",
]

#: The twelve :func:`~tcren.topology.surface.surface_complementarity` outputs, under catalogue
#: names. ``sc_cells`` and ``sc_coverage`` are bookkeeping -- how much surface the comparison
#: actually reached -- and are carried so a caller can tell a low complementarity from a thin one.
_SURF_MAP: dict[str, str] = {
    "shape_r": "sc_shape", "charge_r": "sc_charge", "phobic_r": "sc_phobic",
    "charge_product": "sc_charge_prod", "phobic_product": "sc_phobic_prod",
    "gap_mean": "sc_gap_mean", "gap_sd": "sc_gap_sd",
    # the gap integrated over the contact plane, the two signs apart: void against interlock
    "gap_vol": "sc_gap_vol", "interlock": "sc_interlock", "gap_index": "sc_gap_index",
    # and the same field resolved by sign: how much of it interlocks, how deep, how high the rest
    # stands off, and the balance of the two volumes
    "interlock_frac": "sc_interlock_frac", "gap_depth": "sc_gap_depth",
    "gap_height": "sc_gap_height", "gap_asym": "sc_gap_asym",
    "d_h": "sc_dh", "d_charge": "sc_dcharge", "d_phobic": "sc_dphobic",
    "n_cells": "sc_cells", "coverage": "sc_coverage",
}
SURFACE_FEATURES: tuple[str, ...] = tuple(_SURF_MAP.values())

#: Everything this module emits, in catalogue order.
LITERATURE_FEATURES: tuple[str, ...] = SURFACE_FEATURES + (
    "co_pep", "co_mhc", "partcoef_tcr", "partcoef_pmhc")

_CELL_LOOPS: tuple[str, ...] = ("TRA:CDR1", "TRA:CDR2", "TRA:CDR3",
                                "TRB:CDR1", "TRB:CDR2", "TRB:CDR3")


def surface_features(structure: Structure) -> dict[str, float]:
    """The twelve surface-complementarity quantities, or NaN throughout if the maps cannot be built.

    Args:
        structure: a chain-typed, MHC-annotated TCR-pMHC structure.

    Returns:
        ``{name: value}`` over :data:`SURFACE_FEATURES`. NaN rather than 0 wherever a map is
        unreachable -- an interface whose surface could not be rasterised has no complementarity,
        and a 0 would rank it as perfectly anti-complementary.
    """
    from .surface import surface_complementarity, surface_map

    row = dict.fromkeys(SURFACE_FEATURES, float("nan"))
    try:
        sc = surface_complementarity(surface_map(structure, side="pmhc"),
                                     surface_map(structure, side="tcr"))
    except Exception:  # noqa: BLE001 - an unmappable groove is a NaN row, not a dead batch
        return row
    return {**row, **{v: float(sc.get(k, float("nan"))) for k, v in _SURF_MAP.items()}}


def _contact_frame(structure: Structure, cutoff: float) -> pl.DataFrame | None:
    """Long form of the CDR-loop x pMHC contact set, keeping which target each edge reaches.

    :func:`tcren.topology.graph._biadjacency` builds the same edge set but collapses peptide and
    MHC into one column axis, which is exactly the distinction both descriptors below turn on:
    the participation coefficient partitions a residue's edges by target, and contact order is
    normalised by the target's own length.
    """
    cm = ContactMap.from_structure(structure, cutoff=cutoff)
    frames = []
    for iface, target in (("tcr_peptide", "peptide"), ("tcr_mhc", "mhc")):
        d = cm.interface(iface)
        if d.is_empty():
            continue
        frames.append(d.select(
            pl.concat_str([pl.col("chain.type.from"), pl.col("region.type.from")],
                          separator=":").alias("loop"),
            pl.col("chain.id.from").alias("cf"),
            pl.col("residue.index.from").cast(pl.Int64).alias("rf"),
            pl.col("residue.index.to").cast(pl.Int64).alias("rt"),
            pl.lit(target).alias("target")))
    if not frames:
        return None
    t = (pl.concat(frames, how="vertical")
         .filter(pl.col("loop").is_in(list(_CELL_LOOPS)))
         .unique(subset=["cf", "rf", "target", "rt"]))
    return None if t.is_empty() else t


def participation_coefficient(t: pl.DataFrame) -> dict[str, float]:
    r"""Di Paola's :math:`P_i = 1 - \sum_s (k_{si} / k_i)^2`, averaged over the engaged residues.

    The modules :math:`s` are the two things a TCR residue can touch -- peptide and MHC -- on the
    receptor side, and the six CDR loops on the pMHC side. :math:`k_{si}` is residue *i*'s edge
    count into module *s* and :math:`k_i` its total degree, so :math:`P_i` is 0 for a residue whose
    contacts all land in one module and approaches :math:`1 - 1/S` for one that spreads them
    evenly over all *S*. Higher is more shared.

    Args:
        t: the long contact frame from :func:`_contact_frame`.

    Returns:
        ``{"partcoef_tcr": float, "partcoef_pmhc": float}``, NaN where no residue is engaged.
    """
    out: dict[str, float] = {"partcoef_tcr": float("nan"), "partcoef_pmhc": float("nan")}
    # TCR side: does an engaged loop residue reach the peptide, the MHC, or both?
    a = (t.group_by(["cf", "rf", "target"]).len()
         .group_by(["cf", "rf"]).agg(pl.col("len").sum().alias("k"),
                                     (pl.col("len") ** 2).sum().alias("sq")))
    if a.height:
        out["partcoef_tcr"] = float((1 - a["sq"] / a["k"] ** 2).mean())
    # pMHC side: does an engaged target residue read one CDR loop or several?
    b = (t.group_by(["target", "rt", "loop"]).len()
         .group_by(["target", "rt"]).agg(pl.col("len").sum().alias("k"),
                                         (pl.col("len") ** 2).sum().alias("sq")))
    if b.height:
        out["partcoef_pmhc"] = float((1 - b["sq"] / b["k"] ** 2).mean())
    return out


def contact_order(t: pl.DataFrame) -> dict[str, float]:
    r"""Plaxco's length-normalised sequence spread, per target.

    For each CDR loop, the mean absolute sequence separation :math:`|i - j|` between the distinct
    target residues it reaches; averaged over loops and divided by the target's own span *L* in
    residues, which is the normalisation Plaxco's result rests on. Higher means one loop reaches
    residues far apart in the target sequence.

    Args:
        t: the long contact frame from :func:`_contact_frame`.

    Returns:
        ``{"co_pep": float, "co_mhc": float}``, NaN where the target spans under two residues or
        no loop reaches more than one of its residues.
    """
    out = {"co_pep": float("nan"), "co_mhc": float("nan")}
    for target, key in (("peptide", "co_pep"), ("mhc", "co_mhc")):
        g = t.filter(pl.col("target") == target)
        if not g.height:
            continue
        span = int(g["rt"].max()) - int(g["rt"].min()) + 1
        if span < 2:
            continue
        seps = []
        for (_loop,), h in g.group_by(["loop"]):
            idx = np.unique(h["rt"].to_numpy())
            if len(idx) > 1:
                d = np.abs(idx[:, None] - idx[None, :])[np.triu_indices(len(idx), 1)]
                seps.append(float(d.mean()))
        if seps:
            out[key] = float(np.mean(seps)) / span
    return out


def literature_features(structure: Structure, *, cutoff: float = 5.0) -> dict[str, float]:
    """Every descriptor in this module, as a flat row.

    Args:
        structure: a chain-typed, CDR-region-annotated, MHC-annotated TCR-pMHC structure.
        cutoff: heavy-atom contact threshold in Angstrom for the graph-derived pair. The surface
            pair has its own length scales, fixed in :mod:`tcren.topology.surface`.

    Returns:
        ``{name: value}`` over :data:`LITERATURE_FEATURES`, NaN wherever the structure gives a
        descriptor no support.
    """
    row = dict.fromkeys(LITERATURE_FEATURES, float("nan"))
    row.update(surface_features(structure))
    t = _contact_frame(structure, cutoff)
    if t is not None:
        row.update(participation_coefficient(t))
        row.update(contact_order(t))
    return row


def _selfcheck() -> None:  # pragma: no cover - exercised by tests/unit/test_literature.py
    """The two graph descriptors on a hand-built contact set with a known answer."""
    # Two TCR residues. The first touches one peptide residue and one MHC residue (P = 1 - 1/2 =
    # 0.5); the second touches two peptide residues only (P = 1 - 1 = 0). Mean 0.25.
    t = pl.DataFrame({
        "loop": ["TRA:CDR3"] * 4,
        "cf": ["A"] * 4, "rf": [1, 1, 2, 2], "rt": [3, 7, 3, 5],
        "target": ["peptide", "mhc", "peptide", "peptide"],
    })
    p = participation_coefficient(t)
    assert abs(p["partcoef_tcr"] - 0.25) < 1e-12, p

    # Peptide residues 3 and 5 are reached by one loop, span = 5 - 3 + 1 = 3, mean separation 2.
    c = contact_order(t)
    assert abs(c["co_pep"] - 2.0 / 3.0) < 1e-12, c

    # A frame with no MHC edges leaves the MHC descriptor NaN, never 0.
    only_pep = t.filter(pl.col("target") == "peptide")
    assert np.isnan(contact_order(only_pep)["co_mhc"])
    print("literature self-check OK")


if __name__ == "__main__":  # pragma: no cover
    _selfcheck()
