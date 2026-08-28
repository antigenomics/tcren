"""Available residue pairs — the sites of the contact-map model, and their integer codes.

The reference state of this model is the set of pairs that **could** have contacted. A receptor
residue and a partner residue whose Cα atoms lie within ``radius`` are a *site*; whether a
heavy-atom contact formed within ``cutoff`` is the response. That is what lets the fields read
"the backbone put you in reach — did your side chain engage?", and it is why the availability mask
is **Cα-only**: a side-chain-aware or residue-type-aware radius would make the field circular.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from ..structure.model import Structure
from .model import AA, CLASSES, DBIN, PottsModel, REGIONS, ROLES

#: MHC chain types carrying a groove. B2M is out: it has no groove region and sits under the
#: floor rather than under the receptor.
MHC_PARTNER: tuple[str, ...] = ("MHCa", "MHCb")

#: The groove regions ``tcren.mhc.regions`` projects onto a mapped MHC chain. A residue outside
#: them has no within-region coordinate, so it cannot carry a sequence offset and is dropped.
GROOVE_REGIONS: tuple[str, ...] = ("HELIX_A1", "HELIX_A2", "HELIX_B1", "GROOVE_FLOOR")

SITE_COLUMNS: tuple[str, ...] = (
    "pdb.id", "aa.rec", "chain.rec", "region.rec", "pos.rec",
    "aa.par", "pos.par", "par.len", "role.par", "cls", "d_heavy", "d_ca", "sigma",
)


def _sided(structure: Structure, partner: tuple[str, ...]) -> pl.DataFrame | None:
    """``pose._interface_layers`` plus the normalised partner key and both annotation joins."""
    from ..contacts.table import residue_annotation
    from ..pose import _interface_layers

    w = _interface_layers(structure, 5.0, partner=partner)
    if w.is_empty() or "d3" not in w.columns:
        return None
    # The partner-side key is whichever end of the pair the normalised receptor key is not. Chain
    # ids are distinct across a complex and the layers are inter-chain, so comparing chain suffices.
    rec_is_from = pl.col("chain.id.from") == pl.col("key.tcr.chain")
    w = w.with_columns(
        pl.when(rec_is_from).then(pl.col("chain.id.to"))
          .otherwise(pl.col("chain.id.from")).alias("key.par.chain"),
        pl.when(rec_is_from).then(pl.col("residue.index.to"))
          .otherwise(pl.col("residue.index.from")).alias("key.par.res"),
    )
    ra = residue_annotation(structure)
    rec_ann = ra.select(
        pl.col("chain.id").alias("key.tcr.chain"), pl.col("residue.index").alias("key.tcr.res"),
        pl.col("chain.type").alias("chain.rec"), pl.col("region.type").alias("region.rec"),
        pl.col("region.start").alias("region.start.rec"))
    par_ann = ra.select(
        pl.col("chain.id").alias("key.par.chain"), pl.col("residue.index").alias("key.par.res"),
        pl.col("region.type").alias("region.par"), pl.col("region.start").alias("region.start.par"))
    return (w.join(rec_ann, on=["key.tcr.chain", "key.tcr.res"], how="left")
             .join(par_ann, on=["key.par.chain", "key.par.res"], how="left")
             .with_columns(pl.lit(structure.pdb_id).alias("pdb.id"),
                           (pl.col("key.tcr.res") - pl.col("region.start.rec")).alias("pos.rec")))


def available_pairs(structure: Structure, partner: str = "peptide", *,
                    radius: float = 15.0, cutoff: float = 5.0) -> pl.DataFrame:
    """Every receptor:partner residue pair inside ``radius``, with whether it contacted.

    Args:
        structure: A chain-typed structure. For ``partner="mhc"`` it must also carry MHC groove
            regions — run :func:`tcren.mhc.annotate_mhc` (or ``annotate_mhc_batch`` over a set)
            first, or every MHC residue is dropped for want of a groove region.
        partner: ``"peptide"`` or ``"mhc"``.
        radius: Availability radius, Å, on the Cα–Cα distance.
        cutoff: Contact definition, Å, on the closest heavy-atom distance.

    Returns:
        One row per site: ``pdb.id``, the two residue types and positions, the receptor chain and
        region, the partner role, the partner class, both distances and ``sigma``. ``pos.par`` is a
        **within-region** offset on both partners (the peptide chain carries one region starting at
        0, so there it is the plain 0-based peptide position), because the coupling kernel is
        defined on sequence offsets inside one loop or one helix.

    Example:
        >>> from tcren.potts import available_pairs           # doctest: +SKIP
        >>> pairs = available_pairs(structure)                # doctest: +SKIP
        >>> pairs["sigma"].mean()                             # doctest: +SKIP
    """
    if partner not in ("peptide", "mhc"):
        raise ValueError(f"partner must be 'peptide' or 'mhc', got {partner!r}")
    w = _sided(structure, ("PEPTIDE",) if partner == "peptide" else MHC_PARTNER)
    if w is None:
        return pl.DataFrame(schema={c: pl.Utf8 for c in SITE_COLUMNS})
    if partner == "peptide":
        peptide = next((c.sequence() for c in structure.chains
                        if c.chain_type == "PEPTIDE"), "")
        try:
            from ..refine.anchors import predict_anchors
            anchors = set(predict_anchors(peptide, structure=structure).anchors)
        except Exception:                       # no anchor call: every position is TCR-facing
            anchors = set()
        w = w.with_columns(
            pl.col("key.par.res").alias("pos.par"),
            pl.lit(len(peptide), dtype=pl.Int64).alias("par.len"),
            pl.when(pl.col("key.par.res").is_in(list(anchors))).then(pl.lit("anchor"))
              .otherwise(pl.lit("tcr_facing")).alias("role.par"),
            pl.lit("peptide").alias("cls"))
    else:
        w = w.filter(pl.col("region.par").is_in(list(GROOVE_REGIONS)))
        if w.is_empty():
            return pl.DataFrame(schema={c: pl.Utf8 for c in SITE_COLUMNS})
        lens = {r.region_type: len(r.residues) for c in structure.chains for r in c.regions
                if r.region_type in GROOVE_REGIONS}
        w = w.with_columns(
            (pl.col("key.par.res") - pl.col("region.start.par")).alias("pos.par"),
            pl.col("region.par").replace_strict(lens, default=0)
              .cast(pl.Int64).alias("par.len"),
            pl.col("region.par").alias("role.par"),
            pl.lit("mhc").alias("cls"))
    return (w.rename({"aa.tcr": "aa.rec", "aa.pep": "aa.par",
                      "d1": "d_heavy", "d3": "d_ca"})
            .with_columns((pl.col("d_heavy").is_not_null() & (pl.col("d_heavy") <= cutoff))
                          .cast(pl.Float64).alias("sigma"))
            .filter(pl.col("d_ca").is_not_null() & (pl.col("d_ca") <= radius)
                    & pl.col("aa.rec").is_in(list(AA)) & pl.col("aa.par").is_in(list(AA)))
            .select(list(SITE_COLUMNS))
            .sort("pdb.id", maintain_order=True))


def site_codes(sites: pl.DataFrame, model: PottsModel | None = None, *,
               radius: float = 15.0, dbin: float = DBIN):
    """Integer code arrays per one-body block, plus the block sizes and the annotated frame.

    Bin edges are **global** (base 0), never derived from the frame's own minimum: a parameter
    vector fitted on one structure set has to index the same bins on every other, or every score is
    read against the wrong coefficients.

    Returns:
        ``(codes, sizes, frame)`` — ``codes`` is one integer array per block in the order
        receptor residue, partner residue, distance bin, receptor region, partner role, partner
        class; ``sizes`` the number of levels of each; ``frame`` the input with ``loop``,
        ``pchain`` and a global ``sid`` added.
    """
    aa = tuple(model.alphabet) if model else AA
    regions = tuple(model.regions) if model else REGIONS
    roles = tuple(model.roles) if model else ROLES
    classes = tuple(model.classes) if model else CLASSES
    if model is not None:
        radius, dbin = model.radius, model.dbin
    aa_ix = {a: i for i, a in enumerate(aa)}

    q = sites.with_columns(
        (pl.col("chain.rec") + ":" + pl.col("region.rec")).alias("loop")
    ).with_columns(
        pl.when(pl.col("loop").is_in(list(regions))).then(pl.col("loop"))
          .otherwise(pl.lit("other")).alias("loop"),
        # The partner-chain key: one chain on the peptide arm, the groove region on the MHC arm,
        # which is the unit `pos.par` is numbered within.
        pl.when(pl.col("cls") == "peptide").then(pl.lit("PEP"))
          .otherwise(pl.col("role.par")).alias("pchain"),
    ).with_row_index("sid")

    nb = int(np.floor(radius / dbin)) + 1
    codes = [
        np.array([aa_ix[a] for a in q["aa.rec"]], dtype=np.int64),
        np.array([aa_ix[b] for b in q["aa.par"]], dtype=np.int64),
        np.clip(np.floor(q["d_ca"].to_numpy() / dbin).astype(np.int64), 0, nb - 1),
        np.array([regions.index(r) for r in q["loop"]], dtype=np.int64),
        np.array([roles.index(r) for r in q["role.par"]], dtype=np.int64),
        np.array([classes.index(c) for c in q["cls"]], dtype=np.int64),
    ]
    sizes = [len(aa), len(aa), nb, len(regions), len(roles), len(classes)]
    return codes, sizes, q


def eta(codes, model: PottsModel) -> np.ndarray:
    """The one-body log-odds of a contact at every site, under ``model``."""
    J = model.coupling_array()
    blocks = [np.asarray(model.h_rec), np.asarray(model.h_par), np.asarray(model.g_dist),
              np.asarray(model.g_region), np.asarray(model.g_role), np.asarray(model.g_class)]
    out = np.full(len(codes[0]), float(model.alpha))
    for ck, b in zip(codes, blocks):
        if ck.size and int(ck.max()) >= len(b):
            raise ValueError("a site falls outside the model's blocks — the model was fitted with "
                             "a different radius, alphabet or level set")
        out += b[ck]
    return out + J[codes[0], codes[1]]
