"""DFIRE reference states at residue granularity, and the corrections they imply for a
plain contact potential such as TCRen.

TCRen counts a residue pair as *in contact* when its closest heavy atoms sit within 5 Å,
and every contact then counts the same. Two pieces of information the coordinates already
carry are discarded by that: **how far apart** the two residues actually are, and **how
they are turned relative to each other**. DFIRE (Zhou & Zhou, *Protein Sci* 2002) supplies
the first, its dipolar successor DFIRE2 (Yang & Zhou, *Protein Sci* 2008) the second.

The DFIRE reference state is the observation that in a *finite* globular system the number
of residue pairs at separation :math:`r` grows as :math:`r^\\alpha` with
:math:`\\alpha \\approx 1.61`, not as the ideal-gas :math:`r^2`. Its energy is

.. math::
    u(a, b, r) = -\\ln \\frac{N(a, b, r)}
                             {(r/r_c)^{\\alpha}\\,(\\Delta r/\\Delta r_c)\\, N(a, b, r_c)}

in units of :math:`RT`, matching the natural-log convention of
:func:`tcren.potential.derive_tcren`. DFIRE2 adds an orientation coordinate: with
:math:`\\hat{u}_a` the Cα→Cβ direction of residue *a* and :math:`\\hat{r}` the interresidue
direction, the pair is described by :math:`(\\cos\\theta_a, \\cos\\theta_b)`, whose isotropic
reference density is uniform.

**Why the corrections are transferable, and the potential is not.** A distance- and
orientation-resolved 20×20 potential needs 20 × 20 × (bins) × (orientation cells)
occupancies. The 374 reference crystals carry roughly 8,000 TCR:peptide contacts, which is
about two observations per cell — no estimate at all. The two *corrections*, by contrast,
are one number per amino-acid pair, and they are properties of packing geometry rather than
of TCR biology, so they can be estimated on every inter-chain residue pair of every
interface in the same crystals and then added to the sparse TCR:peptide derivation. That
is what :func:`corrections` does and what its ``scope`` argument selects.

The decomposition returned by :func:`corrections` is, per amino-acid pair,

``E0``
    the orientation-free DFIRE energy averaged over the pair's observed contacts — the
    stochastic term, the direct analogue of one TCRen cell;
``C_dist``
    the change in TCRen's own log-odds when each contact is weighted by the DFIRE volume
    element :math:`(r/r_c)^{-\\alpha}` instead of counting one — the distance correction;
``C_rot``
    :math:`-\\mathrm{KL}\\!\\left(P(\\cos\\theta_a, \\cos\\theta_b \\mid a, b)\\,\\|\\,
    \\mathrm{uniform}\\right)`, the orientational free energy a contact-only count cannot
    see — the rotation correction. It is ``<= 0`` by construction, and its magnitude is the
    orientational information the pair carries, in nats.

``DFIRE2 = E0 + C_rot`` exactly, so the shipped DFIRE2 matrix and the corrected TCRen are
built from the same three columns.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl
from scipy.spatial import cKDTree

from ..structure.model import Structure
from .model import AA20, Potential
from .._provenance import not_in_tcren2

#: Finite-size radial exponent of the DFIRE reference state (Zhou & Zhou 2002).
ALPHA = 1.61
#: Radial cutoff of the reference state (Å). Contacts beyond it are defined to be at zero energy.
RC = 14.5
#: Radial bin width (Å).
DR = 0.5
#: Contact cutoff (Å) on the closest heavy-atom distance — the TCRen contact definition.
CONTACT_CUTOFF = 5.0
#: Bin edges in cos θ. Three bins per residue gives nine orientation cells per pair.
COS_EDGES = (-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0)
#: Contacts a pair needs before its rotation term is read rather than set to zero. The
#: Miller-Madow correction is first order in 1/N and leaves a positive residual: simulating
#: an isotropic null over the nine cells, its 99th percentile is 0.60 nats at N = 9, 0.13 at
#: N = 50 and 0.029 at N = 200 — the last an order of magnitude below the smallest correction
#: worth reading, against an observed range down to -2.2.
MIN_ORIENTED = 200


@dataclass(frozen=True)
class DfireDecomposition:
    """Per-amino-acid-pair decomposition of the DFIRE2 energy of a contact.

    Attributes:
        table: Long table with ``residue.aa.from``, ``residue.aa.to``, ``E0``, ``C_dist``,
            ``C_rot``, ``n_contacts`` and ``n_oriented``.
        radial: The distance-resolved DFIRE energies ``u(a, b, r)``.
        scope: Which residue pairs the corrections were estimated on.
        n_structures: How many structures contributed.
    """

    table: pl.DataFrame
    radial: pl.DataFrame
    scope: str
    n_structures: int

    def dfire2(self) -> Potential:
        """The DFIRE2 matrix: ``E0 + C_rot``, one energy per amino-acid pair."""
        long = self.table.select(
            "residue.aa.from", "residue.aa.to",
            (pl.col("E0") + pl.col("C_rot")).alias("value"),
        )
        return Potential(name="DFIRE2", matrix=long, alphabet=AA20)


@not_in_tcren2('As pair_geometry.')
def geometry_set(struct_dir, rc: float = RC, on_error: str = "skip") -> tuple[pl.DataFrame, int]:
    """Stack :func:`pair_geometry` over a folder of structures.

    Uses the batched annotation path, so the whole of ``Native2026`` costs one arda call
    rather than 374.

    Args:
        struct_dir: Folder of PDB/mmCIF structures.
        rc: Radial cutoff (Å).
        on_error: Forwarded to the annotation pass.

    Returns:
        ``(geometry, n_structures)``.
    """
    from ..paper.helpers import iter_annotated_set

    frames, n = [], 0
    for s in iter_annotated_set(struct_dir, on_error=on_error):
        g = pair_geometry(s, rc=rc)
        if g.height:
            frames.append(g.with_columns(pl.lit(s.pdb_id).alias("pdb.id")))
            n += 1
    return (pl.concat(frames) if frames else pl.DataFrame()), n


@not_in_tcren2('DFIRE reference states are an independent line of work, evaluated against TCRen2 rather than folded into it.')
def pair_geometry(structure: Structure, rc: float = RC) -> pl.DataFrame:
    """Distance and mutual orientation of every inter-chain residue pair within ``rc``.

    Separation is measured between representative atoms (Cβ, or Cα for glycine and any
    residue missing its Cβ). Orientation is the pair of cosines between each residue's
    Cα→Cβ direction and the interresidue direction; glycine has no such direction and its
    cosine is null, which drops that pair from the rotation term but not from the radial one.

    Args:
        structure: A parsed, chain-typed structure.
        rc: Radial cutoff (Å).

    Returns:
        Columns ``chain.type.from``/``.to``, ``residue.aa.from``/``.to``, ``dist``,
        ``cos.from``, ``cos.to``, ``contact`` (closest heavy atom within
        :data:`CONTACT_CUTOFF`), one row per unordered inter-chain residue pair.
    """
    rep, direction, aa, ctype, cidx = [], [], [], [], []
    for ci, chain in enumerate(structure.chains):
        for res in chain.residues:
            r = res.cb_or_ca
            if r is None:
                continue
            rep.append(r)
            ca, cb = res.ca, res.cb
            direction.append(
                (cb - ca) / np.linalg.norm(cb - ca)
                if ca is not None and cb is not None and np.linalg.norm(cb - ca) > 1e-6
                else np.full(3, np.nan)
            )
            aa.append(res.aa)
            ctype.append(chain.chain_type or chain.chain_id)
            cidx.append(ci)

    schema = {
        "chain.type.from": pl.Utf8, "chain.type.to": pl.Utf8,
        "residue.aa.from": pl.Utf8, "residue.aa.to": pl.Utf8,
        "dist": pl.Float64, "cos.from": pl.Float64, "cos.to": pl.Float64,
        "contact": pl.Boolean,
    }
    if len(rep) < 2:
        return pl.DataFrame(schema=schema)

    rep = np.asarray(rep, dtype=np.float64)
    direction = np.asarray(direction, dtype=np.float64)
    aa = np.asarray(aa, dtype=object)
    ctype = np.asarray(ctype, dtype=object)
    cidx = np.asarray(cidx, dtype=np.int64)

    pairs = cKDTree(rep).query_pairs(r=rc, output_type="ndarray")
    if len(pairs) == 0:
        return pl.DataFrame(schema=schema)
    i, j = pairs[:, 0], pairs[:, 1]
    keep = cidx[i] != cidx[j]
    i, j = i[keep], j[keep]
    if len(i) == 0:
        return pl.DataFrame(schema=schema)

    delta = rep[j] - rep[i]
    dist = np.linalg.norm(delta, axis=1)
    rhat = delta / dist[:, None]
    cos_from = np.einsum("ij,ij->i", direction[i], rhat)
    cos_to = -np.einsum("ij,ij->i", direction[j], rhat)

    return pl.DataFrame(
        {
            "chain.type.from": ctype[i], "chain.type.to": ctype[j],
            "residue.aa.from": aa[i], "residue.aa.to": aa[j],
            "dist": dist, "cos.from": cos_from, "cos.to": cos_to,
            "contact": _heavy_atom_contact(structure, cidx, i, j, rep),
        },
        schema=schema,
    )


def _heavy_atom_contact(
    structure: Structure, cidx: np.ndarray, i: np.ndarray, j: np.ndarray, rep: np.ndarray
) -> np.ndarray:
    """Flag the residue pairs whose closest heavy atoms are within :data:`CONTACT_CUTOFF`.

    Runs one atom-level neighbour search over the whole structure and looks the answer up
    per pair, rather than measuring each pair — the representative-atom pairs already
    within :data:`RC` are the only candidates, and a heavy-atom contact implies a
    representative separation well inside that.
    """
    coords, owner = [], []
    slot = 0
    for chain in structure.chains:
        for res in chain.residues:
            if res.cb_or_ca is None:
                continue
            for atom in res.atoms:
                if atom.element != "H":
                    coords.append(atom.coord)
                    owner.append(slot)
            slot += 1
    if not coords:
        return np.zeros(len(i), dtype=bool)
    owner = np.asarray(owner, dtype=np.int64)
    apairs = cKDTree(np.asarray(coords, dtype=np.float64)).query_pairs(
        r=CONTACT_CUTOFF, output_type="ndarray"
    )
    if len(apairs) == 0:
        return np.zeros(len(i), dtype=bool)
    a, b = owner[apairs[:, 0]], owner[apairs[:, 1]]
    n = len(rep)
    touching = set((int(min(x, y)) * n + int(max(x, y))) for x, y in zip(a, b) if x != y)
    keys = np.minimum(i, j) * n + np.maximum(i, j)
    return np.fromiter((int(k) in touching for k in keys), dtype=bool, count=len(keys))


_TCR = {"TRA", "TRB", "TRD", "TRG"}
# ``classify_chains`` alone types the class-I heavy chain "MHC"; ``annotate_mhc`` refines it
# to MHCa/MHCb. Accept both so the scope filter does not silently return nothing.
_MHC = {"MHC", "MHCa", "MHCb"}
_SCOPES = {
    "tcr_peptide": (_TCR, {"PEPTIDE"}),
    "tcr_mhc": (_TCR, _MHC),
    "peptide_mhc": ({"PEPTIDE"}, _MHC),
}


@not_in_tcren2('As corrections.')
def select_scope(geom: pl.DataFrame, scope: str = "all") -> pl.DataFrame:
    """Restrict a :func:`pair_geometry` table to one interface, orienting each pair.

    ``"all"`` keeps every inter-chain pair as measured. A named interface additionally
    swaps the two sides where needed so that ``from`` is always the first chain class of
    the pair — otherwise a TCR:peptide contact would be counted in whichever orientation
    the residue indices happened to fall, and the derived matrix would be half-symmetrised
    by accident.
    """
    if scope == "all":
        return geom
    if scope not in _SCOPES:
        raise ValueError(f"scope must be 'all' or one of {sorted(_SCOPES)}, got {scope!r}")
    left, right = _SCOPES[scope]
    fwd = pl.col("chain.type.from").is_in(list(left)) & pl.col("chain.type.to").is_in(list(right))
    rev = pl.col("chain.type.from").is_in(list(right)) & pl.col("chain.type.to").is_in(list(left))
    flip = {f"{s}.from": f"{s}.to" for s in ("chain.type", "residue.aa", "cos")}
    flip |= {v: k for k, v in flip.items()}
    swapped = geom.filter(rev).rename(flip).select(geom.columns)
    return pl.concat([geom.filter(fwd), swapped])


@not_in_tcren2('As pair_geometry.')
def radial_potential(
    geom: pl.DataFrame, rc: float = RC, dr: float = DR, pseudocount: float = 1.0
) -> pl.DataFrame:
    """Distance-resolved DFIRE energies ``u(a, b, r)`` from a pair-geometry table.

    The reference occupancy of a bin is the last bin's count scaled by
    :math:`(r/r_c)^{\\alpha}`, which is the finite-size correction to the ideal-gas
    :math:`r^2`; the bin widths are uniform, so their ratio drops out.

    Returns:
        Columns ``residue.aa.from``, ``residue.aa.to``, ``r`` (bin centre), ``n``,
        ``n_ref``, ``u``.
    """
    edges = np.arange(0.0, rc + dr / 2, dr)
    binned = geom.filter(pl.col("dist") < rc).with_columns(
        ((pl.col("dist") / dr).floor() * dr + dr / 2).alias("r")
    )
    counts = binned.group_by("residue.aa.from", "residue.aa.to", "r").agg(
        pl.len().alias("n")
    )
    grid = pl.DataFrame(
        [(a, b, float(r)) for a in AA20 for b in AA20 for r in (edges[:-1] + dr / 2)],
        schema=["residue.aa.from", "residue.aa.to", "r"],
        orient="row",
    )
    merged = grid.join(counts, on=["residue.aa.from", "residue.aa.to", "r"], how="left").with_columns(
        pl.col("n").fill_null(0).cast(pl.Float64) + pseudocount
    )
    last = float(edges[-2] + dr / 2)
    merged = merged.with_columns(
        pl.col("n").filter(pl.col("r") == last).first()
        .over("residue.aa.from", "residue.aa.to").alias("n_rc")
    )
    return merged.with_columns(
        (pl.col("n_rc") * (pl.col("r") / last) ** ALPHA).alias("n_ref")
    ).with_columns(
        (-(pl.col("n") / pl.col("n_ref")).log()).alias("u")
    ).drop("n_rc")


def _logodds(counts: pl.DataFrame, pseudocount: float = 1.0) -> pl.DataFrame:
    """TCRen's log-odds over a weighted 20×20 pair-count table, on the same convention.

    Mirrors :func:`tcren.potential.derive_tcren` exactly — natural log, marginal product
    reference, one pseudocount — so a difference of two of these is a pure re-weighting
    effect and carries no formula change.
    """
    grid = pl.DataFrame(
        [(a, b) for a in AA20 for b in AA20],
        schema=["residue.aa.from", "residue.aa.to"], orient="row",
    )
    merged = (
        grid.join(counts, on=["residue.aa.from", "residue.aa.to"], how="left")
        .with_columns(pl.col("count").fill_null(0.0) + pseudocount)
        .with_columns(
            pl.col("count").sum().over("residue.aa.from").alias("total.from"),
            pl.col("count").sum().over("residue.aa.to").alias("total.to"),
            pl.col("count").sum().alias("total"),
        )
    )
    return merged.select(
        "residue.aa.from", "residue.aa.to",
        (-(pl.col("count") * pl.col("total") / pl.col("total.to") / pl.col("total.from")).log())
        .alias("value"),
    )


def _rotation_term(
    contacts: pl.DataFrame, edges=COS_EDGES, min_oriented: int = MIN_ORIENTED
) -> pl.DataFrame:
    """``C_rot = -KL(observed orientation ‖ isotropic)`` per amino-acid pair, in ``RT``.

    The isotropic reference is uniform in cos θ, so with equal-width bins in cos θ the
    reference cell probability is ``1/K`` and the divergence is ``ln K - H``. Glycine has no
    Cβ and hence no direction, so its cosines are not finite and its contacts are excluded
    here — ``n_oriented`` records how many survived. The plug-in
    entropy is biased low, hence the divergence biased high, by ``(K-1)/(2N)`` to first
    order; the Miller–Madow term subtracts that, and the result is clipped at zero because
    a negative divergence is an estimation artefact rather than a signal. Below
    ``min_oriented`` contacts the residual is comparable to the effect and the term is set to
    zero outright, so a thin cell is not credited with orientational order it has not shown.
    """
    k = len(edges) - 1
    cells = k * k
    cut = [pl.col(c).cut(list(edges[1:-1]), labels=[str(x) for x in range(k)])
           for c in ("cos.from", "cos.to")]
    per_cell = (
        contacts.filter(pl.col("cos.from").is_finite() & pl.col("cos.to").is_finite())
        .with_columns(cut[0].alias("bf"), cut[1].alias("bt"))
        .group_by("residue.aa.from", "residue.aa.to", "bf", "bt")
        .agg(pl.len().alias("n"))
    )
    return (
        per_cell.with_columns(pl.col("n").sum().over("residue.aa.from", "residue.aa.to").alias("N"))
        .with_columns((pl.col("n") / pl.col("N")).alias("p"))
        .group_by("residue.aa.from", "residue.aa.to")
        .agg(
            (pl.lit(float(np.log(cells))) + (pl.col("p") * pl.col("p").log()).sum()).alias("kl_raw"),
            pl.col("N").first().alias("n_oriented"),
        )
        .with_columns(
            pl.when(pl.col("n_oriented") >= min_oriented)
            .then(pl.col("kl_raw") - (cells - 1) / (2 * pl.col("n_oriented")))
            .otherwise(0.0)
            .clip(lower_bound=0.0).alias("kl")
        )
        .select("residue.aa.from", "residue.aa.to", (-pl.col("kl")).alias("C_rot"), "n_oriented")
    )


@not_in_tcren2('The distance and rotation corrections are measured against TCRen2, not applied to the shipped matrix. Which interface they are estimated on decides their sign: transferred from peptide:MHC they improve CPL, pooled over all interfaces they harm it.')
def corrections(
    geom: pl.DataFrame, scope: str = "all", n_structures: int = 0, rc: float = RC,
    min_oriented: int = MIN_ORIENTED,
) -> DfireDecomposition:
    """Decompose the DFIRE2 energy of a contact into ``E0``, ``C_dist`` and ``C_rot``.

    Args:
        geom: Stacked :func:`pair_geometry` output over a structure set.
        scope: Interface to estimate on — ``"all"`` (every inter-chain pair, the widest
            sample and the default), or ``tcr_peptide`` / ``tcr_mhc`` / ``peptide_mhc``.
        n_structures: Recorded on the result for provenance.
        rc: Radial cutoff (Å), which must match the one :func:`pair_geometry` used.
        min_oriented: Contacts a pair needs before its rotation term is read (see
            :data:`MIN_ORIENTED`).

    Returns:
        A :class:`DfireDecomposition`. Pairs with no observed contact carry zero for all
        three terms, so adding the corrections to a potential leaves those cells untouched
        rather than moving them on no evidence. That has to be forced for ``C_dist``: it is
        a difference of two pseudocounted log-odds, which is non-zero even where nothing
        was seen.
    """
    geom = select_scope(geom, scope)
    radial = radial_potential(geom, rc=rc)
    contacts = geom.filter(pl.col("contact"))

    e0 = (
        contacts.with_columns(((pl.col("dist") / DR).floor() * DR + DR / 2).alias("r"))
        .join(radial.select("residue.aa.from", "residue.aa.to", "r", "u"),
              on=["residue.aa.from", "residue.aa.to", "r"], how="left")
        .group_by("residue.aa.from", "residue.aa.to")
        .agg(pl.col("u").mean().alias("E0"), pl.len().alias("n_contacts"))
    )

    flat = _logodds(contacts.group_by("residue.aa.from", "residue.aa.to")
                    .agg(pl.len().cast(pl.Float64).alias("count")))
    weighted = _logodds(
        contacts.with_columns(((pl.col("dist") / rc) ** (-ALPHA)).alias("w"))
        .group_by("residue.aa.from", "residue.aa.to")
        .agg(pl.col("w").sum().alias("count"))
    )
    c_dist = weighted.join(flat, on=["residue.aa.from", "residue.aa.to"], suffix="_flat").select(
        "residue.aa.from", "residue.aa.to",
        (pl.col("value") - pl.col("value_flat")).alias("C_dist"),
    )

    grid = pl.DataFrame([(a, b) for a in AA20 for b in AA20],
                        schema=["residue.aa.from", "residue.aa.to"], orient="row")
    table = (
        grid.join(e0, on=["residue.aa.from", "residue.aa.to"], how="left")
        .join(c_dist, on=["residue.aa.from", "residue.aa.to"], how="left")
        .join(_rotation_term(contacts, min_oriented=min_oriented),
              on=["residue.aa.from", "residue.aa.to"], how="left")
        .with_columns(
            pl.col("E0").fill_null(0.0), pl.col("C_dist").fill_null(0.0),
            pl.col("C_rot").fill_null(0.0), pl.col("n_contacts").fill_null(0),
            pl.col("n_oriented").fill_null(0),
        )
        .with_columns(
            pl.when(pl.col("n_contacts") > 0).then(pl.col("C_dist")).otherwise(0.0)
            .alias("C_dist")
        )
    )
    return DfireDecomposition(table=table, radial=radial, scope=scope, n_structures=n_structures)


@not_in_tcren2('As corrections.')
def apply_corrections(
    potential: Potential, decomposition: DfireDecomposition,
    terms: tuple[str, ...] = ("dist", "rot"), name: str | None = None,
) -> Potential:
    """Add the DFIRE corrections to an existing contact potential.

    ``C_dist`` re-references the pair's contacts by the finite-size volume element and
    ``C_rot`` credits its orientational specificity; both are estimated on the wide
    inter-chain sample, so a sparse TCR:peptide derivation gains them without needing the
    occupancies itself.

    Args:
        potential: The potential to correct, e.g. the shipped TCRen2 matrix.
        decomposition: Output of :func:`corrections`.
        terms: Which corrections to add — any of ``"dist"``, ``"rot"``.
        name: Name of the result (default: the input's, suffixed with the terms applied).

    Returns:
        A new :class:`Potential` over the input's own cells; cells the decomposition does
        not cover are returned unchanged.
    """
    cols = {"dist": "C_dist", "rot": "C_rot"}
    unknown = set(terms) - set(cols)
    if unknown:
        raise ValueError(f"unknown correction term(s): {sorted(unknown)}")
    delta = decomposition.table.select(
        "residue.aa.from", "residue.aa.to",
        sum((pl.col(cols[t]) for t in terms), start=pl.lit(0.0)).alias("_d"),
    )
    long = (
        potential.matrix.join(delta, on=["residue.aa.from", "residue.aa.to"], how="left")
        .with_columns((pl.col("value") + pl.col("_d").fill_null(0.0)).alias("value"))
        .drop("_d")
    )
    return Potential(
        name=name or f"{potential.name}+{'+'.join(terms)}",
        matrix=long, alphabet=potential.alphabet,
    )
