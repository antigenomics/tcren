"""Derivation of the TCRen statistical potential from observed contact maps.

This is a direct port of the R derivations in ``code_paper/2_TCRen_derivation.Rmd``
(``variant="classic"``) and ``tcren_am/tcren_am.Rmd`` (``variant="am"``). The classic
variant reproduces ``TCRen_potential.csv``; the alignment-matrix variant reproduces
``tcren_am/tcren.txt``.
"""

from __future__ import annotations

from itertools import product

import polars as pl

from .model import AA20, AA21, Potential
from .._provenance import not_in_tcren2

# Classic derivation enumerates the 20 standard amino acids (Cys included in the
# grid; dropped from the *from* axis only after the log-odds are computed) — the
# same tuple/order as ``model.AA20``.


def symmetrize_counts(counts: pl.DataFrame) -> pl.DataFrame:
    """Fold a directed aa-pair count table onto its transpose: ``N + Nᵀ``.

    TCRen counts are **directed** — ``from`` is a TCR residue and ``to`` a peptide residue — so
    ``N[a,b]`` and ``N[b,a]`` are different observations and the derived matrix is asymmetric.
    Adding the transpose treats each contact as an *unordered* pair, which is the convention
    Miyazawa–Jernigan uses. Diagonal cells double, as they must: a C–C contact is one unordered
    pair observed from both sides.

    Symmetrising here — on the **raw counts, before the log-odds** — is not the same as averaging
    the finished potential. The marginals (``total.from`` / ``total.to``) are recomputed from the
    folded counts, so the *expected* term of the log-odds changes too; averaging the energies
    afterwards leaves the asymmetric background in place. On the Native2026 derivation set the two
    disagree by 0.29 on average (max 0.82), so the distinction is not cosmetic.

    **Cysteine.** The classic directed derivation drops ``from == "C"`` because free Cys is
    essentially absent from CDR loops — on Native2026 only **4 of 8062** contacts (0.05 %) have a
    TCR-side Cys, against 32 (0.40 %) on the peptide side. Folding *grafts* those peptide-side
    observations onto the Cys row instead of discarding the column, so the symmetric matrix keeps
    a full 20×20 alphabet at no cost: the row that would have been dropped for having no data
    inherits the data the other axis did have.

    Args:
        counts: Long table with ``residue.aa.from``, ``residue.aa.to`` and ``count``.

    Returns:
        The folded table, with one row per unordered pair-cell (still stored in both
        orientations, so it is a full symmetric matrix).

    Example:
        >>> import polars as pl
        >>> c = pl.DataFrame({"residue.aa.from": ["A"], "residue.aa.to": ["W"], "count": [3.0]})
        >>> out = symmetrize_counts(c)
        >>> sorted((r["residue.aa.from"], r["residue.aa.to"], r["count"]) for r in out.iter_rows(named=True))
        [('A', 'W', 3.0), ('W', 'A', 3.0)]
    """
    swapped = counts.select(
        pl.col("residue.aa.to").alias("residue.aa.from"),
        pl.col("residue.aa.from").alias("residue.aa.to"),
        pl.col("count"),
    )
    return (
        pl.concat([counts.select(swapped.columns), swapped])
        .group_by("residue.aa.from", "residue.aa.to")
        .agg(pl.col("count").sum())
    )


def derive_tcren(
    contacts: pl.DataFrame,
    include: list[str] | None = None,
    exclude: list[str] | None = None,
    pseudocount: int = 1,
    variant: str = "classic",
    beta: float = 44.0,
    drop_cys: bool | None = None,
    weights: dict[str, float] | None = None,
    weight_col: str | None = None,
    symmetric: bool = False,
    smooth_beta: float = 0.0,
    smooth_matrix: str = "BLOSUM62",
    impute_min_count: int = 0,
    impute_donors: int = 1,
) -> Potential:
    """Derive a TCRen potential from a table of residue contacts.

    Args:
        contacts: Long table of TCR↔peptide contacts with at least
            ``residue.aa.from``, ``residue.aa.to`` and (for filtering) ``pdb.id``.
        include: If given, keep only contacts whose ``pdb.id`` is in this list.
        exclude: If given, drop contacts whose ``pdb.id`` is in this list.
        pseudocount: Added to every amino-acid pair count (default 1).
        variant: ``"classic"`` (natural-log log-odds over 20 aa, Cys dropped from the
            "from" axis) or ``"am"`` (log2/``beta`` over 21 symbols including a gap,
            Cys retained).
        beta: Temperature divisor used by the ``"am"`` variant.
        drop_cys: Override the per-variant default for dropping ``from == "C"`` rows.
            Forced to ``False`` when ``symmetric`` is set (dropping one axis would
            un-symmetrise the result).
        weights: Optional per-structure weights ``{pdb.id: weight}``. When given, each
            structure's contributions to the aa-pair counts are multiplied by its weight
            (rows whose ``pdb.id`` is absent from the map default to weight ``1.0``);
            this down-weights redundancy while keeping all data (see
            :func:`tcren.potential.redundancy.cluster_weights`). ``None`` (default) is
            unweighted and byte-identical to the legacy derivation.
        weight_col: Name of a **per-contact** weight column in ``contacts``, multiplied with the
            per-structure ``weights``. This is how a contact is down-weighted rather than dropped:
            excluding backbone-only pairs outright removes 46 % of the observations and empties 69
            of 380 cells, whereas giving them a fractional vote keeps every cell populated. A
            contact whose two residues are merely co-located, and will only sample an interacting
            geometry some of the time, is exactly a fractional observation.
        symmetric: Fold the raw counts onto their transpose (:func:`symmetrize_counts`)
            before the log-odds, yielding a **symmetric** ``value[a,b] == value[b,a]``
            potential over an unordered amino-acid pair — the same convention as the
            bundled Miyazawa–Jernigan matrix, and therefore directly comparable to it.
            Default ``False`` keeps the directed TCR→peptide potential, which is the
            shipped ``TCRen_potential.csv``.
        smooth_beta: Substitution-matrix pseudocount weight
            (:func:`tcren.potential.smoothing.smooth_counts`), applied to the pair counts before
            the log-odds. A cell holding ``smooth_beta`` observations is pulled halfway to the
            prior its chemically similar cells imply; a well-observed cell is left alone. This is
            aimed at the rare residues -- tryptophan, cysteine, methionine -- whose cells are
            otherwise set by the flat pseudocount. ``0.0`` (default) is off and byte-identical to
            the unsmoothed derivation.
        smooth_matrix: Substitution matrix behind that prior, and behind the imputation.
        impute_min_count: Rebuild cells holding fewer than this many observations from their
            nearest substitutable neighbours
            (:func:`tcren.potential.smoothing.impute_thin_cells`), leaving every other cell
            untouched. ``0`` disables it. Applied after ``smooth_beta`` when both are given, so
            the imputation sees the smoothed counts. NOT USED FOR TCRen2.
        impute_donors: How many nearest donor cells that imputation averages over.

    Returns:
        The derived :class:`Potential`. For ``"am"`` the long matrix additionally
        carries a ``count`` column.
    """
    if variant not in ("classic", "am"):
        raise ValueError(f"unknown variant {variant!r}")

    df = contacts
    if include is not None:
        df = df.filter(pl.col("pdb.id").is_in(include))
    if exclude is not None:
        df = df.filter(~pl.col("pdb.id").is_in(exclude))

    alphabet = AA20 if variant == "classic" else AA21
    if drop_cys is None:
        drop_cys = variant == "classic"
    if symmetric:
        drop_cys = False  # dropping the "from" Cys row would break the symmetry we just built

    if weights is None and weight_col is None:
        # Unweighted: one row = one count (byte-identical to the legacy path).
        n_contacts = df.height
        counts = df.group_by("residue.aa.from", "residue.aa.to").agg(
            pl.len().alias("count")
        )
    else:
        # Weighted: each row contributes its structure's weight times its own (default 1.0 each).
        w = (df["pdb.id"].replace_strict(weights, default=1.0, return_dtype=pl.Float64)
             if weights is not None else pl.Series([1.0] * df.height))
        if weight_col is not None:
            w = w * df[weight_col].cast(pl.Float64)
        df = df.with_columns(w.alias("_w"))
        n_contacts = float(df["_w"].sum())
        counts = df.group_by("residue.aa.from", "residue.aa.to").agg(
            pl.col("_w").sum().alias("count")
        )
    if symmetric:
        # Fold before the log-odds so the marginals — and hence the expected term — are
        # recomputed from the folded counts. Each contact now appears twice, so the contact
        # total doubles alongside them.
        counts = symmetrize_counts(counts)
        n_contacts = n_contacts * 2
    if smooth_beta or impute_min_count:
        if variant != "classic":
            raise ValueError("smoothing is defined over the 20 standard residues (variant='classic')")
        from .smoothing import impute_thin_cells, smooth_counts
        if smooth_beta:
            counts = smooth_counts(counts, beta=smooth_beta, matrix=smooth_matrix)
        if impute_min_count:
            counts = impute_thin_cells(counts, min_count=impute_min_count,
                                       donors=impute_donors, matrix=smooth_matrix)
    if variant == "am":
        # The gap/gap cell is seeded with the total number of contacts, mirroring the
        # rbind(tibble("-","-", count = nrow(res))) line in tcren_am.Rmd.
        counts = pl.concat(
            [
                counts,
                pl.DataFrame(
                    {"residue.aa.from": ["-"], "residue.aa.to": ["-"], "count": [n_contacts]}
                ).with_columns(pl.col("count").cast(counts["count"].dtype)),
            ]
        )

    grid = pl.DataFrame(
        list(product(alphabet, alphabet)),
        schema=["residue.aa.from", "residue.aa.to"],
        orient="row",
    )
    merged = (
        grid.join(counts, on=["residue.aa.from", "residue.aa.to"], how="left")
        .with_columns(pl.col("count").fill_null(0) + pseudocount)
        .with_columns(
            pl.col("count").sum().over("residue.aa.from").alias("total.from"),
            pl.col("count").sum().over("residue.aa.to").alias("total.to"),
            pl.col("count").sum().alias("total"),
        )
    )

    odds = (
        pl.col("count") * pl.col("total") / pl.col("total.to") / pl.col("total.from")
    )
    if variant == "classic":
        value = -odds.log()
    else:
        value = -odds.log(base=2) / beta
    merged = merged.with_columns(value.alias("TCRen"))

    if drop_cys:
        merged = merged.filter(pl.col("residue.aa.from") != "C")

    out_cols = ["residue.aa.from", "residue.aa.to", "TCRen"]
    if variant == "am":
        out_cols.append("count")  # the am table keeps observed counts alongside energies
    long = merged.select(out_cols).rename({"TCRen": "value"})

    out_alphabet = AA20 if variant == "classic" else AA21
    out_alphabet = tuple(
        a for a in out_alphabet if a in set(long["residue.aa.from"]) | set(long["residue.aa.to"])
    )
    return Potential(name="TCRen", matrix=long, alphabet=out_alphabet)


@not_in_tcren2('Leave-one-out derivation, for testing how much any single structure moves the matrix. Diagnostic, not a production path.')
def derive_tcren_loo(
    contacts: pl.DataFrame,
    pdb_ids: list[str],
    **kwargs,
) -> pl.DataFrame:
    """Leave-one-out TCRen: derive once per structure, excluding it each time.

    Args:
        contacts: Contact table (see :func:`derive_tcren`).
        pdb_ids: Structures to leave out one at a time (also the inclusion set).
        **kwargs: Forwarded to :func:`derive_tcren`.

    Returns:
        Long table ``residue.aa.from, residue.aa.to, TCRen.LOO, pdb.id`` stacking the
        per-structure potentials.
    """
    frames = []
    for pid in pdb_ids:
        pot = derive_tcren(contacts, include=pdb_ids, exclude=[pid], **kwargs)
        frames.append(
            pot.matrix.select("residue.aa.from", "residue.aa.to", "value")
            .rename({"value": "TCRen.LOO"})
            .with_columns(pl.lit(pid).alias("pdb.id"))
        )
    return pl.concat(frames)


def derive_tcren_by_type(
    contacts: pl.DataFrame,
    *,
    min_count: int = 30,
    **kwargs,
) -> tuple[dict[str, Potential], pl.DataFrame]:
    """Derive one potential per contact type, plus the occupancy report that says whether to trust it.

    The review's suggestion: a contact potential scores a residue pair by identity alone, so it gives
    the same energy to a Lys–Asp salt bridge and a Lys–Asp pair that merely drifts within 5 Å.
    Conditioning the counts on :mod:`tcren.contact_types` separates them. The review also names the
    risk, and it is the real one: splitting a fixed set of contacts across eight types multiplies the
    sparsity of a 20×20 matrix by eight.

    So this returns the report alongside the potentials rather than only the matrices. Read the
    report first: ``n_contacts`` per type, and ``frac_cells_ge_min`` — the share of the 400 cells that
    reach ``min_count`` observations. A type where that is near zero has a matrix made mostly of
    pseudocount, whatever its numbers look like.

    **Measured on Canonical2026** (8002 typed TCR:peptide contacts, 370 structures): the concern is
    the correct one. No type reaches 5% cell occupancy at ``min_count=30``. ``polar``, the largest
    bucket at 3221 contacts, populates 4.75% of cells with a median of 6.5 observations each;
    ``salt_bridge`` (136 contacts) reaches 11 cells of 400; ``stacking`` 13. Correlation with the
    pooled matrix tracks the count and nothing else — polar 0.57, hydrophobic 0.28, cation_pi 0.03 —
    which is what noise looks like, not distinct chemistry. (The pipeline itself is fine: the pooled
    re-derivation reproduces the shipped ``TCRen_potential.csv`` at r = +0.85.)

    On a set this size, use the type to **filter** contacts instead
    (:func:`tcren.contact_types.type_weights`); this function is here so the decision can be re-taken
    against a larger set rather than argued about.

    Args:
        contacts: a contact table carrying ``contact.type`` — from
            :func:`tcren.paper.helpers.contact_table` with ``contact_types=True``, or
            :func:`tcren.contact_types.residue_pair_types`.
        min_count: observations a cell needs before it counts as populated in the report.
        **kwargs: passed through to :func:`derive_tcren`.

    Returns:
        ``(potentials, report)`` — a ``{contact_type: Potential}`` mapping and a polars frame with
        ``contact.type``, ``n_contacts``, ``n_cells_observed``, ``frac_cells_ge_min``,
        ``median_count``.

    Raises:
        ValueError: if ``contacts`` has no ``contact.type`` column.
    """
    if "contact.type" not in contacts.columns:
        raise ValueError("contacts must carry a 'contact.type' column; build the table with "
                         "contact_table(..., contact_types=True) or residue_pair_types()")

    potentials, rows = {}, []
    for ctype in sorted(contacts["contact.type"].drop_nulls().unique().to_list()):
        sub = contacts.filter(pl.col("contact.type") == ctype)
        cells = sub.group_by("residue.aa.from", "residue.aa.to").agg(pl.len().alias("count"))
        rows.append({
            "contact.type": ctype,
            "n_contacts": sub.height,
            "n_cells_observed": cells.height,
            "frac_cells_ge_min": float((cells["count"] >= min_count).sum()) / 400.0,
            "median_count": float(cells["count"].median()) if cells.height else 0.0,
        })
        potentials[ctype] = derive_tcren(sub, **kwargs)
        potentials[ctype].name = f"TCRen[{ctype}]"

    report = pl.DataFrame(rows, schema={
        "contact.type": pl.Utf8, "n_contacts": pl.Int64, "n_cells_observed": pl.Int64,
        "frac_cells_ge_min": pl.Float64, "median_count": pl.Float64,
    }).sort("n_contacts", descending=True)
    return potentials, report
