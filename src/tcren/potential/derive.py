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

# Classic derivation enumerates the 20 standard amino acids (Cys included in the
# grid; dropped from the *from* axis only after the log-odds are computed).
_AA20_CLASSIC: tuple[str, ...] = (
    "L", "F", "I", "M", "V", "W", "Y", "C", "H", "A",
    "G", "P", "T", "S", "Q", "N", "D", "E", "R", "K",
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
        weights: Optional per-structure weights ``{pdb.id: weight}``. When given, each
            structure's contributions to the aa-pair counts are multiplied by its weight
            (rows whose ``pdb.id`` is absent from the map default to weight ``1.0``);
            this down-weights redundancy while keeping all data (see
            :func:`tcren.potential.redundancy.cluster_weights`). ``None`` (default) is
            unweighted and byte-identical to the legacy derivation.

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

    alphabet = _AA20_CLASSIC if variant == "classic" else AA21
    if drop_cys is None:
        drop_cys = variant == "classic"

    if weights is None:
        # Unweighted: one row = one count (byte-identical to the legacy path).
        n_contacts = df.height
        counts = df.group_by("residue.aa.from", "residue.aa.to").agg(
            pl.len().alias("count")
        )
    else:
        # Weighted: each row contributes its structure's weight (default 1.0).
        w = df["pdb.id"].replace_strict(
            weights, default=1.0, return_dtype=pl.Float64
        )
        df = df.with_columns(w.alias("_w"))
        n_contacts = float(df["_w"].sum())
        counts = df.group_by("residue.aa.from", "residue.aa.to").agg(
            pl.col("_w").sum().alias("count")
        )
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
