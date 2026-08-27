"""Per-structure metadata that travels **with** a structure set.

A structure set is a directory (or archive) of TCR-pMHC complexes. Everything a descriptor cannot
be computed from — the binding label, the epitope and allele, and above all the **generator's own
confidence** (``iptm``, ``plddt``, ``ranking_confidence``) — lives beside the structures in a
``metadata.tsv``, keyed by the same id ``tcren features`` writes into ``complex.id``.

The rule this enforces: *a set that ships models ships their confidences*. Without it every
analysis has to rediscover where the confidences went, and they went somewhere different for each
set — which is how a benchmark ends up joined on a hash that is not unique.

Layout, beside the structures::

    <set>/
        metadata.tsv        id + whatever is known
        1ao7.pdb.gz
        ...

``id`` matches :func:`tcren.structure.structure_stem` — the file stem with structure suffixes
removed — so a table produced by ``tcren features`` joins on ``complex.id`` with no massaging.

Reserved column names, all optional except ``id``:

==========================  ================================================================
``id``                      structure stem; the join key. **Required.**
``y``                       binding label, 1/0, when the set carries one
``epitope``                 peptide sequence
``mhc``                     allele, e.g. ``HLA-A*02:01``
``iptm``                    AlphaFold/TCRmodel2 interface pTM
``plddt``                   mean pLDDT
``ptm``                     pTM
``ranking_confidence``      the generator's own ranking scalar
``provenance``              free text: what produced this structure
==========================  ================================================================

Any further column is passed through untouched.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

__all__ = ["METADATA_FILE", "RESERVED", "find_metadata", "read_metadata", "join_metadata"]

METADATA_FILE = "metadata.tsv"
RESERVED = ("id", "y", "epitope", "mhc", "iptm", "plddt", "ptm", "ranking_confidence",
            "provenance")


def find_metadata(structures: str | Path) -> Path | None:
    """The ``metadata.tsv`` governing ``structures``, or ``None``.

    ``structures`` may be the set directory, a file inside it, or an archive; the search walks up
    from the path so that pointing at one structure still finds its set's table.
    """
    p = Path(structures)
    for d in ([p] if p.is_dir() else []) + list(p.parents)[:3]:
        f = d / METADATA_FILE
        if f.is_file():
            return f
    return None


def read_metadata(structures: str | Path) -> pl.DataFrame | None:
    """Read the set's metadata table, or ``None`` if it has none.

    Raises:
        ValueError: if the file exists but has no ``id`` column, which would make it unjoinable.
    """
    f = find_metadata(structures)
    if f is None:
        return None
    t = pl.read_csv(f, separator="\t", infer_schema_length=None)
    if "id" not in t.columns:
        raise ValueError(f"{f} has no 'id' column; it cannot be joined to a feature table")
    if t["id"].n_unique() != t.height:
        raise ValueError(f"{f} has duplicate ids; the join key must be unique")
    return t


def join_metadata(table: pl.DataFrame, structures: str | Path, *, on: str = "complex.id",
                  columns: tuple[str, ...] | None = None) -> pl.DataFrame:
    """Left-join a set's metadata onto a feature table.

    Args:
        table: a ``tcren features`` output, or anything keyed by ``on``.
        structures: the structure set, so its ``metadata.tsv`` can be found.
        on: the key column in ``table``. Default ``complex.id``.
        columns: restrict to these metadata columns (``id`` is always kept). ``None`` takes all.

    Returns:
        ``table`` with the metadata columns appended. Unchanged if the set has no metadata, so this
        is safe to call unconditionally.
    """
    meta = read_metadata(structures)
    if meta is None:
        return table
    if columns is not None:
        meta = meta.select(["id", *(c for c in columns if c in meta.columns and c != "id")])
    clash = (set(meta.columns) & set(table.columns)) - {"id"}
    if clash:
        meta = meta.rename({c: f"meta.{c}" for c in clash})
    return table.join(meta.rename({"id": on}), on=on, how="left")
