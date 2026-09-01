"""Helpers for the Nat Comput Sci 2022 reproduction notebooks.

``contact_table`` replaces the legacy mir ``extract_contact_map`` (it returns the same
TCR↔peptide contact columns the R analyses consume, computed through the tcren pipeline).
``compare`` is the small regression utility behind ``07_compare_legacy.ipynb``.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from ..contactmap import ContactMap
from ..structure.model import Structure

# The ContactMap.tcr_peptide() columns the R benchmarks consume, in order.
_CONTACT_COLS = [
    "pdb.id", "chain.type.from", "region.type.from",
    "residue.index.from", "residue.index.to",
    "pos.from", "pos.to", "residue.aa.from", "residue.aa.to",
]


def contact_table(
    structure: Structure, cutoff: float = 5.0, count_atoms: bool = False,
    contact_types: bool = False,
) -> pl.DataFrame:
    """TCR↔peptide contact table for an annotated structure (the mir-replacement).

    The structure must already be chain-typed (``classify_chains``) and MHC-annotated
    (``annotate_mhc``). Returns the columns the R benchmarks use:
    ``pdb.id, chain.type.from, region.type.from, residue.index.from, residue.index.to,
    pos.from, pos.to, residue.aa.from, residue.aa.to``.

    When ``count_atoms`` is set, an extra ``n_atom_contacts`` column (the heavy-atom-pair
    count per residue pair) is carried through for atomic-weighted scoring. Default
    ``False`` keeps the schema byte-identical to the legacy output.

    ``contact_types`` adds ``contact.type`` from :func:`tcren.contact_types.residue_pair_types`.
    Without it a cached contact table cannot be typed after the fact — ``atom.from``, ``atom.to``
    and ``dist`` are all dropped here — which is what blocked a type-aware potential derivation.
    """
    tp = ContactMap.from_structure(
        structure, cutoff=cutoff, count_atoms=count_atoms
    ).tcr_peptide()
    cols = list(_CONTACT_COLS)
    if count_atoms:
        cols.append("n_atom_contacts")
    out = tp.select(cols).unique()
    if contact_types:
        from ..contact_types import residue_pair_types
        typed = residue_pair_types(structure, "tcr_peptide", cutoff=cutoff).select(
            "chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to",
            "contact.type")
        keys = ["residue.index.from", "residue.index.to"]
        out = out.join(typed.select(keys + ["contact.type"]).unique(subset=keys),
                       on=keys, how="left")
    return out












def _read_any(path: str | Path) -> pl.DataFrame:
    """Read a CSV/TSV, transparently handling ``.gz`` and tab vs comma."""
    path = Path(path)
    name = path.name[:-3] if path.suffix == ".gz" else path.name
    sep = "\t" if name.endswith((".tsv", ".txt")) else ","
    return pl.read_csv(path, separator=sep, infer_schema_length=2000)


def compare(
    old_path: str | Path,
    new_path: str | Path,
    keys: list[str],
    value_cols: list[str] | None = None,
    tol: float = 1e-6,
) -> dict:
    """Compare two tables on ``keys`` and report row-set + max numeric differences.

    Returns ``{rows_old, rows_new, matched, only_old, only_new, max_abs_diff, status}``
    where ``status`` is ``"pass"`` when the key sets agree and every shared numeric column
    differs by ≤ ``tol``.
    """
    old, new = _read_any(old_path), _read_any(new_path)
    ko = set(map(tuple, old.select(keys).rows()))
    kn = set(map(tuple, new.select(keys).rows()))
    only_old, only_new = ko - kn, kn - ko

    max_abs = 0.0
    if value_cols is None:
        value_cols = [
            c for c in old.columns
            if c in new.columns and c not in keys and old[c].dtype.is_numeric()
        ]
    if value_cols and not only_old and not only_new:
        joined = old.join(new, on=keys, how="inner", suffix="__new")
        for c in value_cols:
            diff = (joined[c] - joined[f"{c}__new"]).abs().max()
            if diff is not None:
                max_abs = max(max_abs, float(diff))

    status = "pass" if not only_old and not only_new and max_abs <= tol else "FAIL"
    return {
        "rows_old": old.height, "rows_new": new.height,
        "matched": len(ko & kn), "only_old": len(only_old), "only_new": len(only_new),
        "max_abs_diff": max_abs, "status": status,
    }


#: Backwards-compatible alias for the pre-2.3.0 private name.


# --- moved to `tcren.annotation.batch` -----------------------------------------------------------
# The batched annotation helpers are infrastructure, not paper code, and five modules below this one
# in the stack were importing them from here. Re-exported so nothing written against the old
# location breaks.
from ..annotation.batch import (  # noqa: E402,F401
    _annotate,
    _batch_annotate,
    annotate_batch,
    annotate_structure_set,
    iter_annotated_set,
    iter_typed,
    mhc_annotation,
)
