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

# ContactMap.tcr_peptide() column -> R analysis column name.
_CONTACT_RENAME = {
    "pdb.id": "pdb.id",
    "chain.type.from": "chain.type.from",
    "region.type.from": "region.type.from",
    "residue.index.from": "residue.index.from",
    "residue.index.to": "residue.index.to",
    "pos.from": "pos.from",
    "pos.to": "pos.to",
    "residue.aa.from": "residue.aa.from",
    "residue.aa.to": "residue.aa.to",
}


def contact_table(structure: Structure, cutoff: float = 5.0) -> pl.DataFrame:
    """TCR↔peptide contact table for an annotated structure (the mir-replacement).

    The structure must already be chain-typed (``classify_chains``) and MHC-annotated
    (``annotate_mhc``). Returns the columns the R benchmarks use:
    ``pdb.id, chain.type.from, region.type.from, residue.index.from, residue.index.to,
    pos.from, pos.to, residue.aa.from, residue.aa.to``.
    """
    tp = ContactMap.from_structure(structure, cutoff=cutoff).tcr_peptide()
    return tp.select(list(_CONTACT_RENAME)).unique()


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
