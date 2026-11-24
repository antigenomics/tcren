"""Helpers for the Nat Comput Sci 2022 reproduction notebooks.

``contact_table`` replaces the legacy mir ``extract_contact_map`` (it returns the same
TCR↔peptide contact columns the R analyses consume, computed through the tcren pipeline).
``compare`` is the small regression utility behind ``07_compare_legacy.ipynb``.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from ..contactmap import contact_table as contact_table

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
