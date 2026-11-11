"""Explicit metadata-to-structure identities, independent of assay outcomes."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ..structure import is_structure_file


@dataclass(frozen=True)
class MetadataInput:
    identity: str
    path: Path
    member: str | None
    row: int


def metadata_inputs(path, *, path_column="structure_path", id_column="id", member_column=None):
    """Validate all identities/paths before writing; stream selected archive members later.

    Structure paths resolve relative to the metadata file. For an archive, an
    explicitly named member column identifies the structure inside it. Extra
    metadata columns, including experimental outcomes, do not affect computation.
    """
    path = Path(path)
    try:
        if path.suffix == ".parquet":
            table = pl.read_parquet(path)
        elif str(path).endswith((".tsv", ".tsv.gz", ".csv", ".csv.gz")):
            table = pl.read_csv(path, separator="," if ".csv" in path.suffixes else "\t",
                                infer_schema=False)
        else:
            raise ValueError("metadata must be TSV[.gz], CSV[.gz] or Parquet")
    except pl.exceptions.PolarsError as exc:
        raise ValueError(f"cannot read metadata {path}: {exc}") from exc
    required = [id_column, path_column] + ([member_column] if member_column else [])
    if len(required) != len(set(required)):
        raise ValueError("identity, path and member columns must be distinct")
    missing = set(required) - set(table.columns)
    if missing:
        raise ValueError(f"missing metadata columns: {', '.join(sorted(missing))}")
    if not table.height:
        raise ValueError("metadata contains no structures")
    seen, output = set(), []
    for ordinal, record in enumerate(table.select(required).iter_rows(named=True), 1):
        identity, filename = record[id_column], record[path_column]
        if identity is None or not str(identity).strip() or str(identity) != str(identity).strip():
            raise ValueError(f"metadata row {ordinal}: identity must be nonempty without outer whitespace")
        identity = str(identity)
        if identity in seen:
            raise ValueError(f"duplicate metadata identity: {identity}")
        seen.add(identity)
        if filename is None or not str(filename).strip():
            raise ValueError(f"metadata row {ordinal}: empty structure path")
        target = Path(str(filename))
        if not target.is_absolute():
            target = path.parent / target
        target = target.resolve()
        if not target.is_file():
            raise ValueError(f"metadata row {ordinal}: structure file missing: {target}")
        member = record.get(member_column) if member_column else None
        member = str(member) if member is not None and str(member) else None
        if member:
            if not str(target).endswith((".tar", ".tar.gz", ".tgz")):
                raise ValueError(f"metadata row {ordinal}: member requires a tar archive")
            if not is_structure_file(member) or Path(member).name.startswith("._"):
                raise ValueError(f"metadata row {ordinal}: invalid structure member: {member}")
        elif not is_structure_file(target) or target.name.startswith("._"):
            raise ValueError(f"metadata row {ordinal}: archive requires an explicit --member-column, or path is not a structure")
        output.append(MetadataInput(identity, target, member, ordinal))
    return output
