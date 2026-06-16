"""Filesystem locations for tcren's reference data.

The library's runtime dataset lives in the repo ``data/`` directory (or ``$TCREN_DATA_DIR``):
the canonical ``Native2026`` structure set (HF ``isalgo/tcren_structures``, gitignored),
``PDB_date.tsv`` and ``orient_metadata.json``. Structures are fetched lazily; nothing here is
bundled into the installed package.
"""

from __future__ import annotations

import os
from pathlib import Path

from .structure.io import STRUCTURE_SUFFIXES

_REPO = Path(__file__).resolve().parents[2]
NATIVE2026 = "Native2026"
# The canonical reference structures (and full Native2026 set) live in this HF dataset.
HF_REPO = "isalgo/tcren_structures"


def data_dir() -> Path:
    """Root of the runtime dataset: ``$TCREN_DATA_DIR`` or the repo ``data/`` directory."""
    env = os.environ.get("TCREN_DATA_DIR")
    return Path(env) if env else _REPO / "data"


def native_dir() -> Path:
    """Directory holding the canonical ``Native2026`` structures (``data/Native2026``)."""
    return data_dir() / NATIVE2026


def _local_reference(pdb_id: str) -> Path | None:
    """A Native2026 structure file for ``pdb_id`` under the data dir, if present locally."""
    base = native_dir()
    for suffix in STRUCTURE_SUFFIXES:
        for name in (f"{pdb_id}{suffix}", f"{pdb_id}{suffix}.gz"):
            cand = base / name
            if cand.exists():
                return cand
    return None


def _fetch_reference_from_hf(pdb_id: str, folder: str = NATIVE2026) -> Path | None:
    """Download (and cache) a single reference structure from the HF dataset.

    Returns the cached file path, or ``None`` if ``huggingface_hub`` is missing or the file
    cannot be fetched. ``hf_hub_download`` caches under the HF cache, so repeat lookups are
    local (no network). This is what lets an installed library/CLI orient a brand-new,
    non-canonical structure without a populated repo ``data/``.
    """
    try:
        from huggingface_hub import hf_hub_download  # noqa: PLC0415
    except ImportError:
        return None
    for suffix in (".pdb.gz", ".cif.gz", ".pdb", ".cif"):
        try:
            path = hf_hub_download(HF_REPO, f"{folder}/{pdb_id}{suffix}", repo_type="dataset")
            return Path(path)
        except Exception:  # noqa: BLE001 - try the next suffix / fall through to None
            continue
    return None


def reference_structure_path(pdb_id: str) -> Path:
    """Resolve a canonical reference structure by id (plain/gzipped PDB/mmCIF).

    Looks under ``data/Native2026`` first; if absent (e.g. a pip-installed library with no
    repo ``data/``), lazily downloads it from the HF dataset into the HF cache. This makes
    orienting a new, non-canonical structure work out of the box for both the library and CLI.

    Raises ``FileNotFoundError`` if it is neither local nor fetchable.
    """
    local = _local_reference(pdb_id)
    if local is not None:
        return local
    fetched = _fetch_reference_from_hf(pdb_id)
    if fetched is not None:
        return fetched
    raise FileNotFoundError(
        f"{pdb_id} not found in {native_dir()} and could not be fetched from {HF_REPO}. "
        f"Populate Native2026 (`tcren paper bootstrap`) or install `huggingface_hub`."
    )
