"""Filesystem locations for tcren's reference data.

The library's runtime dataset lives under :func:`tcren_home` --- the source checkout when
tcren is run from one, ``$TCREN_HOME`` when set, and a user cache directory otherwise. It holds
the canonical ``Native2026`` structure set (HF ``isalgo/tcren_structures``, gitignored),
``PDB_date.tsv`` and the built MHC allele reference. Structures are fetched lazily; nothing here
is bundled into the installed package, except ``Canonical2026``'s ``orient_metadata.json``, which
rides in ``tcren/data/`` so an installed ``superimpose`` can describe the database it fetched.
"""

from __future__ import annotations

import os
from pathlib import Path

from .structure.io import STRUCTURE_SUFFIXES

NATIVE2026 = "Native2026"
# The canonical reference structures (and full Native2026 set) live in this HF dataset.
HF_REPO = "isalgo/tcren_structures"


def tcren_home() -> Path:
    """Root of tcren's on-disk reference data.

    ``$TCREN_HOME`` when set; otherwise the source checkout, recognised by its
    ``pyproject.toml``. An installed wheel has no checkout above it --- ``parents[2]`` is then
    ``site-packages``' parent --- so it falls back to a user cache directory, which is both
    writable and stable across upgrades.
    """
    env = os.environ.get("TCREN_HOME")
    if env:
        return Path(env)
    checkout = Path(__file__).resolve().parents[2]
    if (checkout / "pyproject.toml").exists():
        return checkout
    cache = os.environ.get("XDG_CACHE_HOME")
    return (Path(cache) if cache else Path.home() / ".cache") / "tcren"


def data_dir() -> Path:
    """Root of the runtime dataset: ``$TCREN_DATA_DIR`` or ``data/`` under :func:`tcren_home`."""
    env = os.environ.get("TCREN_DATA_DIR")
    return Path(env) if env else tcren_home() / "data"


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
