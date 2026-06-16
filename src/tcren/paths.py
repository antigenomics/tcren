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


def data_dir() -> Path:
    """Root of the runtime dataset: ``$TCREN_DATA_DIR`` or the repo ``data/`` directory."""
    env = os.environ.get("TCREN_DATA_DIR")
    return Path(env) if env else _REPO / "data"


def native_dir() -> Path:
    """Directory holding the canonical ``Native2026`` structures (``data/Native2026``)."""
    return data_dir() / NATIVE2026


def reference_structure_path(pdb_id: str) -> Path:
    """Resolve a Native2026 structure file by id, trying plain and gzipped PDB/mmCIF.

    Raises ``FileNotFoundError`` (with a bootstrap hint) if no candidate exists — run
    ``tcren paper bootstrap`` (or fetch the HF dataset) to populate ``data/Native2026``.
    """
    base = native_dir()
    for suffix in STRUCTURE_SUFFIXES:
        for name in (f"{pdb_id}{suffix}", f"{pdb_id}{suffix}.gz"):
            cand = base / name
            if cand.exists():
                return cand
    raise FileNotFoundError(
        f"{pdb_id} not found in {base}. Populate Native2026 (e.g. `tcren paper bootstrap`)."
    )
