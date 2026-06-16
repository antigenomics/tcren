"""Access to a local TCR3D native-structures database.

A :class:`NativeDatabase` points at a directory holding the TCR3D download:

* ``cif/<pdb_id>_renumbered.cif`` — renumbered complex structures,
* ``tcr_chain_data.tsv`` — per-chain TCR annotation (V/J genes, CDR sequences),
* ``tcr_complexes_data.tsv`` — per-complex annotation (MHC allele, epitope, geometry),
* ``version.json`` / ``manifest.tsv`` — provenance written by the bootstrap step.

The default root is ``$TCREN_NATIVE_DIR`` or ``<repo>/data/native``. Point a database at
any other directory (e.g. a previous TCR3D release you have on disk) to use a custom set.
"""

from __future__ import annotations

import json
import os
from functools import cached_property
from pathlib import Path

import polars as pl

_REPO = Path(__file__).resolve().parents[3]
CIF_SUFFIX = "_renumbered.cif"
CHAIN_DATA = "tcr_chain_data.tsv"
COMPLEX_DATA = "tcr_complexes_data.tsv"
VERSION_FILE = "version.json"
MANIFEST_FILE = "manifest.tsv"


def default_native_root() -> Path:
    """Resolve the default native-database root (``$TCREN_NATIVE_DIR`` or repo data dir)."""
    env = os.environ.get("TCREN_NATIVE_DIR")
    return Path(env) if env else _REPO / "data" / "native"


class NativeDatabase:
    """A view over a local TCR3D native-structures directory."""

    def __init__(self, root: str | Path | None = None):
        self.root = Path(root) if root is not None else default_native_root()

    # --- paths ------------------------------------------------------------
    @property
    def cif_dir(self) -> Path:
        return self.root / "cif"

    @property
    def chain_data_path(self) -> Path:
        return self.root / CHAIN_DATA

    @property
    def complex_data_path(self) -> Path:
        return self.root / COMPLEX_DATA

    @property
    def version_path(self) -> Path:
        return self.root / VERSION_FILE

    @property
    def manifest_path(self) -> Path:
        return self.root / MANIFEST_FILE

    # --- presence / provenance -------------------------------------------
    def is_present(self) -> bool:
        """True if the core files (cif dir + both tables) exist."""
        return (
            self.cif_dir.is_dir()
            and self.chain_data_path.exists()
            and self.complex_data_path.exists()
        )

    def version(self) -> dict:
        """Return the stored version metadata (``{}`` if not bootstrapped)."""
        if self.version_path.exists():
            return json.loads(self.version_path.read_text())
        return {}

    def manifest(self) -> pl.DataFrame:
        """Return the per-CIF manifest (raises if absent)."""
        return pl.read_csv(self.manifest_path, separator="\t")

    # --- data -------------------------------------------------------------
    @cached_property
    def chain_data(self) -> pl.DataFrame:
        return pl.read_csv(self.chain_data_path, separator="\t")

    @cached_property
    def complex_data(self) -> pl.DataFrame:
        return pl.read_csv(self.complex_data_path, separator="\t")

    def cif_files(self) -> list[Path]:
        """All ``*_renumbered.cif`` files, sorted by pdb id."""
        return sorted(self.cif_dir.glob(f"*{CIF_SUFFIX}"))

    def pdb_ids(self) -> list[str]:
        """PDB ids available as CIF files."""
        return [p.name[: -len(CIF_SUFFIX)] for p in self.cif_files()]

    def cif_for(self, pdb_id: str) -> Path:
        """Path to a single complex CIF (raises ``FileNotFoundError`` if absent)."""
        path = self.cif_dir / f"{pdb_id}{CIF_SUFFIX}"
        if not path.exists():
            raise FileNotFoundError(f"{pdb_id} not in native database at {self.cif_dir}")
        return path
