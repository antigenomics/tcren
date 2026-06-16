"""Re-derive TCRen potentials from the native (TCR3D) structures.

Runs the tcren pipeline (parse → chain typing → contacts) over the native CIFs to build
TCR–peptide contact maps, then derives a statistical potential from them. The contact
table can be cached so re-derivation (e.g. with a different variant or pseudocount) is
cheap.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from ..annotation import classify_chains
from ..contactmap import ContactMap
from ..potential import Potential, derive_tcren
from ..structure import parse_structure
from .database import NativeDatabase

_TCR_PEPTIDE_COLUMNS = [
    "pdb.id",
    "chain.type.from",
    "region.type.from",
    "residue.index.from",
    "residue.index.to",
    "residue.aa.from",
    "residue.aa.to",
]


def _organism_map(db: NativeDatabase) -> dict[str, str]:
    return {
        r["PDB_ID"]: ("mouse" if r["TCR_organism"] == "Mouse" else "human")
        for r in db.complex_data.iter_rows(named=True)
    }


def native_contact_table(
    db: NativeDatabase,
    pdb_ids: list[str] | None = None,
    cutoff: float = 5.0,
    on_error: str = "skip",
) -> pl.DataFrame:
    """Compute TCR↔peptide contacts across native structures.

    Args:
        db: The native database.
        pdb_ids: Subset of pdb ids (defaults to all CIFs present).
        cutoff: Contact distance cutoff (Å).
        on_error: ``"skip"`` to drop structures that fail to process, ``"raise"``.

    Returns:
        Stacked TCR↔peptide contact rows (one table over all structures).
    """
    organisms = _organism_map(db)
    ids = pdb_ids if pdb_ids is not None else db.pdb_ids()
    frames = []
    for pdb_id in ids:
        try:
            s = parse_structure(db.cif_for(pdb_id), pdb_id=pdb_id)
            classify_chains(s, organism=organisms.get(pdb_id, "human"))
            cm = ContactMap.from_structure(s, cutoff=cutoff)
            tp = cm.tcr_peptide()
            if tp.height:
                frames.append(tp.select(_TCR_PEPTIDE_COLUMNS))
        except Exception:
            if on_error == "raise":
                raise
    return pl.concat(frames) if frames else pl.DataFrame(schema={c: pl.Utf8 for c in _TCR_PEPTIDE_COLUMNS})


def cache_path(db: NativeDatabase) -> Path:
    return db.root / "derived" / "native_contact_maps.parquet"


def precompute_contacts(db: NativeDatabase, cutoff: float = 5.0) -> Path:
    """Compute and cache the native contact table under ``<root>/derived/``."""
    table = native_contact_table(db, cutoff=cutoff)
    out = cache_path(db)
    out.parent.mkdir(parents=True, exist_ok=True)
    table.write_parquet(out)
    return out


def derive_native_potential(
    db: NativeDatabase,
    variant: str = "classic",
    pseudocount: int = 1,
    include: list[str] | None = None,
    use_cache: bool = True,
) -> Potential:
    """Derive a TCRen potential from the native structures.

    Uses the cached contact table when present (and ``use_cache``), otherwise computes
    it on the fly.
    """
    cache = cache_path(db)
    if use_cache and cache.exists():
        contacts = pl.read_parquet(cache)
    else:
        contacts = native_contact_table(db)
    return derive_tcren(contacts, include=include, variant=variant, pseudocount=pseudocount)
