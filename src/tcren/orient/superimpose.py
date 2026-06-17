"""Superimpose query structures onto a canonical database by MHC.

Unlike :func:`tcren.orient.run_folder` (which *builds* a canonical set from native complexes
using the per-class derived frame), :func:`superimpose` brings a **new** structure into the
canonical frame defined by an existing database (``data/Canonical2026`` by default).

How it works, per query structure:

1. Chain-type + MHC-annotate the query; read its MHC **class** (MHCI/MHCII) and **species**.
2. Select every database structure with the *same* class and species (from the database's
   ``orient_metadata.json``).
3. Superpose the query's conserved groove Cα onto each selected database structure (sequence
   alignment establishes the residue correspondence, so alleles/numbering differ freely).
4. Average the resulting rigid transforms — translations by mean, rotations by the chordal
   (SVD-orthonormalised) mean — into one consensus placement, and apply it.

Because every database member already sits in the same canonical frame, each superposition
independently yields a canonical placement; averaging over the whole matching ensemble cancels
the per-structure groove variation.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..paths import data_dir
from ..structure.model import Structure
from .align import OrientationResult, _matched_anchors, apply_transform
from .chains import rename_chains

# Annotated database references, keyed by (db_dir, class, species) — built once per process.
_DB_CACHE: dict[tuple[str, str, str], list[Structure]] = {}


def _metadata_path(db_dir: Path) -> Path:
    """``orient_metadata.json`` inside the database dir, else alongside it (repo ``data/``)."""
    for cand in (db_dir / "orient_metadata.json", db_dir.parent / "orient_metadata.json"):
        if cand.exists():
            return cand
    raise FileNotFoundError(f"no orient_metadata.json for canonical database {db_dir}")


def _matching_ids(db_dir: Path, mhc_class: str, species: str) -> list[str]:
    """Database ids whose metadata matches ``mhc_class`` and ``species`` (case-insensitive)."""
    meta = json.loads(_metadata_path(db_dir).read_text())
    sp = species.lower()
    return [r["pdb.id"] for r in meta
            if r.get("status") == "ok" and r.get("mhc.class") == mhc_class
            and str(r.get("species", "")).lower() == sp]


def _canonical_references(db_dir: Path, mhc_class: str, species: str,
                          organism: str) -> list[Structure]:
    """Load + batch-annotate the matching-class/species database structures (cached)."""
    key = (str(db_dir), mhc_class, species.lower())
    if key in _DB_CACHE:
        return _DB_CACHE[key]

    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..mhc import annotate_mhc_batch
    from ..paper.helpers import _batch_annotate
    from ..structure.io import structure_paths

    by_id = {p.name.split(".")[0]: p for p in structure_paths(db_dir)}
    refs: list[Structure] = []
    for rid in _matching_ids(db_dir, mhc_class, species):
        if rid in by_id:
            from ..structure import parse_structure
            refs.append(parse_structure(by_id[rid], pdb_id=rid))
    records = _batch_annotate(refs, _import_arda())
    for s, recs in zip(refs, records):
        classify_chains(s, organism=organism, precomputed_records=recs)
    annotate_mhc_batch(refs)
    _DB_CACHE[key] = refs
    return refs


def _average_transform(transforms: list[tuple[np.ndarray, np.ndarray]]) -> tuple[np.ndarray, np.ndarray]:
    """Chordal-mean rotation (SVD-orthonormalised) + mean translation over ``transforms``."""
    rot_mean = np.mean([R for R, _ in transforms], axis=0)
    u, _s, vt = np.linalg.svd(rot_mean)
    rot = u @ vt
    if np.linalg.det(rot) < 0:  # reflect back to a proper rotation
        u[:, -1] *= -1
        rot = u @ vt
    tran = np.mean([t for _, t in transforms], axis=0)
    return rot, tran


def _query_class_species(structure: Structure) -> tuple[str, str]:
    mhc_class = next((c.chain_supertype for c in structure.chains
                      if c.chain_type in ("MHCa", "MHCb")), "MHCI")
    return mhc_class, structure.complex_species or "Human"


def superimpose(
    structure: Structure,
    db_dir: str | Path | None = None,
    organism: str = "human",
    annotate: bool = True,
) -> tuple[Structure, OrientationResult]:
    """Superimpose ``structure`` onto a canonical database by MHC (see module docstring).

    ``structure`` is chain-typed + MHC-annotated here unless ``annotate=False`` (used by the
    threaded batch driver, which annotates the whole input set in one mmseqs pass first). The
    ensemble alignment itself is mmseqs-free, so it is the part safe to run on a thread pool.
    ``db_dir`` defaults to ``data/Canonical2026``. Returns the oriented, A–E renamed structure
    and the consensus :class:`~tcren.orient.align.OrientationResult` (averaged over the ensemble).
    """
    from Bio.SVDSuperimposer import SVDSuperimposer

    if annotate:
        from ..annotation import classify_chains
        from ..mhc import annotate_mhc

        classify_chains(structure, organism=organism)
        annotate_mhc(structure)
    mhc_class, species = _query_class_species(structure)

    db_dir = Path(db_dir) if db_dir is not None else data_dir() / "Canonical2026"
    refs = _canonical_references(db_dir, mhc_class, species, organism)
    if not refs:
        raise ValueError(f"no {mhc_class}/{species} structures in canonical database {db_dir}")

    transforms, rmsds = [], []
    for ref in refs:
        mob, rp = _matched_anchors(structure, ref)
        if len(mob) < 3:
            continue
        sup = SVDSuperimposer()
        sup.set(rp, mob)
        sup.run()
        transforms.append(sup.get_rotran())
        rmsds.append(float(sup.get_rms()))
    if not transforms:
        raise ValueError(f"too few matched groove anchors to superimpose {structure.pdb_id}")

    rot, tran = _average_transform(transforms)
    result = OrientationResult(
        rotation=rot, translation=tran, rmsd=float(np.mean(rmsds)),
        n_anchor_atoms=len(transforms), reference_id=f"ensemble:{mhc_class}:{species}",
    )
    oriented = apply_transform(structure, result)
    oriented.cell_type = structure.cell_type
    oriented, chain_map = rename_chains(oriented)
    return oriented, result
