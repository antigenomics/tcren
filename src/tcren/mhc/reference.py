"""Build and load the curated MHC reference under ``database/mhc/``.

The reference is a single FASTA (``alleles.aa.fasta``) whose headers encode the
metadata (``allele|locus|mhc_class|chain_role|species``) plus a ``metadata.tsv`` mirror. It is
built on demand from IMGT by ``tcren build-mhc-ref`` and written under
:func:`tcren.paths.tcren_home`, not bundled into the wheel.
The mmseqs search index is built on demand into a gitignored cache (mirroring arda's
commit-FASTA / build-index-on-demand split).
"""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from ..paths import tcren_home
from . import imgt
from .imgt import MhcAllele

DATABASE_DIR = tcren_home() / "database" / "mhc"
CACHE_DIR = tcren_home() / "data" / "mhc_cache"

_META_FIELDS = ("allele", "locus", "mhc_class", "chain_role", "species")


def _header(allele: MhcAllele) -> str:
    return "|".join(
        (allele.allele, allele.locus, allele.mhc_class, allele.chain_role, allele.species)
    )


def build(
    species: tuple[str, ...] = ("human", "mouse"),
    cache_dir: Path = CACHE_DIR,
    out_dir: Path = DATABASE_DIR,
    force_download: bool = False,
) -> Path:
    """Download, curate and write the committed MHC reference.

    Args:
        species: Which species to include.
        cache_dir: Where raw downloads are cached (gitignored).
        out_dir: Where the curated ``alleles.aa.fasta`` + ``metadata.tsv`` are written.
        force_download: Re-download even if cached files exist.

    Returns:
        Path to the written ``alleles.aa.fasta``.
    """
    alleles: list[MhcAllele] = []
    if "human" in species:
        hla = imgt.download_human(cache_dir, force=force_download)
        alleles += imgt.parse_human(hla)
    if "mouse" in species:
        mouse, human_b2m = imgt.download_mouse(cache_dir, force=force_download)
        alleles += imgt.parse_mouse(mouse, human_b2m)

    out_dir.mkdir(parents=True, exist_ok=True)
    # Write via a temp file + os.replace so a concurrent reader never sees a half-written
    # reference (its presence is the "already built" gate in reference_fasta()).
    fasta = out_dir / "alleles.aa.fasta"
    tmp_fasta = fasta.with_name(fasta.name + ".tmp")
    with tmp_fasta.open("w") as fh:
        for al in alleles:
            fh.write(f">{_header(al)}\n{al.sequence}\n")
    os.replace(tmp_fasta, fasta)

    meta = out_dir / "metadata.tsv"
    tmp_meta = meta.with_name(meta.name + ".tmp")
    pl.DataFrame(
        {f: [getattr(a, f) for a in alleles] for f in _META_FIELDS}
    ).write_csv(tmp_meta, separator="\t")
    os.replace(tmp_meta, meta)
    return fasta


def reference_fasta(out_dir: Path = DATABASE_DIR) -> Path:
    """Path to the committed reference FASTA (raise if the reference is not built)."""
    fasta = out_dir / "alleles.aa.fasta"
    if not fasta.exists():
        raise FileNotFoundError(
            f"MHC reference not found at {fasta}; run `tcren.mhc.reference.build()` "
            "or `tcren build-mhc-ref`"
        )
    return fasta


def reference_db(cache_dir: Path = CACHE_DIR) -> Path:
    """Path to a compiled, **pre-indexed** mmseqs DB of the allele reference (built once, cached).

    `mmseqs easy-search` otherwise rebuilds the target DB *and* its k-mer prefilter index from the
    ~28k-allele FASTA on every call. Caching `createdb` saves little; the dominant cost is the
    prefilter index, so we also run `createindex` once. Reusing this DB cuts a single-structure MHC
    search from ~4.5 s to ~0.9 s. Built into the gitignored ``data/mhc_cache`` when missing or older
    than the FASTA.

    The build is serialized through :func:`arda._locking.build_lock`: tcren is routinely run
    concurrently against the same cache (one process per SLURM-array task / Nextflow sample), and
    an unguarded ``createdb`` + ``createindex`` into the shared path would let every other process
    search a half-written index. The ``createindex`` marker is written last, so its freshness gates
    completeness.
    """
    import tempfile

    import arda.mmseqs as mmseqs
    from arda._locking import build_lock

    fasta = reference_fasta()
    db = cache_dir / "alleles_db"
    db_marker = db.with_name(db.name + ".dbtype")        # createdb output
    idx_marker = db.with_name(db.name + ".idx.dbtype")   # createindex output (written last)

    def _fresh() -> bool:
        fasta_mtime = fasta.stat().st_mtime
        return (db_marker.exists() and db_marker.stat().st_mtime >= fasta_mtime
                and idx_marker.exists() and idx_marker.stat().st_mtime >= fasta_mtime)

    if _fresh():
        return db
    cache_dir.mkdir(parents=True, exist_ok=True)
    with build_lock(cache_dir / ".alleles_db.lock", done=_fresh) as ours:
        if ours:
            mmseqs.createdb(fasta, db, dbtype=1)
            with tempfile.TemporaryDirectory() as tmp:
                mmseqs.run(["createindex", str(db), tmp, "--search-type", "1"])
    return db


def parse_header(header: str) -> dict[str, str]:
    """Parse a reference FASTA header back into its metadata fields."""
    return dict(zip(_META_FIELDS, header.split("|")))
