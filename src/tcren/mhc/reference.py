"""Build and load the curated MHC reference shipped under ``database/mhc/``.

The committed reference is a single FASTA (``alleles.aa.fasta``) whose headers encode the
metadata (``allele|locus|mhc_class|chain_role|species``) plus a ``metadata.tsv`` mirror.
The mmseqs search index is built on demand into a gitignored cache (mirroring arda's
commit-FASTA / build-index-on-demand split).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from . import imgt
from .imgt import MhcAllele

_REPO = Path(__file__).resolve().parents[3]
DATABASE_DIR = _REPO / "database" / "mhc"
CACHE_DIR = _REPO / "data" / "mhc_cache"

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
    fasta = out_dir / "alleles.aa.fasta"
    with fasta.open("w") as fh:
        for al in alleles:
            fh.write(f">{_header(al)}\n{al.sequence}\n")

    pl.DataFrame(
        {f: [getattr(a, f) for a in alleles] for f in _META_FIELDS}
    ).write_csv(out_dir / "metadata.tsv", separator="\t")
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
    """
    import tempfile

    import arda.mmseqs as mmseqs

    fasta = reference_fasta()
    db = cache_dir / "alleles_db"
    db_marker = db.with_name(db.name + ".dbtype")        # createdb output
    idx_marker = db.with_name(db.name + ".idx.dbtype")   # createindex output
    fasta_mtime = fasta.stat().st_mtime
    stale = (not db_marker.exists() or db_marker.stat().st_mtime < fasta_mtime
             or not idx_marker.exists() or idx_marker.stat().st_mtime < fasta_mtime)
    if stale:
        cache_dir.mkdir(parents=True, exist_ok=True)
        mmseqs.createdb(fasta, db, dbtype=1)
        with tempfile.TemporaryDirectory() as tmp:
            mmseqs.run(["createindex", str(db), tmp, "--search-type", "1"])
    return db


def parse_header(header: str) -> dict[str, str]:
    """Parse a reference FASTA header back into its metadata fields."""
    return dict(zip(_META_FIELDS, header.split("|")))
