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
    import tempfile
    from arda._locking import build_lock

    fasta = out_dir / "alleles.aa.fasta"
    # Unique staging prevents concurrent builders from truncating each other's files.
    with tempfile.TemporaryDirectory(prefix=".reference-", dir=out_dir) as tmp:
        stage = Path(tmp)
        with (stage / fasta.name).open("w") as fh:
            for al in alleles:
                fh.write(f">{_header(al)}\n{al.sequence}\n")
        pl.DataFrame({f: [getattr(a, f) for a in alleles] for f in _META_FIELDS}).write_csv(
            stage / "metadata.tsv", separator="\t", quote_style="never")
        with build_lock(out_dir / ".reference.lock", done=lambda: False):
            os.replace(stage / "metadata.tsv", out_dir / "metadata.tsv")
            os.replace(stage / fasta.name, fasta)
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


def _db_complete(db: Path) -> bool:
    """Check MMseqs sequence/header/index data against their offset tables, not just markers."""
    try:
        for suffix in ("", "_h", ".idx"):
            prefix = db.with_name(db.name + suffix)
            if Path(str(prefix) + ".dbtype").stat().st_size != 4:
                return False
            parts = ([prefix] if prefix.is_file() else
                     [p for p in prefix.parent.glob(prefix.name + ".*")
                      if p.name[len(prefix.name) + 1:].isdigit()])
            size = sum(p.stat().st_size for p in parts)
            with Path(str(prefix) + ".index").open() as fh:
                end = max(int(row[1]) + int(row[2]) for row in
                          (line.split() for line in fh))
            if not size or size < end:
                return False
        return True
    except (OSError, ValueError, IndexError):
        return False


def reference_db(cache_dir: Path = CACHE_DIR) -> Path:
    """Return an immutable, content/version-bound pre-indexed MHC database.

    Writers publish a complete generation before replacing its pointer. Existing
    searches retain their old generation, including during repair or reference updates.
    Legacy flat indexes remain untouched because their source identity is unknown.
    """
    import hashlib
    import json
    import shutil
    import tempfile
    import uuid

    import arda.mmseqs as mmseqs
    from arda._locking import build_lock
    from ..provenance import searched_index_identity

    fasta = reference_fasta()
    # Snapshot the small reference input so a concurrent atomic FASTA replacement
    # cannot change the bytes consumed by createdb after the identity was computed.
    source = fasta.read_bytes()
    identity = {"source_sha256": hashlib.sha256(source).hexdigest(),
                "mmseqs_version": mmseqs.version_key(mmseqs.version())}
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    pointer = cache_dir / (key + ".json")
    generations = cache_dir / "generations"

    def current():
        try:
            record = json.loads(pointer.read_text())
            name = record["generation"]
            if not isinstance(name, str) or Path(name).name != name or name in (".", ".."):
                return None
            db = generations / name / "alleles_db"
            if (record["identity"] == identity and _db_complete(db)
                    and record["files_sha256"] == searched_index_identity(db)["files_sha256"]):
                return db
        except (OSError, ValueError, KeyError, TypeError):
            pass
        return None

    ready = current()
    if ready is not None:
        return ready
    cache_dir.mkdir(parents=True, exist_ok=True)
    with build_lock(cache_dir / ".alleles_db.lock", done=lambda: current() is not None) as ours:
        if ours:
            generations.mkdir(exist_ok=True)
            with tempfile.TemporaryDirectory(prefix=".build-", dir=generations) as tmp:
                stage = Path(tmp)
                snapshot = stage / "source.fasta"
                snapshot.write_bytes(source)
                db = stage / "alleles_db"
                mmseqs.createdb(snapshot, db, dbtype=1)
                mmseqs.run(["createindex", str(db), str(stage / "work"),
                            "--search-type", "1", "--threads", "1"])
                if not _db_complete(db):
                    raise RuntimeError("MMseqs produced an incomplete MHC reference database")
                hashes = searched_index_identity(db)["files_sha256"]
                snapshot.unlink()
                if (stage / "work").exists():
                    shutil.rmtree(stage / "work")
                name = key + "-" + uuid.uuid4().hex
                os.replace(stage, generations / name)
                temporary = cache_dir / (".pointer-" + uuid.uuid4().hex)
                try:
                    temporary.write_text(json.dumps({"generation": name, "identity": identity,
                                                     "files_sha256": hashes}) + "\n")
                    os.replace(temporary, pointer)
                finally:
                    temporary.unlink(missing_ok=True)
    ready = current()
    if ready is None:
        raise RuntimeError("MHC reference database failed publication validation")
    return ready


def parse_header(header: str) -> dict[str, str]:
    """Parse a reference FASTA header back into its metadata fields."""
    return dict(zip(_META_FIELDS, header.split("|")))
