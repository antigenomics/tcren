"""Download and version the TCR3D native-structures database.

Fetches the TCR3D complex CIFs and annotation tables into a :class:`NativeDatabase`
root, records provenance (remote ``Last-Modified`` per file, SHA-256, sizes) in
``version.json`` and a per-CIF ``manifest.tsv``, and exposes update checks so callers
can verify the local copy is current before using it.
"""

from __future__ import annotations

import hashlib
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import polars as pl
import requests

from .database import (
    CHAIN_DATA,
    CIF_SUFFIX,
    COMPLEX_DATA,
    NativeDatabase,
)

BASE_URL = "https://tcr3d.ibbr.umd.edu/static/download"
TARBALL = "TCR_complexes_cif.tar.gz"
FILES = (TARBALL, CHAIN_DATA, COMPLEX_DATA)


def _url(name: str) -> str:
    return f"{BASE_URL}/{name}"


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def remote_metadata(timeout: float = 30.0) -> dict[str, dict]:
    """HEAD each remote file and return its ``Last-Modified`` and size."""
    meta = {}
    for name in FILES:
        resp = requests.head(_url(name), allow_redirects=True, timeout=timeout)
        resp.raise_for_status()
        meta[name] = {
            "last_modified": resp.headers.get("Last-Modified"),
            "bytes": int(resp.headers.get("Content-Length", 0)),
        }
    return meta


def needs_update(db: NativeDatabase, timeout: float = 30.0) -> list[str]:
    """Return the files whose remote ``Last-Modified`` differs from the stored version.

    An empty list means the local copy is current. Performs network HEAD requests.
    """
    stored = db.version().get("files", {})
    if not stored:
        return list(FILES)
    remote = remote_metadata(timeout=timeout)
    changed = []
    for name in FILES:
        if remote[name]["last_modified"] != stored.get(name, {}).get("last_modified"):
            changed.append(name)
    return changed


def _download_file(url: str, dest: Path, timeout: float) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=timeout) as resp:
        resp.raise_for_status()
        last_modified = resp.headers.get("Last-Modified")
        with dest.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
    return last_modified


def _write_manifest(db: NativeDatabase) -> None:
    chain_ids = set(db.chain_data["pdb_id"].to_list())
    complex_ids = set(db.complex_data["PDB_ID"].to_list())
    rows = []
    for cif in db.cif_files():
        pdb_id = cif.name[: -len(CIF_SUFFIX)]
        rows.append(
            {
                "pdb_id": pdb_id,
                "filename": cif.name,
                "bytes": cif.stat().st_size,
                "sha256": _sha256(cif),
                "in_chain_data": pdb_id in chain_ids,
                "in_complex_data": pdb_id in complex_ids,
            }
        )
    pl.DataFrame(rows).write_csv(db.manifest_path, separator="\t")


def bootstrap(
    db: NativeDatabase | None = None,
    force: bool = False,
    timeout: float = 300.0,
    now: str | None = None,
) -> NativeDatabase:
    """Download, extract and version the native database.

    Args:
        db: Target database (defaults to the standard root).
        force: Re-download even if files already exist.
        timeout: Per-request timeout in seconds.
        now: Override the recorded download timestamp (ISO 8601); defaults to now.

    Returns:
        The populated :class:`NativeDatabase`.
    """
    db = db or NativeDatabase()
    db.root.mkdir(parents=True, exist_ok=True)

    # Best-effort remote Last-Modified so provenance is recorded even when files already
    # exist locally (e.g. fetched out of band); falls back to the stored value offline.
    try:
        remote = remote_metadata(timeout=timeout)
    except requests.RequestException:
        remote = {}
    stored_files = db.version().get("files", {})

    file_meta: dict[str, dict] = {}
    for name in FILES:
        dest = db.root / name
        if force or not dest.exists():
            last_modified = _download_file(_url(name), dest, timeout)
        else:
            last_modified = remote.get(name, {}).get("last_modified") or stored_files.get(
                name, {}
            ).get("last_modified")
        file_meta[name] = {
            "url": _url(name),
            "last_modified": last_modified,
            "bytes": dest.stat().st_size,
            "sha256": _sha256(dest),
        }

    # (Re)extract the CIF tarball into <root>/cif/.
    tar_path = db.root / TARBALL
    if force or not db.cif_dir.is_dir() or not any(db.cif_dir.glob(f"*{CIF_SUFFIX}")):
        with tarfile.open(tar_path) as tar:
            tar.extractall(db.root, filter="data")

    db.__dict__.pop("chain_data", None)  # drop cached_property values
    db.__dict__.pop("complex_data", None)
    _write_manifest(db)

    version = {
        "source": "tcr3d.ibbr.umd.edu",
        "downloaded_at": now or datetime.now(timezone.utc).isoformat(),
        "files": file_meta,
        "n_cif": len(db.cif_files()),
        "n_complexes": db.complex_data.height,
        "n_chain_rows": db.chain_data.height,
    }
    db.version_path.write_text(json.dumps(version, indent=2))
    return db


def ensure(
    db: NativeDatabase | None = None,
    auto_update: bool = False,
    offline: bool = False,
    timeout: float = 30.0,
) -> NativeDatabase:
    """Ensure the native database is present (and optionally current) before use.

    Args:
        db: Target database (defaults to the standard root).
        auto_update: If True, HEAD the remote files and re-download any that changed.
        offline: Skip all network access; only verify local presence.

    Returns:
        The ready :class:`NativeDatabase`.

    Raises:
        FileNotFoundError: If the database is absent and ``offline`` is set.
    """
    db = db or NativeDatabase()
    if not db.is_present():
        if offline:
            raise FileNotFoundError(
                f"native database missing at {db.root} and offline=True; "
                "run `tcren native bootstrap`"
            )
        return bootstrap(db)
    if auto_update and not offline:
        changed = needs_update(db, timeout=timeout)
        if changed:
            return bootstrap(db, force=True)
    return db
