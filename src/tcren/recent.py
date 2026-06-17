"""Fetch recent TCR-pMHC structures from the RCSB PDB into ``data/pdb_recent`` (gitignored).

Two entry points, both gzipping validated mmCIF into the destination:

* :func:`fetch_ids` — download specific PDB ids (e.g. the Native2026 set) straight from RCSB.
* :func:`discover_similar` — RCSB full-text search for TCR:peptide:MHC entries (optionally
  released after a date / excluding ids we already have), to surface *new* structures.

Robustness notes baked in per the PDB's current state:
* IDs may be **longer than 4 characters** (extended ``pdb_0000XXXX`` accessions) — handled.
* The PDB is **deprecating split ``.pdb`` files** for large structures, so we always pull
  **mmCIF** (``.cif.gz``), which tcren reads natively.

Every kept structure is annotated (batched, one mmseqs pass) and must have all **5 required
chains** — MHCα, β2m *or* MHCβ, peptide, and the two TCR chains (TRA/TRB or TRG/TRD) — else it
is dropped. ``huggingface_hub``/network are not involved; this uses ``requests`` + the RCSB APIs.
"""

from __future__ import annotations

from pathlib import Path

import requests

from .paths import data_dir

RCSB_FILE = "https://files.rcsb.org/download/{pid}.cif.gz"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"

# The 5 required roles for a complete TCR-pMHC complex (β2m or MHCβ satisfies the MHC second chain).
_TCR = ("TRA", "TRB", "TRG", "TRD")


def recent_dir() -> Path:
    """The gitignored destination for fetched structures (``data/pdb_recent``)."""
    return data_dir() / "pdb_recent"


def _download_cif_gz(pid: str, dest: Path, timeout: float = 60.0, force: bool = False) -> Path | None:
    """Download ``{pid}.cif.gz`` from RCSB into ``dest`` (works for 4-char and extended ids)."""
    out = dest / f"{pid.lower()}.cif.gz"
    if out.exists() and not force:
        return out
    try:
        r = requests.get(RCSB_FILE.format(pid=pid), timeout=timeout)
        r.raise_for_status()
    except requests.RequestException:
        return None
    dest.mkdir(parents=True, exist_ok=True)
    out.write_bytes(r.content)  # RCSB already serves gzip-compressed mmCIF
    return out


def _has_required_chains(structure) -> bool:
    """True if the annotated structure has MHCα + (β2m|MHCβ) + peptide + a TCR pair."""
    types = {c.chain_type for c in structure.chains}
    has_mhc2 = ("B2M" in types) or ("MHCb" in types)
    has_tcr_pair = ({"TRA", "TRB"} <= types) or ({"TRG", "TRD"} <= types)
    return "MHCa" in types and has_mhc2 and "PEPTIDE" in types and has_tcr_pair


def _validate_complete(paths: list[Path], organism: str = "human") -> dict[str, bool]:
    """Batch-annotate downloaded structures; return ``{path_stem: is_complete}``.

    Uses one batched arda pass + one batched MHC pass (no per-structure mmseqs).
    """
    from .annotation import classify_chains
    from .annotation.arda_adapter import _import_arda
    from .mhc import annotate_mhc_batch
    from .paper.helpers import _batch_annotate
    from .structure import import_structure, structure_id_from_path

    structs = []
    for p in paths:
        try:
            structs.append((p, import_structure(p, pdb_id=structure_id_from_path(p))))
        except Exception:  # noqa: BLE001 - unparseable download
            structs.append((p, None))
    ok = [s for _, s in structs if s is not None]
    records = _batch_annotate(ok, _import_arda())
    for s, recs in zip(ok, records):
        classify_chains(s, organism=organism, precomputed_records=recs)
    annotate_mhc_batch(ok)
    return {p.name: (s is not None and _has_required_chains(s)) for p, s in structs}


def discover_similar(after_date: str | None = None, limit: int = 200, timeout: float = 30.0) -> list[str]:
    """RCSB full-text search for TCR:peptide:MHC structures; return candidate PDB ids.

    ``after_date`` (``YYYY-MM-DD``) restricts to entries released on/after it (find *new*
    structures). This is the keyword-driven discovery step; the returned ids are downloaded
    and then strictly validated by :func:`_has_required_chains` (an agent can further curate
    the keyword set or the returned list before a fetch).
    """
    nodes = [{
        "type": "terminal", "service": "full_text",
        "parameters": {"value": '"T cell receptor" AND "MHC" AND "peptide"'},
    }]
    if after_date:
        nodes.append({
            "type": "terminal", "service": "text",
            "parameters": {"attribute": "rcsb_accession_info.initial_release_date",
                           "operator": "greater_or_equal", "value": after_date},
        })
    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}},
    }
    try:
        r = requests.post(RCSB_SEARCH, json=query, timeout=timeout)
        r.raise_for_status()
    except requests.RequestException:
        return []
    return [hit["identifier"] for hit in r.json().get("result_set", [])]


def fetch_ids(ids: list[str], dest: Path | None = None, organism: str = "human") -> dict:
    """Download ``ids`` from RCSB into ``dest`` as ``.cif.gz``, keep only complete complexes.

    Incomplete (missing one of the 5 required chains) or unparseable downloads are removed.
    Returns a summary dict ``{requested, downloaded, complete, kept}``.
    """
    dest = dest or recent_dir()
    dest.mkdir(parents=True, exist_ok=True)
    downloaded = [p for pid in ids if (p := _download_cif_gz(pid, dest)) is not None]
    complete = _validate_complete(downloaded, organism=organism)
    kept = 0
    for p in downloaded:
        if complete.get(p.name):
            kept += 1
        else:
            p.unlink(missing_ok=True)  # drop incomplete / unannotatable
    return {"requested": len(ids), "downloaded": len(downloaded),
            "complete": sum(complete.values()), "kept": kept, "dest": str(dest)}


def native2026_ids() -> list[str]:
    """PDB ids of the local Native2026 set (the seed for a ``pdb_recent`` refresh)."""
    from .paths import native_dir
    from .structure import structure_id_from_path

    d = native_dir()
    return sorted({structure_id_from_path(p) for p in d.glob("*") if p.is_file()}) if d.is_dir() else []
