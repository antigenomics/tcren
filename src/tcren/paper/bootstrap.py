"""Bootstrap data for the Nat Comput Sci reproduction notebooks.

Fetches the structure sets from the Hugging Face dataset ``isalgo/tcren_structures`` and
the pinned vdjdb release into the shared ``notebooks/data/`` dir, and stages the paper's
non-structure data there plus the legacy (mir/R) text outputs used as the regression
oracle under ``notebooks/natcompsci2022/results_legacy/``. Network access uses ``requests``
only (no ``huggingface_hub`` dependency). All copied csv/tsv/txt are gzipped; the
downloaded HF structures are gitignored.
"""

from __future__ import annotations

import gzip
import io
import shutil
import time
import zipfile
from pathlib import Path

import requests

from ..structure.io import is_structure_file

_REPO = Path(__file__).resolve().parents[3]
PAPER_DIR = _REPO / "notebooks" / "natcompsci2022"
# Shared data lives one level up (notebooks/data) so sibling notebooks can reuse it; the
# regression results stay under the paper notebook dir (results_legacy/, results_new/).
DATA_DIR = _REPO / "notebooks" / "data"

HF_REPO = "isalgo/tcren_structures"
# Native2022 = original 2022 paper PDB set (oracle); Native2026 = the comprehensive 2026
# TCR:pMHC set the new TCRen is derived from; Canonical2026 = Native2026 after canonical
# re-orientation (tcren orient). All structures on the Hub are gzipped (*.pdb.gz / *.cif.gz).
HF_FOLDERS = ("Native2022", "Native2026", "PolyV2022", "Bobisse", "Bigot")
CANONICAL_FOLDER = "Canonical2026"
_HF_TREE = "https://huggingface.co/api/datasets/{repo}/tree/main/{folder}?recursive=true"
_HF_RESOLVE = "https://huggingface.co/datasets/{repo}/resolve/main/{path}"

VDJDB_DATE = "2022-03-30"
_VDJDB_URL = "https://github.com/antigenomics/vdjdb-db/releases/download/{date}/vdjdb-{date}.zip"

# Externally published inputs the new pipeline is allowed to consume (→ notebooks/data/).
# Strictly: MJ/Keskin potentials, the Birnbaum yeast-display set, IEDB, and the per-paper
# candidate-epitope lists (Bobisse/Bigot). vdjdb + PDB release dates are fetched live.
_EXTERNAL_INPUTS = [
    "MJ_Keskin_potentials.csv", "Birnbaum.tsv",
    "benchmark_candidate_epitopes_IEDB.txt", "iedb_slim.csv",
    "Bobisse/Bobisse_peptides.tsv", "Bobisse/Bobisse_candidate_epitopes.txt",
    "Bigot/Bigot_candidate_epitopes.txt", "Bigot/Bigot_cognate_epitopes.csv",
]
# Legacy mir/R outputs — COMPARISON BASELINES ONLY (→ data_legacy/, never pipeline inputs).
_LEGACY_DERIVED = [
    "contact_maps_PDB.csv", "summary_PDB_structures.csv", "contacts_PDB.csv",
    "TCRen_potential.csv", "PDB_MHC_annotation.csv",
]
# Legacy mir annotation output directories (per structure set) → data_legacy/annotation/.
_LEGACY_ANNOTATION = {
    "Native2022": "data/output_TCRen/structures_annotation",
    "PolyV2022": "data/TCRpMHCmodels/output_TCRen",
    "Bobisse": "data/Bobisse/output_TCRen_TCRpMHCmodels/structures_annotation",
    "Bigot": "data/Bigot/output_TCRen",
}


def _gzip_copy(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with src.open("rb") as fi, gzip.open(dest, "wb") as fo:
        shutil.copyfileobj(fi, fo)


def _hf_files(folder: str, timeout: float) -> list[str]:
    resp = requests.get(_HF_TREE.format(repo=HF_REPO, folder=folder), timeout=timeout)
    resp.raise_for_status()
    return [e["path"] for e in resp.json() if e.get("type") == "file"]


def _download_with_requests(path: str, dest: Path, timeout: float, retries: int = 5) -> None:
    """Robust single-file download (handles transient SSL/CDN drops with backoff)."""
    url = _HF_RESOLVE.format(repo=HF_REPO, path=path)
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                with tmp.open("wb") as fo:
                    for chunk in r.iter_content(1 << 20):
                        fo.write(chunk)
                tmp.replace(dest)
            return
        except (requests.RequestException, OSError):
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def fetch_hf_structures(
    data_dir: Path, folders: tuple[str, ...] = HF_FOLDERS, force: bool = False,
    timeout: float = 120.0,
) -> dict[str, int]:
    """Download the HF structure folders into ``<data_dir>/structures/`` (gitignored).

    Uses ``huggingface_hub.snapshot_download`` when available (robust LFS download with
    resume + retries); otherwise falls back to a retrying ``requests`` loop. Neither is a
    hard package dependency — this is an optional reproduction tool.
    """
    out = data_dir  # per-set folders land directly under data_dir (no "structures/" wrapper)
    out.mkdir(parents=True, exist_ok=True)

    def _present(folder: str) -> bool:
        d = out / folder
        return d.is_dir() and any(p.is_file() and is_structure_file(p) for p in d.glob("*"))

    # Skip-if-present: only hit the network for folders that are missing locally. (Otherwise
    # snapshot_download HEADs every one of ~800 files each run — the dominant bootstrap cost.)
    todo = list(folders) if force else [f for f in folders if not _present(f)]
    if todo:
        try:
            from huggingface_hub import snapshot_download  # noqa: PLC0415

            snapshot_download(
                repo_id=HF_REPO, repo_type="dataset", local_dir=str(out),
                allow_patterns=[f"{f}/*" for f in todo],
            )
        except ImportError:
            for folder in todo:
                for path in _hf_files(folder, timeout):
                    dest = out / path
                    if dest.exists() and not force:
                        continue
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    _download_with_requests(path, dest, timeout)
    return {
        f: sum(1 for p in (out / f).glob("*") if p.is_file() and is_structure_file(p))
        for f in folders if (out / f).is_dir()
    }


def fetch_vdjdb(data_dir: Path, date: str = VDJDB_DATE, timeout: float = 300.0) -> Path:
    """Download the pinned vdjdb release and gzip its slim table into ``<data_dir>/vdjdb/``."""
    out = data_dir / "vdjdb" / f"vdjdb-{date}.slim.txt.gz"
    if out.exists():
        return out
    resp = requests.get(_VDJDB_URL.format(date=date), timeout=timeout)
    resp.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        slim = next(n for n in zf.namelist() if n.endswith(".slim.txt"))
        out.parent.mkdir(parents=True, exist_ok=True)
        with zf.open(slim) as fi, gzip.open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return out


_RCSB_DATES_API = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"


def fetch_pdb_dates(data_dir: Path, pdb_ids: list[str], timeout: float = 30.0) -> Path:
    """Fetch PDB initial-release dates from RCSB into ``<data_dir>/PDB_date.csv.gz``.

    RCSB release dates are external published metadata (needed for the holdout date-split).
    Queried per entry via the RCSB Data API; cached, so re-runs only fetch missing ids.
    """
    import csv

    out = data_dir / "PDB_date.csv.gz"
    have: dict[str, str] = {}
    if out.exists():
        with gzip.open(out, "rt") as fh:
            for row in csv.DictReader(fh):
                have[row["pdb.id"]] = row["release.date"]
    for pdb in pdb_ids:
        pid = pdb.lower()
        if pid in have:
            continue
        try:
            r = requests.get(_RCSB_DATES_API.format(pdb=pid), timeout=timeout)
            if r.ok:
                have[pid] = (r.json().get("rcsb_accession_info", {})
                             .get("initial_release_date", "") or "")[:10]
        except requests.RequestException:
            continue
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as fh:
        fh.write("pdb.id,release.date\n")
        for pid in sorted(have):
            fh.write(f"{pid},{have[pid]}\n")
    return out


def copy_external_inputs(data_dir: Path, repo_data: Path | None = None) -> int:
    """Stage the allowed externally published inputs (gzipped) into ``<data_dir>/``.

    These are the only non-structure inputs the reproduced pipeline may consume: the
    MJ/Keskin potentials, the Birnbaum set, IEDB, and the Bobisse/Bigot candidate lists.
    """
    repo_data = repo_data or (_REPO / "data")
    n = 0
    for rel in _EXTERNAL_INPUTS:
        src = repo_data / rel
        if src.exists():
            _gzip_copy(src, data_dir / (rel + ".gz"))
            n += 1
    return n


def copy_legacy_results(paper_dir: Path = PAPER_DIR, repo_data: Path | None = None) -> int:
    """Stage the legacy mir/R outputs into ``data_legacy/`` (gzipped) — comparison only.

    Holds the 2022 baselines the new tcren results are measured against: the published
    TCRen matrix, mir contacts, the old non-redundancy ``summary``, the legacy MHC
    annotation, the mir annotation oracle, paper ``source_data/`` and the other-tools
    (TITAN/ERGO-II) outputs. Never consumed as a pipeline input.
    """
    repo_data = repo_data or (_REPO / "data")
    out = paper_dir / "data_legacy"
    n = 0
    for fig in (repo_data / "source_data").glob("*.csv"):
        _gzip_copy(fig, out / "source_data" / (fig.name + ".gz"))
        n += 1
    for name in _LEGACY_DERIVED:
        src = repo_data / name
        if src.exists():
            _gzip_copy(src, out / (name + ".gz"))
            n += 1
    for src in (repo_data / "benchmark_other_tools").rglob("*"):
        if src.is_file():
            rel = src.relative_to(repo_data / "benchmark_other_tools").as_posix()
            _gzip_copy(src, out / "benchmark_other_tools" / (rel + ".gz"))
            n += 1
    for label, rel in _LEGACY_ANNOTATION.items():
        adir = _REPO / rel
        # general/markup/resmarkup = the annotation regression oracle; raw atomdist contacts
        # are redundant with the committed contact_maps_PDB.csv, so they are not staged.
        for tsv in ("general.txt", "markup.txt", "resmarkup.txt"):
            src = adir / tsv
            if src.exists():
                _gzip_copy(src, out / "annotation" / label / (tsv + ".gz"))
                n += 1
    return n


def bootstrap(
    paper_dir: Path = PAPER_DIR, data_dir: Path = DATA_DIR, structures: bool = True,
    canonical: bool = False, **_legacy_flags,
) -> dict:
    """Fetch the HF structure sets into ``data_dir`` (``notebooks/data`` by default).

    Each set lands in its own folder directly under ``data_dir`` (e.g. ``notebooks/data/
    Native2026/``) — no ``structures/`` wrapper — and is gitignored. The non-structure inputs
    (vdjdb, Birnbaum, MJ/Keskin, IEDB, epitope lists, PDB dates) and the legacy comparison
    baselines are already committed under ``<paper_dir>/data_legacy/``, so they are not
    re-fetched here. Pass ``canonical=True`` to also fetch the re-oriented ``Canonical2026`` set.
    """
    summary: dict = {}
    if structures:
        folders = HF_FOLDERS + ((CANONICAL_FOLDER,) if canonical else ())
        summary["structures"] = fetch_hf_structures(data_dir, folders=folders)
    return summary
