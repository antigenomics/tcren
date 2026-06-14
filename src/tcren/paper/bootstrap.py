"""Bootstrap data for the Nat Comput Sci 2022 reproduction (``notebooks/natcompsci2022/``).

Fetches the structure sets from the Hugging Face dataset ``isalgo/tcren_structures`` and
the pinned vdjdb release, and stages the paper's non-structure data plus the legacy
(mir/R) text outputs used as the regression oracle. Network access uses ``requests`` only
(no ``huggingface_hub`` dependency). All copied csv/tsv/txt are gzipped; the downloaded HF
structures are gitignored.
"""

from __future__ import annotations

import gzip
import io
import shutil
import time
import zipfile
from pathlib import Path

import requests

_REPO = Path(__file__).resolve().parents[3]
PAPER_DIR = _REPO / "notebooks" / "natcompsci2022"

HF_REPO = "isalgo/tcren_structures"
HF_FOLDERS = ("Native2022", "PolyV2022", "Bobisse", "Bigot")
_HF_TREE = "https://huggingface.co/api/datasets/{repo}/tree/main/{folder}?recursive=true"
_HF_RESOLVE = "https://huggingface.co/datasets/{repo}/resolve/main/{path}"

VDJDB_DATE = "2022-03-30"
_VDJDB_URL = "https://github.com/antigenomics/vdjdb-db/releases/download/{date}/vdjdb-{date}.zip"

# Paper non-structure data: repo-relative path -> destination name under data/.
_PAPER_DATA = [
    "Birnbaum.tsv", "MJ_Keskin_potentials.csv", "PDB_MHC_annotation.csv", "PDB_date.csv",
    "benchmark_candidate_epitopes_IEDB.txt", "iedb_slim.csv", "TCRen_potential.csv",
    "Bobisse/Bobisse_peptides.tsv", "Bobisse/Bobisse_candidate_epitopes.txt",
    "Bigot/Bigot_candidate_epitopes.txt", "Bigot/Bigot_cognate_epitopes.csv",
]
# Legacy mir annotation output directories (per structure set).
_LEGACY_ANNOTATION = {
    "Native2022": "data/output_TCRen/structures_annotation",
    "PolyV2022": "data/TCRpMHCmodels/output_TCRen",
    "Bobisse": "data/Bobisse/output_TCRen_TCRpMHCmodels/structures_annotation",
    "Bigot": "data/Bigot/output_TCRen",
}
_LEGACY_DERIVED = [
    "contact_maps_PDB.csv", "summary_PDB_structures.csv", "contacts_PDB.csv",
    "TCRen_potential.csv",
]


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
    out = data_dir / "structures"
    out.mkdir(parents=True, exist_ok=True)
    try:
        from huggingface_hub import snapshot_download  # noqa: PLC0415

        snapshot_download(
            repo_id=HF_REPO, repo_type="dataset", local_dir=str(out),
            allow_patterns=[f"{f}/*" for f in folders],
        )
    except ImportError:
        for folder in folders:
            for path in _hf_files(folder, timeout):
                dest = out / path
                if dest.exists() and not force:
                    continue
                dest.parent.mkdir(parents=True, exist_ok=True)
                _download_with_requests(path, dest, timeout)
    return {f: len(list((out / f).glob("*.pdb"))) for f in folders}


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


def copy_paper_data(data_dir: Path, repo_data: Path | None = None) -> int:
    """Copy + gzip the paper's non-structure data + ``source_data/`` into ``<data_dir>/``."""
    repo_data = repo_data or (_REPO / "data")
    n = 0
    for rel in _PAPER_DATA:
        src = repo_data / rel
        if src.exists():
            _gzip_copy(src, data_dir / (rel + ".gz"))
            n += 1
    for sub in ("benchmark_other_tools", "source_data", "TCRpMHCmodels"):
        for src in (repo_data / sub).rglob("*"):
            if src.is_file() and (sub != "TCRpMHCmodels" or src.name.endswith("-complex-templates.csv")):
                _gzip_copy(src, data_dir / sub / (src.relative_to(repo_data / sub).as_posix() + ".gz"))
                n += 1
    return n


def copy_legacy_results(paper_dir: Path = PAPER_DIR, repo_data: Path | None = None) -> int:
    """Stage the legacy mir/R text outputs (regression oracle) into ``results_legacy/`` (gzipped)."""
    repo_data = repo_data or (_REPO / "data")
    out = paper_dir / "results_legacy"
    n = 0
    for fig in (repo_data / "source_data").glob("*.csv"):
        _gzip_copy(fig, out / "source_data" / (fig.name + ".gz"))
        n += 1
    for name in _LEGACY_DERIVED:
        src = repo_data / name
        if src.exists():
            _gzip_copy(src, out / (name + ".gz"))
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
    paper_dir: Path = PAPER_DIR, structures: bool = True, vdjdb: bool = True,
    data: bool = True, legacy: bool = True,
) -> dict:
    """Run the full bootstrap; return a summary dict."""
    data_dir = paper_dir / "data"
    summary: dict = {}
    if structures:
        summary["structures"] = fetch_hf_structures(data_dir)
    if vdjdb:
        summary["vdjdb"] = str(fetch_vdjdb(data_dir))
    if data:
        summary["paper_data_files"] = copy_paper_data(data_dir)
    if legacy:
        summary["legacy_files"] = copy_legacy_results(paper_dir)
    return summary
