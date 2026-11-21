"""Fetch and verify the pinned static dashboard template."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen

ASSETS = Path(__file__).parent
REPOSITORY = "antigenomics/tcren-dashboard"


def _lock():
    return json.loads((ASSETS / "template-lock.json").read_text(encoding="utf-8"))


def verify_template(source: Path = ASSETS) -> None:
    """Verify the five assets against the packaged lock."""
    for name, expected in _lock()["sha256"].items():
        if hashlib.sha256((source / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"dashboard template checksum mismatch: {name}")


def bootstrap_template(output: Path, source: Path | None = None) -> str:
    """Download the exact upstream revision to a new directory."""
    if output.exists() or output.is_symlink():
        raise ValueError("template output already exists")
    lock = _lock()
    revision = lock["revision"]
    base = f"https://raw.githubusercontent.com/{REPOSITORY}/{revision}/assets"
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temporary:
        staged = Path(temporary) / "assets"
        staged.mkdir()
        for name in lock["sha256"]:
            if source is None:
                with urlopen(f"{base}/{name}", timeout=30) as response:
                    payload = response.read(2**20 + 1)
            else:
                payload = (source / name).read_bytes()
            if len(payload) > 2**20:
                raise ValueError(f"dashboard template asset too large: {name}")
            (staged / name).write_bytes(payload)
        verify_template(staged)
        if output.exists() or output.is_symlink():
            raise ValueError("template output already exists")
        staged.rename(output)
    return revision
