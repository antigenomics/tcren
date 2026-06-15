"""Upload oriented structures to a HuggingFace dataset folder (optional dependency)."""

from __future__ import annotations

from pathlib import Path


def push_oriented(out_dir: str | Path, repo_id: str, folder: str = "Native2026") -> None:
    """Upload the oriented PDBs in ``out_dir`` to ``repo_id`` under ``folder/`` on HF.

    Uses ``huggingface_hub`` lazily (not a hard dependency); requires a logged-in token.
    """
    from huggingface_hub import HfApi  # noqa: PLC0415

    api = HfApi()
    api.upload_folder(
        folder_path=str(out_dir),
        path_in_repo=folder,
        repo_id=repo_id,
        repo_type="dataset",
    )
