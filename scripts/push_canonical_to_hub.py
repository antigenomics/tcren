#!/usr/bin/env python
"""INTERNAL dev tool — upload an oriented structure set to the HF dataset.

Not part of the user-facing CLI: maintainers run this to refresh
``isalgo/tcren_structures`` after rebuilding a canonical set with ``tcren orient``.
Requires a logged-in HF token with write access.

    python scripts/push_canonical_to_hub.py data/Canonical2026 --folder Canonical2026
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir", type=Path, help="local dir of oriented structures to upload")
    ap.add_argument("--repo", default="isalgo/tcren_structures", help="HF dataset repo id")
    ap.add_argument("--folder", default="Canonical2026", help="path_in_repo subfolder")
    args = ap.parse_args()

    from huggingface_hub import HfApi

    HfApi().upload_folder(
        folder_path=str(args.out_dir),
        path_in_repo=args.folder,
        repo_id=args.repo,
        repo_type="dataset",
    )
    print(f"pushed {args.out_dir} -> {args.repo}/{args.folder}")


if __name__ == "__main__":
    main()
