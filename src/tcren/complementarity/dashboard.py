"""Write compact offline dashboards with compressed scenes loaded on demand."""
from __future__ import annotations

import gzip
import json
from pathlib import Path


def write_dashboard(scenes, output: str | Path) -> Path:
    """Stream scenes to a new directory, preserving existing deliverables.

    Serve the output using ``python -m http.server --directory OUTPUT``. No CDN,
    remote scripts or network service is needed. Scene IDs are labels, never paths.
    """
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    entries = []
    for i, scene in enumerate(scenes):
        name = f"scene-{i}.json.gz"
        payload = json.dumps(scene, allow_nan=False, separators=(",", ":")).encode()
        (output / name).write_bytes(gzip.compress(payload, mtime=0))
        entries.append(dict(id=scene["id"], reference=scene["frame"]["reference"],
                            label=scene.get("group_label", scene["frame"]["reference"]),
                            search=" ".join(c["sequence"] for c in scene.get("components", [])),
                            file=name, warnings=scene["warnings"]))
    (output / "index.json").write_text(json.dumps(entries, allow_nan=False), encoding="utf-8")
    template = Path(__file__).with_name("dashboard.html").read_text(encoding="utf-8")
    (output / "index.html").write_text(template, encoding="utf-8")
    return output / "index.html"
