"""Write compact offline dashboards with compressed scenes loaded on demand."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

from .tables import SceneTables


def write_dashboard(scenes, output: str | Path) -> Path:
    """Stream scenes to a new directory, preserving existing deliverables.

    Serve the output using ``python -m http.server --directory OUTPUT``. No CDN,
    remote scripts or network service is needed. Scene IDs are labels, never paths.
    """
    output = Path(output)
    template = Path(__file__).with_name("dashboard.html").read_text(encoding="utf-8")
    viewer = Path(__file__).with_name("viewer.js").read_text(encoding="utf-8")
    output.mkdir(parents=True, exist_ok=False)
    entries = []
    tables = SceneTables(output / "tables")
    seen = set()
    for i, scene in enumerate(scenes):
        if scene.get("schema") != "tcren.complementarity/2":
            raise ValueError("dashboard requires tcren.complementarity/2 scenes")
        if scene["id"] in seen:
            raise ValueError(f"duplicate structure ID: {scene['id']}")
        seen.add(scene["id"])
        name = f"scene-{i}.json.gz"
        payload = json.dumps(scene, allow_nan=False, separators=(",", ":")).encode()
        (output / name).write_bytes(gzip.compress(payload, mtime=0))
        tables.add(scene)
        entries.append(dict(id=scene["id"], reference=scene["frame"]["reference"],
                            label=scene.get("group_label", scene["frame"]["reference"]),
                            search=" ".join(c["sequence"] for c in scene.get("components", [])),
                            file=name, warnings=scene["warnings"]))
    (output / "index.json").write_text(json.dumps(entries, allow_nan=False), encoding="utf-8")
    tables.close()
    (output / "index.html").write_text(template, encoding="utf-8")
    (output / "viewer.js").write_text(viewer, encoding="utf-8")
    return output / "index.html"
