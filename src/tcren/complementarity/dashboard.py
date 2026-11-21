"""Write compact offline dashboards with compressed scenes loaded on demand."""
from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

from .tables import SceneTables


def write_browser_groups(output, entries, *, manifest=None, summaries=None, fixed=None, ranking=None):
    """Write a small directory and one metadata file per reference group."""
    groups = {}
    for entry in entries:
        groups.setdefault(entry["reference"], []).append(entry)
    (output / "groups").mkdir(exist_ok=True)
    index = []
    for number, (reference, models) in enumerate(groups.items()):
        file = f"groups/{number}.json"
        payload = {"entries": models}
        if summaries is not None:
            payload["summaries"] = [summaries[model["id"]] for model in models]
        if fixed and reference in fixed:
            payload["fixed"] = fixed[reference]
        (output / file).write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
        label = models[0].get("group_label") or models[0].get("label") or reference
        index.append(dict(id=reference, label=label, file=file))
    directory = dict(groups=index, count=len(entries))
    if manifest is not None:
        directory["manifest"] = manifest
        directory["ranking"] = ranking
    (output / "groups.json").write_text(json.dumps(directory, allow_nan=False), encoding="utf-8")


def install_viewer_assets(output: Path, *, tables: bool = False) -> None:
    """Install the packaged page and scripts into a prepared static site."""
    output = Path(output)
    html, script = ("views.html", "views-browser.js") if tables else ("dashboard.html", "viewer.js")
    for source, target in ((html, "index.html"), (script, script), ("view-controls.js", "view-controls.js")):
        shutil.copyfile(Path(__file__).with_name(source), output / target)


def write_dashboard(scenes, output: str | Path) -> Path:
    """Stream scenes to a new directory, preserving existing deliverables.

    Serve the output using ``python -m http.server --directory OUTPUT``. No CDN,
    remote scripts or network service is needed. Scene IDs are labels, never paths.
    """
    output = Path(output)
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
        identity = scene["id"]
        short_id = identity[:8] if len(identity) > 24 else identity
        peptides = [c["sequence"] for c in scene.get("components", []) if c["role"] == "peptide"]
        display_label = scene.get("display_label") or " · ".join([*peptides, short_id])
        entries.append(dict(id=identity, display_label=display_label, reference=scene["frame"]["reference"],
                            label=scene.get("group_label", scene["frame"]["reference"]),
                            search=" ".join(c["sequence"] for c in scene.get("components", [])),
                            file=name, warnings=scene["warnings"]))
    (output / "index.json").write_text(json.dumps(entries, allow_nan=False), encoding="utf-8")
    write_browser_groups(output, entries)
    tables.close()
    install_viewer_assets(output)
    return output / "index.html"
