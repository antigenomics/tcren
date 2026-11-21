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
        (output / file).write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(",", ":")), encoding="utf-8")
        label = models[0].get("group_label") or models[0].get("label") or reference
        index.append(dict(id=reference, label=label, file=file))
    directory = dict(groups=index, count=len(entries))
    if manifest is not None:
        directory["manifest"] = manifest
        directory["ranking"] = ranking
    (output / "groups.json").write_text(json.dumps(directory, ensure_ascii=False, allow_nan=False, separators=(",", ":")), encoding="utf-8")


STYLESHEET_LINK = '<link rel="stylesheet" href="dashboard.css">'


def install_dashboard_style(output: Path) -> None:
    """Copy the shared stylesheet without following a prepared asset symlink."""
    target = output / "dashboard.css"
    target.unlink(missing_ok=True)
    shutil.copyfile(Path(__file__).with_name("dashboard.css"), target)


def install_viewer_assets(output: Path, *, tables: bool = False, template: Path | None = None) -> None:
    """Install the packaged page and scripts into a prepared static site."""
    output = Path(output)
    html, script = ("views.html", "views-browser.js") if tables else ("dashboard.html", "viewer.js")
    assets = Path(__file__).parent
    page = (assets / html).read_text(encoding="utf-8")
    if template is not None:
        prepared = Path(template).read_text(encoding="utf-8")
        controls = page[page.index('<nav id="view-controls"'):page.index('</script>') + len('</script>')]
        if "</style>" not in prepared or "<h1>" not in prepared:
            raise ValueError("prepared viewer template requires a style block and heading")
        page = prepared.replace("<style", STYLESHEET_LINK + "<style", 1)
        page = page.replace("<h1>", controls + "<h1>", 1)
    install_dashboard_style(output)
    (output / "index.html").write_text(page, encoding="utf-8")
    for source in (script, "view-controls.js"):
        shutil.copyfile(assets / source, output / source)


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


def write_benchmark_page(output: Path, rows, *, title: str, description: str,
                         score_label: str, pairing_labels: tuple[str, str]) -> None:
    """Render a static comparison from recorded rows in a caller-selected stratum."""
    from html import escape

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    cells = []
    epitopes = set()
    for row in rows:
        epitope = row["epitope"]
        epitopes.add(epitope)
        pairing = row["y"]
        if pairing not in {"0", "1"}:
            raise ValueError("benchmark pairing must be 0 or 1")
        unavailable = row["score_unavailable"]
        if row["structural_excluded"].lower() == "true":
            status = "Excluded by structural review"
        elif row["no_interface"].lower() == "true":
            status = "No TCR–peptide interface"
        else:
            status = unavailable or "Scored"
        values = (epitope, row["id"], pairing_labels[int(pairing)], row["iptm"], row["screen_score"])
        cells.append(
            f'<tr data-pairing="{pairing}" data-epitope="{escape(epitope)}">'
            + "".join(f"<td>{escape(value or '—')}</td>" for value in values)
            + f'<td title="{escape(unavailable)}">{escape(status)}</td></tr>'
        )
    options = "".join(f"<option>{escape(epitope)}</option>" for epitope in sorted(epitopes))
    (output / "index.html").write_text(f"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
{STYLESHEET_LINK}
<body class="benchmark-view">
<nav id="view-controls" aria-label="Viewer controls"></nav>
<h1>{escape(title)}</h1><p>{escape(description)}</p>
<label>Find a record <input id="benchmark-search" type="search" placeholder="Epitope or model ID"></label>
<label>Pairing <select id="benchmark-pairing"><option value="">All pairings</option>
<option value="1">{escape(pairing_labels[1])}</option><option value="0">{escape(pairing_labels[0])}</option>
</select></label>
<label>Epitope <select id="benchmark-epitope"><option value="">All epitopes</option>{options}</select></label>
<p id="benchmark-visible" role="status">{len(cells)} rows</p>
<div class="table-scroll"><table id="benchmark-table"><thead><tr>
<th scope="col">Epitope</th><th scope="col">Model ID</th><th scope="col">Pairing</th>
<th scope="col" aria-sort="none"><button data-sort="3">ipTM</button></th>
<th scope="col" aria-sort="none"><button data-sort="4">{escape(score_label)}</button></th>
<th scope="col">Status</th></tr></thead><tbody>{''.join(cells)}</tbody></table></div>
<script type="module" src="view-controls.js"></script></body></html>
""", encoding="utf-8")
    install_dashboard_style(output)
    shutil.copyfile(Path(__file__).with_name("view-controls.js"), output / "view-controls.js")


def serve_static(output: Path, *, port: int = 8772, asset_roots=()) -> None:
    """Serve a static site with explicitly allowed prepared asset roots on loopback."""
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    output = Path(output).resolve()
    roots = (output, *(Path(root).resolve() for root in asset_roots))

    class Handler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=output, **kwargs)

        def translate_path(self, path):
            target = Path(super().translate_path(path)).resolve()
            if not any(target.is_relative_to(root) for root in roots):
                return "/dev/null/unavailable"
            return str(target)

        def list_directory(self, path):
            self.send_error(404)
            return None

    with ThreadingHTTPServer(("127.0.0.1", port), Handler) as server:
        print(f"TCRen static viewer at http://127.0.0.1:{server.server_port}/", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
