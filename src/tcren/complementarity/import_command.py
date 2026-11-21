"""Build a static dashboard for one explicitly declared scoring cohort."""

import json
import math
import shutil
from html import escape
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory
from urllib.parse import urlencode

import polars as pl
import typer

from ..screen import MIN_COHORT
from ..structure.io import structure_stem
from .command import complementarity_command
from .dashboard import STYLESHEET_LINK, install_dashboard_style
from .inputs import metadata_inputs
from .template import bootstrap_template, verify_template
from .view_export import export_views

dashboard_app = typer.Typer(help="Build a static dashboard from model metadata.")


def _metadata(path, binary):
    models = metadata_inputs(path, id_column="structure_id")
    table = pl.read_csv(path, separator="\t", infer_schema=False)
    columns = ("iptm", "plddt", "label")
    if missing := set(columns) - set(table.columns):
        raise ValueError(f"missing metadata columns: {', '.join(sorted(missing))}")
    for model, row in zip(models, table.select(columns).iter_rows(named=True), strict=True):
        for key, bounds in zip(columns, ((0, 1), (0, 100), None), strict=True):
            try:
                value = float(row[key])
            except (TypeError, ValueError):
                value = math.nan
            if (not math.isfinite(value) or (bounds and not bounds[0] <= value <= bounds[1])
                    or (key == "label" and binary and value not in (0, 1))):
                raise ValueError(f"metadata row {model.row}: invalid {key}")
    metadata = table.select("structure_id", *columns).cast({
        "iptm": pl.Float64, "plddt": pl.Float64, "label": pl.Float64,
    })
    return models, metadata


def _write_index(output, rows, task, kind, references):
    cells = []
    for row in rows.sort("screen_score", descending=True, nulls_last=True).iter_rows(named=True):
        identity = row["structure_id"]
        query = escape(urlencode({"group": references[identity], "id": identity}), quote=True)
        links = f'<a href="maps/?{query}">Molecular</a> · <a href="views/?{query}">Table views</a>'
        values = (row["iptm"], row["plddt"], row["label"],
                  row["screen_score"], row["score_unavailable"])
        cells.append(f"<tr><td>{escape(identity)}</td><td>{links}</td>" + "".join(
            f"<td>{escape(str(value)) if value is not None else ''}</td>" for value in values
        ) + "</tr>")
    install_dashboard_style(output)
    copyfile(Path(__file__).with_name("view-controls.js"), output / "view-controls.js")
    (output / "index.html").write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width">'
        '<meta name="color-scheme" content="light dark">'
        f'<title>TCRen dashboard</title>{STYLESHEET_LINK}'
        '<body class="ranking-home">'
        '<nav id="view-controls" aria-label="Viewer controls" data-theme-only></nav>'
        '<script type="module" src="view-controls.js"></script>'
        f'<h1>TCRen {task} ranking</h1>'
        f"<p>{rows.height} models · {kind} labels · <a href='scores.tsv'>Download scores</a></p>"
        '<div class="score-scroll" role="region" aria-label="Model scores" tabindex="0">'
        '<table><thead><tr><th scope="col">Structure ID</th><th scope="col">Views</th>'
        '<th scope="col">ipTM</th>'
        '<th scope="col">pLDDT</th><th scope="col">Label</th><th scope="col">Score</th>'
        '<th scope="col">Diagnostic</th></tr></thead>'
        f"<tbody>{''.join(cells)}</tbody></table></div></body></html>",
        encoding="utf-8",
    )
    for directory, other, label in (("maps", "views", "Table views"), ("views", "maps", "Molecular")):
        (output / directory / "navigation.json").write_text(json.dumps({"links": [
            {"label": "Ranking", "href": "../"},
            {"label": label, "href": f"../{other}/", "model": True},
        ]}), encoding="utf-8")


@dashboard_app.command("bootstrap")
def bootstrap_dashboard(
    output: Path = typer.Option(..., "--output", "-o", help="New template asset directory"),
    source: Path = typer.Option(None, help="Assets directory from a checkout of the pinned revision"),
):
    """Fetch and verify the viewer template pinned by this TCRen version."""
    try:
        revision = bootstrap_template(output, source)
    except (ValueError, OSError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(f"Wrote dashboard template {revision} to {output}")


@dashboard_app.command("import")
def import_dashboard(
    metadata: Path,
    task: str = typer.Option(..., help="tcr (fixed pMHC) or epitope (fixed receptor)"),
    label_kind: str = typer.Option(..., help="binary (0/1) or continuous"),
    output: Path = typer.Option(..., "--output", "-o", help="New dashboard directory"),
    threads: int = typer.Option(0, "--threads", "-t", min=0),
):
    """Score one cohort from a five-column TSV and export both browser views."""
    from ..cli import screen

    if task not in {"tcr", "epitope"} or label_kind not in {"binary", "continuous"}:
        raise typer.BadParameter("use --task tcr/epitope and --label-kind binary/continuous")
    if output.exists() or output.is_symlink():
        raise typer.BadParameter("output directory already exists")
    if metadata.suffix != ".tsv":
        raise typer.BadParameter("metadata must be a .tsv file")
    try:
        verify_template()
        models, labels = _metadata(metadata, label_kind == "binary")
    except (ValueError, OSError, pl.exceptions.PolarsError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".tcren-import-", dir=output.parent) as temporary:
        work = Path(temporary)
        bundle = work / "dashboard"
        bundle.mkdir()
        inputs = work / "inputs"
        inputs.mkdir()
        score_ids, staged_paths = [], []
        for ordinal, model in enumerate(models):
            extension = model.path.name[len(structure_stem(model.path)):]
            identity = f"model{ordinal:08d}"
            target = inputs / (identity + extension)
            try:
                target.hardlink_to(model.path)
            except OSError:
                shutil.copyfile(model.path, target)
            score_ids.append(identity)
            staged_paths.append(target.relative_to(work).as_posix())
        identities = labels.with_columns(pl.Series("complex.id", score_ids))
        confidence = work / "confidence.tsv"
        identities.select("complex.id", "iptm", "plddt").write_csv(confidence, separator="\t")
        scored = work / "scored.tsv"
        screen(structures=[str(inputs)], out=scored, task=task, confidence=confidence,
               confidence_context="cohort", reference=None, min_cohort=MIN_COHORT,
               on_flag="report", weight=None, template_metadata=None, template_status=None,
               organism="human", threads=threads)
        scores = pl.read_csv(scored, separator="\t", infer_schema_length=None)
        if set(scores["complex.id"]) != set(score_ids):
            raise typer.BadParameter("screen results do not cover every metadata model")
        observed = identities.join(scores.drop("iptm", "plddt", strict=False),
                                   on="complex.id", validate="1:1", maintain_order="left").drop("complex.id")
        observed.write_csv(bundle / "scores.tsv", separator="\t")

        scene_inputs = work / "structures.tsv"
        labels.select("structure_id").with_columns(
            pl.Series("structure_path", staged_paths)
        ).write_csv(scene_inputs, separator="\t")
        maps = bundle / "maps"
        complementarity_command(
            structures=None, output=maps, task=task, origin="center-of-mass", cutoff=5.,
            threads=threads, legacy=False, scratch=None, exclude_id=None, metadata=scene_inputs,
            path_column="structure_path", id_column="structure_id", member_column=None, portable=True,
        )
        scenes = json.loads((maps / "index.json").read_text(encoding="utf-8"))
        if {scene["id"] for scene in scenes} != set(labels["structure_id"]):
            raise typer.BadParameter("molecular scenes do not cover every metadata model")
        references = {scene["id"]: scene["reference"] for scene in scenes}
        export_views(maps, bundle / "views")
        _write_index(bundle, observed, task, label_kind, references)
        if output.exists() or output.is_symlink():
            raise typer.BadParameter("output directory already exists")
        bundle.rename(output)
    typer.echo(f"Wrote dashboard for {len(models)} models to {output}")
