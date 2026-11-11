"""Bounded exports for table-owned topography and contact-map consumers."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import polars as pl
import typer

from .tables import write_table_partition
from .views import LOOPS, PARTNERS, SCHEMAS, validate_grid, view_tables

SUMMARY_SCHEMA = {
    **dict.fromkeys(("structure_id", "reference", "task"), pl.String),
    **dict.fromkeys(("grid_cells", "occupied_cells", "occupied_boundary_cells", "varying_residues",
                     "residue_pairs", "resolved_ca_pairs", "contacting_pairs", "atom_contact_witnesses"), pl.Int64),
    "surface_may_be_clipped": pl.Boolean,
    **dict.fromkeys(("occupancy", "surface_z_min", "surface_z_max", "surface_z_mean",
                     "occupied_cell_mean_hydropathy"), pl.Float64),
}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json(path, value):
    Path(path).write_text(json.dumps(value, allow_nan=False, indent=2), encoding="utf-8")


def summarize_tables(tables):
    """Descriptive values computed from exported rows, without assay labels."""
    cells, residues, pairs = (tables[n] for n in
                             ("topography_cells", "topography_residues", "contact_heatmap"))
    occupied = cells.filter(pl.col("occupied"))
    boundary = occupied.filter(
        (pl.col("ix") == 0) | (pl.col("iy") == 0) |
        (pl.col("ix") == cells["ix"].max()) | (pl.col("iy") == cells["iy"].max())).height
    values = occupied.select(
        pl.col("surface_z").min().alias("surface_z_min"),
        pl.col("surface_z").max().alias("surface_z_max"),
        pl.col("surface_z").mean().alias("surface_z_mean"),
        pl.col("hydropathy").mean().alias("occupied_cell_mean_hydropathy"),
    ).row(0, named=True)
    return dict(grid_cells=cells.height, occupied_cells=occupied.height,
                occupied_boundary_cells=boundary, surface_may_be_clipped=boundary > 0,
                occupancy=occupied.height/cells.height if cells.height else None,
                varying_residues=residues.height, residue_pairs=pairs.height,
                resolved_ca_pairs=pairs["has_both_ca"].sum(),
                contacting_pairs=pairs["contacting"].sum(),
                atom_contact_witnesses=pairs["atom_contact_count"].sum(), **values)


def export_views(exports, output, *, grid=(64, 64), extent=(-40., 40., -40., 40.), probe=1.4):
    """Read frozen scenes once; Parquet owns all derived data and serving caches.

    One structure per partition bounds memory and supports selective inspection.
    The cache and summaries are deliberately constructed by reading the written
    Parquet, so neither plotting backend can silently introduce another geometry.
    """
    validate_grid(grid, extent, probe)
    exports, output = Path(exports), Path(output)
    entries = json.loads((exports / "index.json").read_text())
    if (not entries or any(not isinstance(e.get("id"), str) or not e["id"].strip() for e in entries)
            or len({e["id"] for e in entries}) != len(entries)):
        raise ValueError("input index must contain unique nonempty structure IDs")
    output.mkdir(parents=True, exist_ok=False)
    (output / "cache").mkdir()
    for name in SCHEMAS:
        (output / name).mkdir()
    records, files, summaries = [], [], []
    for i, entry in enumerate(entries):
        source = (exports / entry["file"]).resolve()
        if not source.is_relative_to(exports.resolve()):
            raise ValueError("scene path escapes input directory")
        with gzip.open(source, "rt") as stream:
            scene = json.load(stream)
        if scene.get("schema") != "tcren.complementarity/2" or scene["id"] != entry["id"]:
            raise ValueError("scene schema or identity differs from input index")
        tables, facets = view_tables(scene, grid=grid, extent=extent, probe=probe)
        paths, hashes, owning = {}, {}, {}
        for name, table in tables.items():
            path = f"{name}/part-{i:05d}.parquet"
            tsv = write_table_partition(table, output, path)
            paths[name], hashes[name] = path, _sha(output / path)
            owning[name] = pl.read_parquet(output / path)
            files.append(dict(table=name, file=path, tsv_file=tsv, rows=table.height,
                              sha256=hashes[name], schema={k: str(v) for k, v in table.schema.items()}))
        cache = f"cache/structure-{i:05d}.json.gz"
        payload = dict(structure_id=scene["id"], sources=paths, source_sha256=hashes,
                       facets=facets, tables={k: t.to_dicts() for k, t in owning.items()})
        (output / cache).write_bytes(gzip.compress(json.dumps(payload, allow_nan=False,
                                               separators=(",", ":")).encode(), mtime=0))
        summaries.append(dict(structure_id=scene["id"], reference=scene["frame"]["reference"],
                              task=scene["task"], **summarize_tables(owning)))
        records.append(dict(id=scene["id"], display_label=entry.get("display_label", entry.get("label", scene["id"])), reference=scene["frame"]["reference"],
                            task=scene["task"], tables=paths, cache=cache, facets=facets,
                            scene_sha256=_sha(source)))
    summary = pl.DataFrame(summaries, schema=SUMMARY_SCHEMA)
    write_table_partition(summary, output, "summary.parquet")
    _json(output / "summary.json", summaries)
    _json(output / "index.json", records)
    for source_name, target_name in [("views.html", "index.html"), ("views-browser.js", "views-browser.js")]:
        (output / target_name).write_text(Path(__file__).with_name(source_name).read_text(), encoding="utf-8")
    _json(output / "manifest.json", dict(
        schema="tcren.complementarity.views/1", source_index_sha256=_sha(exports / "index.json"),
        grid_ny_nx=grid, extent_xy=extent, probe_angstrom=probe, files=files,
        loop_order=LOOPS, partner_order=PARTNERS,
        units=dict(coordinates="angstrom", ca_distance="angstrom", hydropathy="Kyte-Doolittle"),
        definitions=dict(tcr="lower sphere envelope of TCR variable domains, facing pMHC",
                         epitope="upper sphere envelope of peptide, facing TCR",
                         height="signed canonical z; not an intrinsic bulge above a fitted groove",
                         face_apex_z="extreme atom-center z, without sphere radius or probe",
                         hydropathy="exposed sphere owner's residue KD; unknown is null",
                         mean_hydropathy="area-weighted over occupied cells with known KD",
                         contacts="existing heavy-atom noncovalent witnesses; not CA cutoff",
                         sequence_order="N to C within each facet; beta facet order reversed",
                         empty="all grid cells retained; absent facets explicit; unresolved CA null")))
    return records


def views_command(
    exports: Path = typer.Option(..., "--exports", exists=True, file_okay=False),
    output: Path = typer.Option(..., "--output", "-o"),
    grid: int = typer.Option(64, "--grid", min=1, max=256),
    extent: tuple[float, float, float, float] = typer.Option((-40., 40., -40., 40.), "--extent"),
    probe: float = typer.Option(1.4, "--probe", min=0),
):
    """Export long topography/CA-contact tables from metadata-derived scenes."""
    try:
        records = export_views(exports, output, grid=(grid, grid), extent=extent, probe=probe)
    except (ValueError, OSError, KeyError, pl.exceptions.PolarsError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(f"{len(records)} structures: typed Parquet/TSV.gz, summaries and table-derived JSON caches")
