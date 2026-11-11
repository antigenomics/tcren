"""Render one selected structure using exported tables, never molecular geometry."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import typer

from .views import LOOPS, PARTNERS


def read_structure_tables(directory, identity):
    directory = Path(directory)
    entries = json.loads((directory / "index.json").read_text())
    entry = next((e for e in entries if e["id"] == identity), None)
    if entry is None:
        raise ValueError(f"unknown structure ID: {identity}")
    tables = {}
    for name, filename in entry["tables"].items():
        path = (directory / filename).resolve()
        if not path.is_relative_to(directory.resolve()):
            raise ValueError("table path escapes export directory")
        tables[name] = pl.read_parquet(path).filter(pl.col("structure_id") == identity)
    return tables, entry, json.loads((directory / "manifest.json").read_text())


def plot_views(directory, identity, output, *, contact_mode="all", metric="ca_distance",
               hide_rows=(), hide_columns=()):
    """Return two new PNG paths; absent pairs/cells are masked, never zero-filled."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ValueError("PNG rendering requires the tcren[viz] optional dependencies") from exc

    if contact_mode not in {"all", "pairs"} or metric not in {"ca_distance", "atom_contact_count"}:
        raise ValueError("invalid contact mode or heatmap metric")
    if set(hide_rows) - set(PARTNERS) or set(hide_columns) - set(LOOPS):
        raise ValueError("unknown facet requested")
    columns = [c for c in LOOPS if c not in hide_columns]
    rows = [r for r in PARTNERS if r not in hide_rows]
    if not rows or not columns:
        raise ValueError("at least one row and column must remain visible")
    output = Path(output)
    tables, entry, manifest = read_structure_tables(directory, identity)
    output.mkdir(parents=True, exist_ok=False)
    cells = tables["topography_cells"]
    summary = pl.read_parquet(Path(directory) / "summary.parquet").filter(pl.col("structure_id") == identity)
    clipping = (summary["surface_may_be_clipped"][0]
                if "surface_may_be_clipped" in summary.columns and summary.height else False)
    ny, nx = manifest["grid_ny_nx"]
    extent = manifest["extent_xy"]
    paths = [output / "topography.png", output / "contact-heatmaps.png"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, channel in zip(axes, ["hydropathy", "surface_z"]):
        values = np.full((ny, nx), np.nan)
        for record in cells.iter_rows(named=True):
            if record[channel] is not None:
                values[record["iy"], record["ix"]] = record[channel]
        limits = dict(vmin=-4.5, vmax=4.5, cmap="BrBG_r") if channel == "hydropathy" else dict(cmap="viridis")
        image = ax.imshow(np.ma.masked_invalid(values), origin="lower", extent=extent, **limits)
        ax.set(xlabel="Canonical x (Å), peptide N → C", ylabel="Canonical y (Å)",
               title="Hydropathy (KD)" if channel == "hydropathy" else "Exposed surface z (Å)")
        fig.colorbar(image, ax=ax, shrink=.8)
    warning = " | occupied grid boundary: inspect extent" if clipping else ""
    fig.suptitle(f"{identity}\n{entry['task']} ranking: varying partner, {manifest['definitions'][entry['task']]}{warning}", fontsize=10)
    fig.savefig(paths[0], dpi=160)
    plt.close(fig)
    pairs = tables["contact_heatmap"]
    maximum = max(1, pairs[metric].max() or 1) if metric == "atom_contact_count" else 25
    fig, axes = plt.subplots(len(rows), len(columns), figsize=(3*len(columns), 3.3*len(rows)),
                             squeeze=False, constrained_layout=True)
    for j, partner in enumerate(rows):
        for i, loop in enumerate(columns):
            ax = axes[j, i]
            selected = pairs.filter((pl.col("loop_facet") == loop) & (pl.col("partner_facet") == partner))
            ax.set_title(loop, fontsize=10)
            if selected.is_empty():
                ax.text(.5, .5, "Absent facet", ha="center", va="center", transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            shape = (selected["partner_index"].max()+1, selected["tcr_index"].max()+1)
            values = np.full(shape, np.nan)
            for record in selected.iter_rows(named=True):
                if record[metric] is not None and (contact_mode == "all" or record["contacting"]):
                    values[record["partner_index"], record["tcr_index"]] = record[metric]
            image = ax.imshow(np.ma.masked_invalid(values), origin="upper", aspect="auto",
                              vmin=0, vmax=maximum, cmap="viridis_r" if metric == "ca_distance" else "magma")
            tcr = selected.unique("tcr_index").sort("tcr_index")
            partner_res = selected.unique("partner_index").sort("partner_index")
            for axis, data, prefix in [("x", tcr, "tcr"), ("y", partner_res, "partner")]:
                step = max(1, (data.height+15)//16)
                ticks = list(range(0, data.height, step))
                labels = [f"{data[prefix+'_aa'][k]}{data[prefix+'_number'][k]}{data[prefix+'_insertion'][k]}" for k in ticks]
                getattr(ax, f"set_{axis}ticks")(ticks, labels, rotation=90 if axis == "x" else 0, fontsize=7)
            if i == 0:
                ax.set_ylabel(partner)
            ax.set_xlabel("TCR N → C", fontsize=8)
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    scalar = ScalarMappable(norm=Normalize(0, maximum), cmap="viridis_r" if metric == "ca_distance" else "magma")
    fig.colorbar(scalar, ax=axes.ravel().tolist(), shrink=.6,
                 label="Cα distance (Å), saturated at 25" if metric == "ca_distance" else "Heavy-atom contact witnesses")
    fig.suptitle(f"{identity} | {contact_mode} pairs | grey/blank = missing or masked; within-loop order N → C", fontsize=10)
    for ax in axes.ravel():
        ax.set_facecolor("#eeeeee")
    fig.savefig(paths[1], dpi=160)
    plt.close(fig)
    return paths


def plot_command(
    tables: Path = typer.Option(..., "--tables", exists=True, file_okay=False),
    structure_id: str = typer.Option(..., "--structure-id"),
    output: Path = typer.Option(..., "--output", "-o"),
    contact_mode: str = typer.Option("all", "--contact-mode", help="all or pairs (contacting pairs only)"),
    metric: str = typer.Option("ca_distance", "--metric", help="ca_distance or atom_contact_count"),
    hide_row: list[str] = typer.Option([], "--hide-row"),
    hide_column: list[str] = typer.Option([], "--hide-column"),
):
    """Plot selected-structure topography and contact facets from exported tables."""
    try:
        paths = plot_views(tables, structure_id, output, contact_mode=contact_mode, metric=metric,
                           hide_rows=hide_row, hide_columns=hide_column)
    except (ValueError, OSError, KeyError, pl.exceptions.PolarsError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    for path in paths:
        typer.echo(path)
