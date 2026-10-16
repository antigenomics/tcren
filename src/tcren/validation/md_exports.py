"""CLI exports of existing MD/coordinate APIs; no simulation or fitted score."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shlex
import sys
from enum import Enum

import polars as pl
import numpy as np
import typer

from .. import __version__
from ..runtime import native_threads, resolve_threads, single_threaded
from .md import compare_same_frame_energy_forces, inspect_md_trajectory
from .structures import close_atom_pairs, coordinate_diagnostics

app = typer.Typer(help="Recorded MD and atom-level diagnostics; does not run simulations.")


class MhcClass(str, Enum):
    class_i = "I"
    class_ii = "II"


def _identity(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "sha256": digest.hexdigest(), "bytes": path.stat().st_size}


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _pairs(path, output, cutoff):
    rows = close_atom_pairs(path, cutoff=cutoff)
    schema = {"distance_A": pl.Float64, "interchain": pl.Boolean,
              "possible_disulfide": pl.Boolean, "alternate_location_involved": pl.Boolean,
              "bond_classification": pl.String}
    for side in ("from", "to"):
        for name in ("chain", "residue_index", "author_residue", "insertion_code", "residue", "atom", "altloc"):
            schema[f"{name}.{side}"] = pl.Int64 if name in {"residue_index", "author_residue"} else pl.String
    pl.DataFrame(rows, schema=schema).write_csv(output, separator="\t")


@single_threaded
def export_md_trajectory(trajectory, reference, components_path, potential, output, *,
                         mhc_class, pbc_whole=False, cutoff=5.0, docking=False,
                         organism="human", threads=1):
    """Export exact existing trajectory readouts into a new directory.

    The caller declares PBC preparation and MHC class; neither is inferred. Component
    JSON maps disjoint names to chain-ID lists, including both groove chains for class II.
    Failed runs retain a failed provenance record and cannot be resumed. Source hashes
    are checked again before success, and generated coordinates have no predictor confidence.
    """
    threads = resolve_threads(threads)
    if not pbc_whole or mhc_class not in {"I", "II"}:
        raise ValueError("Explicit PBC-whole acknowledgement and MHC class I/II required")
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("Positive finite contact cutoff in angstroms required")
    sources = {name: _identity(path) for name, path in {
        "trajectory": trajectory, "reference": reference,
        "components": components_path, "potential": potential,
    }.items()}
    components = json.loads(Path(components_path).read_text())
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    record = {"state": "running", "tcren_version": __version__, "sources": sources,
              "mhc_class": mhc_class, "mhc_class_source": "explicit caller declaration",
              "components": components, "cpu_budget": threads, "docking": docking,
              "organism": organism, "contact_cutoff_A": cutoff,
              "command": shlex.join(sys.argv),
              "pbc": "caller declares whole molecules and coherent continuous complex images",
              "units": {"coordinates": "angstrom", "time": "ps", "potential": "kJ/mol"},
              "potential_scope": "whole system; not binding energy",
              "predictor_confidence": "not inherited by evolved coordinates"}
    provenance = output / "provenance.json"
    _json(provenance, record)
    try:
        result = inspect_md_trajectory(trajectory, reference, components, potential,
                                      contact_cutoff=cutoff, docking=docking,
                                      organism=organism, threads=threads)
        schemas = {
            "backbone_atoms": {"atom_index": pl.Int64, "chain.id": pl.String,
                               "residue.index": pl.Int64, "author_residue": pl.Int64,
                               "insertion_code": pl.String, "atom": pl.String,
                               "altloc": pl.String, "component": pl.String},
            "rmsd": {"frame_index": pl.Int64, "time_ps": pl.Float64,
                     "fit_component": pl.String, "measured_component": pl.String,
                     "atoms": pl.Int64, "rmsd_A": pl.Float64},
            "contact_counts": {"frame_index": pl.Int64, "time_ps": pl.Float64,
                               "chain.id.from": pl.String, "chain.id.to": pl.String,
                               "cutoff_A": pl.Float64, "residue_pairs": pl.Int64},
            "contact_occupancy": {"chain.id.from": pl.String, "residue.index.from": pl.Int64,
                                  "chain.id.to": pl.String, "residue.index.to": pl.Int64,
                                  "cutoff_A": pl.Float64, "frames_contacting": pl.Int64,
                                  "minimum_distance_A": pl.Float64,
                                  "maximum_contact_distance_A": pl.Float64,
                                  "frames": pl.Int64, "saved_frame_occupancy": pl.Float64},
            "docking": {"frame_index": pl.Int64, "time_ps": pl.Float64},
        }
        for name, schema in schemas.items():
            # Docking fields are owned by docking_angles/tcr_placement, not redefined here.
            table = pl.DataFrame(result[name]) if name == "docking" and result[name] else pl.DataFrame(result[name], schema=schema)
            table.write_csv(output / f"{name}.tsv", separator="\t")
        for name, key in {"first-frame": "first_frame_pdb", "last-frame": "last_frame_pdb",
                          "minimum-potential-frame": "selected_frame_pdb"}.items():
            path = output / f"{name}.pdb"
            path.write_text(result[key])
        _json(output / "saved-frame-times.json", result["saved_frame_times_ps"])
        _json(output / "minimum-potential-frame.json", result["minimum_potential"])
        if any(_identity(value["path"]) != value for value in sources.values()):
            raise ValueError("MD source changed during inspection")
        record.update(state="complete", frames=len(result["saved_frame_times_ps"]),
                      annotation_provenance=result["annotation_provenance"],
                      initial_image_status=result["initial_image_status"],
                      occupancy_scope="fraction of saved frames; observed residue pairs only",
                      outputs={p.name: _identity(p) for p in sorted(output.iterdir()) if p != provenance})
    except Exception as exc:
        record.update(state="failed", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        _json(provenance, record)
    return record


@app.command("inspect")
def inspect_command(
    trajectory: Path = typer.Option(..., exists=True, dir_okay=False, help="PBC-corrected multi-model PDB with TITLE t= ps"),
    reference: Path = typer.Option(..., exists=True, dir_okay=False, help="Corresponding backbone reference PDB"),
    components: Path = typer.Option(..., exists=True, dir_okay=False, help="JSON component name to chain-ID lists"),
    potential: Path = typer.Option(..., exists=True, dir_okay=False, help="Two-column time ps / whole-system potential kJ/mol XVG"),
    output: Path = typer.Option(..., "-o", "--output", help="New output directory; never resumes/overwrites"),
    mhc_class: MhcClass = typer.Option(..., "--mhc-class", help="Explicit class; include both groove chains for II"),
    pbc_whole: bool = typer.Option(False, "--pbc-whole", help="Acknowledge whole molecules and coherent continuous periodic images"),
    cutoff: float = typer.Option(5.0, min=0.0, help="Heavy-atom residue-contact cutoff in angstroms"),
    docking: bool = typer.Option(False, help="Annotate reference once and export frame-specific docking"),
    organism: str = typer.Option("human"),
    threads: int = typer.Option(1, "-t", "--threads", help="CPU ceiling: 1 serial, 0 all allocated; annotation alone uses extra cores"),
):
    """Partner-relative RMSD, saved-frame minimum, contacts and provenance from existing APIs."""
    try:
        record = export_md_trajectory(trajectory, reference, components, potential, output,
                                      mhc_class=mhc_class.value, pbc_whole=pbc_whole,
                                      cutoff=cutoff, docking=docking, organism=organism, threads=threads)
    except (ValueError, OSError, RuntimeError) as exc:
        typer.echo(f"MD inspection failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"{record['frames']} frames exported; class {mhc_class.value}; {output}")


@app.command("quality")
def quality_command(
    structure: Path = typer.Option(..., "-s", "--structure", exists=True, dir_okay=False),
    output: Path = typer.Option(..., "-o", "--output", help="New output directory"),
    mhc_class: MhcClass = typer.Option(..., "--mhc-class"),
    cutoff: float = typer.Option(2.0, help="Strict heavy-atom proximity cutoff in angstroms"),
):
    """First-model coordinate diagnostics and named nonlocal atom pairs, not a clashscore.

    Retains alternate locations and possible disulfides. Covalent bond classification
    remains unknown; same-residue and adjacent same-chain pairs are excluded heuristically.
    """
    with native_threads():
        if not np.isfinite(cutoff) or cutoff <= 0:
            raise typer.BadParameter("Positive finite proximity cutoff in angstroms required")
        record = {"state": "running", "tcren_version": __version__,
                  "source": _identity(structure), "mhc_class": mhc_class.value,
                  "mhc_class_source": "explicit caller declaration", "cpu_budget": 1,
                  "cutoff_A": cutoff, "command": shlex.join(sys.argv)}
        try:
            output.mkdir(parents=True, exist_ok=False)
        except OSError as exc:
            raise typer.BadParameter(str(exc)) from exc
        provenance = output / "provenance.json"
        try:
            _pairs(structure, output / "close_atom_pairs.tsv", cutoff)
            _json(output / "coordinate_diagnostics.json", coordinate_diagnostics(structure))
            if _identity(structure) != record["source"]:
                raise ValueError("Coordinate source changed during inspection")
            record.update(state="complete", outputs={p.name: _identity(p) for p in sorted(output.iterdir())})
        except (ValueError, OSError, RuntimeError) as exc:
            record.update(state="failed", error=f"{type(exc).__name__}: {exc}")
            typer.echo(record["error"], err=True)
            raise typer.Exit(1) from exc
        finally:
            _json(provenance, record)
        typer.echo(f"Atom-level diagnostics exported; class {mhc_class.value}; {output}")


@app.command("compare-forces")
def compare_forces_command(
    reference_energy: Path = typer.Option(..., exists=True, dir_okay=False),
    reference_forces: Path = typer.Option(..., exists=True, dir_okay=False),
    energy: Path = typer.Option(..., exists=True, dir_okay=False),
    forces: Path = typer.Option(..., exists=True, dir_okay=False),
    atoms: int = typer.Option(..., min=1, help="System-group atom count, identical order in both reruns"),
    topology: Path = typer.Option(..., exists=True, dir_okay=False, help="Declared common rerun topology/TPR, recorded by hash"),
    coordinates: Path = typer.Option(..., exists=True, dir_okay=False, help="Declared common rerun coordinates, recorded by hash"),
    same_system: bool = typer.Option(False, "--same-system", help="Acknowledge identical topology, coordinates and atom order"),
    mhc_class: MhcClass = typer.Option(..., "--mhc-class"),
    output: Path = typer.Option(..., "-o", "--output", help="New output directory"),
):
    """Numerical energy/force fidelity of same-frame reruns, including PME-order comparisons.

    This does not run GROMACS or verify that the supplied XVGs originated from the
    declared common inputs. Preserve both rerun commands/MDPs alongside this export.
    """
    if not same_system:
        raise typer.BadParameter("--same-system acknowledgement required")
    with native_threads():
        sources = {name: _identity(path) for name, path in {
            "reference_energy": reference_energy, "reference_forces": reference_forces,
            "energy": energy, "forces": forces, "topology": topology,
            "coordinates": coordinates,
        }.items()}
        try:
            result = compare_same_frame_energy_forces(reference_energy, reference_forces,
                                                       energy, forces, atoms=atoms)
            if any(_identity(value["path"]) != value for value in sources.values()):
                raise ValueError("Rerun input changed during comparison")
            output.mkdir(parents=True, exist_ok=False)
        except (ValueError, OSError) as exc:
            raise typer.BadParameter(str(exc)) from exc
        _json(output / "comparison.json", result)
        _json(output / "provenance.json", {
            "state": "complete", "tcren_version": __version__, "cpu_budget": 1,
            "mhc_class": mhc_class.value, "mhc_class_source": "explicit caller declaration",
            "sources": sources, "command": shlex.join(sys.argv),
            "correspondence": "caller declares identical topology, coordinates and atom order",
            "units": {"time": "ps", "potential": "kJ/mol", "forces": "kJ/mol/nm"},
            "scope": "numerical fidelity; not binding energy or PME convergence",
            "outputs": {"comparison.json": _identity(output / "comparison.json")},
        })
        typer.echo(f"Same-frame numerical comparison exported; class {mhc_class.value}; {output}")
