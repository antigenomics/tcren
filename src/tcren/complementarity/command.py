"""Bounded structure-to-dashboard command; legacy output is explicitly opt-in."""
from __future__ import annotations

import hashlib
import itertools
import json
import tarfile
import tempfile
from math import isfinite
from pathlib import Path

import typer

from .dashboard import write_dashboard
from .prepare import fixed_component, ranking_frame, reference_frame
from .scene import build_scene


def _archive_structures(source, scratch, identities=None):
    from ..structure import is_structure_file, parse_structure

    found = set()
    with tarfile.open(source, mode="r|*") as archive:
        for member in archive:
            path = Path(member.name)
            if not member.isfile() or path.name.startswith("._") or not is_structure_file(path):
                continue
            if identities is not None and member.name not in identities:
                continue
            if identities is not None:
                found.add(member.name)
            suffix = "".join(path.suffixes)
            selected = identities[member.name] if identities is not None else [(member.name[:-len(suffix)], None)]
            with archive.extractfile(member) as incoming:
                with tempfile.NamedTemporaryFile(suffix=suffix, dir=scratch) as temporary:
                    while chunk := incoming.read(1024 * 1024):
                        temporary.write(chunk)
                    temporary.flush()
                    digest = hashlib.sha256(Path(temporary.name).read_bytes()).hexdigest()
                    for identity, row in selected:
                        structure = parse_structure(temporary.name, pdb_id=identity)
                        structure.annotation_provenance["visualization_input"] = dict(
                            member=member.name, sha256=digest)
                        if identities is not None:
                            structure.annotation_provenance["visualization_input"].update(archive=str(source), metadata_row=row)
                        yield structure
    if identities is not None and (missing := set(identities) - found):
        raise ValueError(f"metadata archive members missing in {source}: {', '.join(sorted(missing)[:5])}")


def _structures(source, scratch, metadata=None):
    from ..structure import iter_structures, parse_structure

    if metadata is not None:
        archives, handled = {}, set()
        for item in metadata:
            if item.member:
                archives.setdefault(item.path, {}).setdefault(item.member, []).append((item.identity, item.row))
        for item in metadata:
            if item.member:
                if item.path not in handled:
                    handled.add(item.path)
                    yield from _archive_structures(item.path, scratch, archives[item.path])
            else:
                structure = parse_structure(item.path, pdb_id=item.identity)
                structure.annotation_provenance["visualization_input"] = dict(
                    path=str(item.path), metadata_row=item.row, sha256=hashlib.sha256(item.path.read_bytes()).hexdigest())
                yield structure
    elif str(source).endswith((".tar.gz", ".tgz", ".tar")):
        yield from _archive_structures(source, scratch)
    else:
        from ..structure import structure_paths

        paths = structure_paths(source)
        for path, (_, structure) in zip(paths, iter_structures(source, full_ids=True,
                                                              importer=parse_structure)):
            structure.annotation_provenance["visualization_input"] = dict(
                path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            yield structure


def complementarity_command(
    structures: Path = typer.Option(None, "--structures", "-s"),
    output: Path = typer.Option(..., "--output", "-o"),
    task: str = typer.Option("tcr", help="tcr: fixed pMHC; epitope: fixed receptor"),
    origin: str = typer.Option("center-of-mass", help="apex or center-of-mass"),
    cutoff: float = typer.Option(5., min=0.1, max=10.),
    threads: int = typer.Option(0, "--threads", "-t", min=0),
    legacy: bool = typer.Option(False, "--legacy", help="Use unchanged historical SVG maps"),
    scratch: Path = typer.Option(None, help="Existing directory for one temporary archive member"),
    exclude_id: list[str] = typer.Option(None, "--exclude-id", help="Explicit source ID to omit; repeatable and recorded"),
    metadata: Path = typer.Option(None, help="TSV[.gz]/CSV[.gz]/Parquet with explicit structure identities and paths"),
    path_column: str = typer.Option("structure_path", help="Metadata path column; relative paths resolve beside metadata"),
    id_column: str = typer.Option("id", help="Metadata structure identity column"),
    member_column: str = typer.Option(None, help="Metadata archive-member column, when paths point to tar files"),
):
    """Export shared-frame side/top contact scenes and a local overlay dashboard."""
    from ..annotation import classify_chains
    from ..annotation.batch import annotate_batch
    from ..mhc import annotate_mhc_batch

    if task not in {"tcr", "epitope"} or origin not in {"apex", "center-of-mass"}:
        raise typer.BadParameter("task must be tcr/epitope and origin apex/center-of-mass")
    if not isfinite(cutoff):
        raise typer.BadParameter("cutoff must be finite")
    if output.exists():
        raise typer.BadParameter("output must be a new directory")
    if scratch is not None and not scratch.is_dir():
        raise typer.BadParameter("scratch must be an existing directory")
    if legacy and not 3 <= cutoff <= 6:
        raise typer.BadParameter("legacy cutoff must be within [3, 6] Angstroms")
    if (structures is None) == (metadata is None):
        raise typer.BadParameter("provide exactly one of --structures or --metadata")
    source_rows, metadata_provenance = None, None
    if metadata is not None:
        from .inputs import metadata_inputs

        try:
            source_rows = metadata_inputs(metadata, path_column=path_column, id_column=id_column,
                                          member_column=member_column)
        except (ValueError, OSError) as exc:
            raise typer.BadParameter(str(exc)) from exc
        metadata_provenance = dict(path=str(metadata.resolve()), sha256=hashlib.sha256(metadata.read_bytes()).hexdigest(),
                                   id_column=id_column, path_column=path_column, member_column=member_column,
                                   order="direct paths in table order; each archive streamed once at its first occurrence")
    references, errors, seen = {}, [], set()
    excluded, input_count = [], 0
    requested_exclusions = set(exclude_id or [])
    implementation = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in sorted(Path(__file__).parent.iterdir())
                      if p.suffix in {".py", ".html", ".js"}}

    def selected_sources():
        nonlocal input_count
        for structure in _structures(structures, scratch, source_rows):
            if metadata_provenance is not None:
                structure.annotation_provenance["visualization_metadata"] = metadata_provenance
            input_count += 1
            if structure.pdb_id in requested_exclusions:
                excluded.append(dict(id=structure.pdb_id, reason="explicit --exclude-id",
                                     input=structure.annotation_provenance.get("visualization_input")))
            else:
                yield structure

    def scenes():
        source = iter(selected_sources())
        while batch := list(itertools.islice(source, 32)):
            records = annotate_batch(batch, organisms=("human", "mouse"), threads=threads)
            for structure, record in zip(batch, records):
                classify_chains(structure, organism="human", autodetect_species=True,
                                precomputed_records=record)
            annotate_mhc_batch(batch, threads=threads)
            for structure in batch:
                try:
                    if structure.pdb_id in seen:
                        raise ValueError("duplicate structure ID")
                    seen.add(structure.pdb_id)
                    if legacy:
                        from ..project2d import contacts_table, project_structure, residue_markup_table
                        from ..viz.svg2d import render_complementarity_map

                        svg = render_complementarity_map(
                            residue_markup_table(structure, project_structure(structure)),
                            contacts=contacts_table(structure, threshold=cutoff))
                        (output / f"legacy-{len(seen)}.svg").write_text(svg, encoding="utf-8")
                    else:
                        types = ("MHCa", "MHCb", "PEPTIDE") if task == "tcr" else ("TRA", "TRB")
                        identity = sorted((c.chain_type, c.sequence()) for c in structure.by_type(*types))
                        key = hashlib.sha256(json.dumps(identity).encode()).hexdigest()
                        if key not in references:
                            frame = reference_frame(structure, origin=origin, task=task)
                            references[key] = fixed_component(structure, task), frame
                        reference, frame = references[key]
                        current = ranking_frame(structure, reference, frame, task=task)
                        scene = build_scene(structure, current, cutoff=cutoff)
                        scene["task"] = task
                        peptides = structure.by_type("PEPTIDE")
                        label = peptides[0].sequence() if task == "tcr" and peptides else "Fixed TCR"
                        scene["group_label"] = f"{label} · {frame.reference_id}"
                        yield scene
                except ValueError as exc:
                    errors.append(dict(id=structure.pdb_id, error=str(exc),
                                       input=structure.annotation_provenance.get("visualization_input"),
                                       chain_roles={c.chain_id: c.chain_type for c in structure.chains}))
            batch.clear()
            records.clear()

    if legacy:
        output.mkdir(parents=True)
        list(scenes())
    else:
        write_dashboard(scenes(), output)
    missing_exclusions = requested_exclusions - {row["id"] for row in excluded}
    if missing_exclusions:
        errors.append(dict(id=None, error="excluded IDs not found: " + ", ".join(sorted(missing_exclusions))))
    (output / "errors.json").write_text(json.dumps(errors, indent=2), encoding="utf-8")
    (output / "excluded.json").write_text(json.dumps(excluded, indent=2), encoding="utf-8")
    from importlib.metadata import version

    (output / "run.json").write_text(json.dumps(dict(
        version=version("tcren"), source=str(metadata if metadata is not None else structures),
        metadata=metadata_provenance, task=task, origin=origin,
        cutoff=cutoff, threads=threads, legacy=legacy, structures=len(seen),
        input_structures=input_count, excluded_structures=len(excluded),
        errors=len(errors), implementation_sha256=implementation), indent=2), encoding="utf-8")
    if not seen:
        raise typer.BadParameter("no structures found")
    typer.echo(f"Wrote {output}: {len(seen)} eligible structures, {len(excluded)} excluded, {len(errors)} errors")
    if errors:
        raise typer.Exit(1)
