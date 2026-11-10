"""Bounded structure-to-dashboard command; legacy output is explicitly opt-in."""
from __future__ import annotations

import hashlib
import itertools
import json
import tarfile
import tempfile
from pathlib import Path

import typer

from .dashboard import write_dashboard
from .prepare import ranking_frame, reference_frame
from .scene import build_scene


def _structures(source, scratch):
    from ..structure import is_structure_file, iter_structures, parse_structure

    if str(source).endswith((".tar.gz", ".tgz", ".tar")):
        with tarfile.open(source, mode="r|*") as archive:
            for member in archive:
                path = Path(member.name)
                if not member.isfile() or path.name.startswith("._") or not is_structure_file(path):
                    continue
                suffix = "".join(path.suffixes)
                identity = member.name[:-len(suffix)]
                with archive.extractfile(member) as incoming:
                    with tempfile.NamedTemporaryFile(suffix=suffix, dir=scratch) as temporary:
                        while chunk := incoming.read(1024 * 1024):
                            temporary.write(chunk)
                        temporary.flush()
                        structure = parse_structure(temporary.name, pdb_id=identity)
                        structure.annotation_provenance["visualization_input"] = dict(
                            member=member.name,
                            sha256=hashlib.sha256(Path(temporary.name).read_bytes()).hexdigest())
                        yield structure
    else:
        from ..structure import structure_paths

        paths = structure_paths(source)
        for path, (_, structure) in zip(paths, iter_structures(source, full_ids=True,
                                                              importer=parse_structure)):
            structure.annotation_provenance["visualization_input"] = dict(
                path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            yield structure


def complementarity_command(
    structures: Path = typer.Option(..., "--structures", "-s"),
    output: Path = typer.Option(..., "--output", "-o"),
    task: str = typer.Option("tcr", help="tcr: fixed pMHC; epitope: fixed receptor"),
    origin: str = typer.Option("center-of-mass", help="apex or center-of-mass"),
    cutoff: float = typer.Option(5., min=0.1, max=10.),
    threads: int = typer.Option(0, "--threads", "-t", min=0),
    legacy: bool = typer.Option(False, "--legacy", help="Use unchanged historical SVG maps"),
    scratch: Path = typer.Option(None, help="Existing directory for one temporary archive member"),
):
    """Export shared-frame side/top contact scenes and a local overlay dashboard."""
    from ..annotation import classify_chains
    from ..annotation.batch import annotate_batch
    from ..mhc import annotate_mhc_batch

    if task not in {"tcr", "epitope"} or origin not in {"apex", "center-of-mass"}:
        raise typer.BadParameter("task must be tcr/epitope and origin apex/center-of-mass")
    if output.exists():
        raise typer.BadParameter("output must be a new directory")
    if scratch is not None and not scratch.is_dir():
        raise typer.BadParameter("scratch must be an existing directory")
    if legacy and not 3 <= cutoff <= 6:
        raise typer.BadParameter("legacy cutoff must be within [3, 6] Angstroms")
    references, errors, seen = {}, [], set()
    implementation = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in sorted(Path(__file__).parent.iterdir())
                      if p.suffix in {".py", ".html", ".js"}}

    def scenes():
        source = iter(_structures(structures, scratch))
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
                            references[key] = structure, frame
                        reference, frame = references[key]
                        current = ranking_frame(structure, reference, frame, task=task)
                        scene = build_scene(structure, current, cutoff=cutoff)
                        scene["task"] = task
                        peptides = structure.by_type("PEPTIDE")
                        label = peptides[0].sequence() if task == "tcr" and peptides else "Fixed TCR"
                        scene["group_label"] = f"{label} · {frame.reference_id}"
                        yield scene
                except ValueError as exc:
                    errors.append(dict(id=structure.pdb_id, error=str(exc)))

    if legacy:
        output.mkdir(parents=True)
        list(scenes())
    else:
        write_dashboard(scenes(), output)
    (output / "errors.json").write_text(json.dumps(errors, indent=2), encoding="utf-8")
    from importlib.metadata import version

    (output / "run.json").write_text(json.dumps(dict(
        version=version("tcren"), source=str(structures), task=task, origin=origin,
        cutoff=cutoff, threads=threads, legacy=legacy, structures=len(seen),
        errors=len(errors), implementation_sha256=implementation), indent=2), encoding="utf-8")
    if not seen:
        raise typer.BadParameter("no structures found")
    typer.echo(f"Wrote {output}: {len(seen)} structures, {len(errors)} errors")
    if errors:
        raise typer.Exit(1)
