"""Memory-bounded raw-coordinate CLI for immutable contact penalties."""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import time

import polars as pl
import typer

from .observations import contact_observations
from .plausibility import comparison_command, fit_command, score_contact_penalties
from ..runtime import resolve_threads, single_threaded


@single_threaded
def penalty_export(sources, resource_path, output, *, chunk=32, threads=1, potentials=("tcren2", "karnaukhov2022"), geometry_only=False, calibration_path=None, normalizer_path=None, surface_context=False, surface_ownership="deepest_sphere", atom_count_path=None):
    """Parse one bounded batch; annotate/search once per batch; release its atoms."""
    from ..annotation.batch import annotate_batch
    from ..annotation import classify_chains
    from ..mhc import annotate_mhc_batch
    from ..structure import parse_structure, structure_paths
    from ..structure.io import structure_stem, resolve_sources
    from ..contactmap import ContactMap
    from ..pipeline import _BUNDLED_POTENTIALS
    from ..exports import _sha256, residue_identity
    from ..provenance import annotation_indexes
    if surface_ownership not in {"deepest_sphere", "accessible_boundary"}:
        raise ValueError("unknown surface ownership method")
    if surface_ownership != "deepest_sphere" and not surface_context:
        raise ValueError("surface ownership requires --surface-context")
    if chunk < 1:
        raise ValueError("positive batch size required")
    threads = resolve_threads(threads)
    paths = [p for source in resolve_sources(sources) for p in structure_paths(source)]
    ids = [structure_stem(p) for p in paths]
    if not paths or len(set(ids)) != len(ids):
        raise ValueError("nonempty input set with globally unique complete stems required")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    resource_path = Path(resource_path)
    resource = json.loads(resource_path.read_text())
    atom_count = json.loads(Path(atom_count_path).read_text()) if atom_count_path else None
    if atom_count and geometry_only:
        raise ValueError("atom-count layer requires score export")
    calibration = json.loads(Path(calibration_path).read_text()) if calibration_path else None
    if calibration and calibration["marginal_resource_sha256"] != _sha256(resource_path):
        raise ValueError("calibration belongs to a different marginal resource")
    normalizer = json.loads(Path(normalizer_path).read_text()) if normalizer_path else None
    if normalizer:
        from .context_normalization import validate_normalizer
        validate_normalizer(normalizer, _sha256(resource_path))
        if geometry_only:
            raise ValueError("normalizer requires score export")
    if surface_context and geometry_only:
        raise ValueError("surface context requires score export")
    if surface_context:
        (output/"surface").mkdir()
    tables = {} if geometry_only else {name: _BUNDLED_POTENTIALS[name]() for name in potentials}
    for name, pot in tables.items():
        pot.matrix.write_csv(output/(name+".potential.tsv"), separator="\t")
    ledger, summaries, start = [], [], time.monotonic()
    lineage = {}
    first_edge = first_residue = True
    with gzip.open(output/("geometry.tsv.gz" if geometry_only else "contacts.tsv.gz"), "wt") as edges, gzip.open(output/"residues.tsv.gz", "wt") as residues, (output/"annotation.jsonl").open("w") as annotation:
        for offset in range(0, len(paths), chunk):
            parsed, surface_atoms, surface_edges = [], [], []
            for path, pid in zip(paths[offset:offset+chunk], ids[offset:offset+chunk]):
                row = dict(id=pid, source=str(path), sha256=_sha256(path), status="error", error="")
                ledger.append(row)
                try:
                    parsed.append((parse_structure(path, pdb_id=pid, keep_modified=True), row))
                except Exception as exc:
                    row["error"] = f"parse: {type(exc).__name__}: {exc}"
            try:
                records = annotate_batch([s for s, _ in parsed], threads=threads)
                typed = []
                for (structure, row), record in zip(parsed, records, strict=True):
                    try:
                        classify_chains(structure, organism="human", autodetect_species=True, precomputed_records=record)
                        typed.append((structure, row))
                    except Exception as exc:
                        row["error"] = f"typing: {type(exc).__name__}: {exc}"
                annotate_mhc_batch([s for s, _ in typed], threads=threads)
            except Exception as exc:
                for _, row in parsed:
                    if not row["error"]:
                        row["error"] = f"annotation: {type(exc).__name__}: {exc}"
                typed = []
            for structure, row in typed:
                try:
                    observations = contact_observations(structure, interface="tcr_peptide", cutoff=5.)
                    if normalizer:
                        from .context_normalization import normalized_edges, apply_normalized
                        current_class = "II" if any(c.chain_supertype == "MHCII" for c in structure.chains) else "I"
                        if current_class != normalizer["mhc_class"]:
                            raise ValueError("coordinate MHC class differs from normalizer")
                    if calibration or normalizer:
                        from .weighted import multiplicity_edges
                        opportunities = multiplicity_edges(observations, resource)
                    if normalizer:
                        increments = normalized_edges(opportunities, normalizer)
                    if calibration:
                        from .weighted import calibrated_weights, apply_weighted
                        weights = calibrated_weights(opportunities, calibration)
                    contacts = None if geometry_only else ContactMap.from_structure(structure, cutoff=5.).interface("tcr_peptide", tcr_regions="all")
                    identity = residue_identity(structure)
                    residues.write(identity.write_csv(separator="\t", include_header=first_residue))
                    first_residue = False
                    if geometry_only:
                        from .mechanism import geometry_diagnostics
                        diagnosis = geometry_diagnostics(observations, resource)
                        edges.write(diagnosis.write_csv(separator="\t", include_header=first_edge))
                        first_edge = False
                        summaries.append(dict(id=structure.pdb_id, n_cdr_edges=diagnosis.height))
                    if surface_context:
                        from .exposure import contact_surface
                        from .surface_reference import surface_reference_score
                        atom_surface, edge_surface = contact_surface(structure, ownership=surface_ownership)
                        surface_atoms.append(atom_surface)
                        surface_edges.append(edge_surface)
                    for pot in tables.values():
                        table, summary = score_contact_penalties(contacts, observations, pot, resource)
                        if atom_count:
                            from .atom_counts import score_atom_counts
                            current_class = "II" if any(c.chain_supertype == "MHCII" for c in structure.chains) else "I"
                            table, atom_summary = score_atom_counts(table, observations, atom_count, summary["original"], mhc_class=current_class)
                            summary.update(atom_summary)
                        if calibration:
                            table, summary = apply_weighted(table, weights, summary)
                        if normalizer:
                            table, summary = apply_normalized(table, increments, summary)
                        if surface_context:
                            summary.update(surface_reference_score(atom_surface, edge_surface, pot))
                        edges.write(table.write_csv(separator="\t", include_header=first_edge))
                        first_edge = False
                        summaries.append(dict(id=structure.pdb_id, mhc_class="II" if any(c.chain_supertype=="MHCII" for c in structure.chains) else "I", **summary))
                    lineage = annotation_indexes([structure], previous=lineage)
                    annotation.write(json.dumps(dict(id=structure.pdb_id, searched_indexes=annotation_indexes([structure])), default=str)+"\n")
                    row["status"] = "ok"
                except Exception as exc:
                    row["error"] = f"score: {type(exc).__name__}: {exc}"
            if surface_context and surface_atoms:
                pl.concat(surface_atoms, how="diagonal_relaxed").write_parquet(output/"surface"/f"atoms-{offset:06d}.parquet")
                pl.concat(surface_edges, how="diagonal_relaxed").write_parquet(output/"surface"/f"edges-{offset:06d}.parquet")
            print(json.dumps(dict(processed=min(offset+chunk, len(paths)), total=len(paths), elapsed_s=time.monotonic()-start)), flush=True)
    pl.DataFrame(ledger).write_csv(output/"accounting.tsv", separator="\t")
    if summaries:
        pl.DataFrame(summaries, infer_schema_length=None).write_csv(output/("geometry_models.tsv" if geometry_only else "scores.tsv"), separator="\t")
    provenance = dict(resource_sha256=_sha256(resource_path), resource_version=resource["version"],
        atom_count_sha256=_sha256(atom_count_path) if atom_count else None,
        calibration_sha256=_sha256(calibration_path) if calibration else None,
        normalizer_sha256=_sha256(normalizer_path) if normalizer else None,
        producer_sha256=_sha256(Path(__file__)), source_count=len(paths), chunk_size=chunk, threads=threads,
        elapsed_s=time.monotonic()-start, potential_sha256={name: _sha256(output/(name+".potential.tsv")) for name in tables},
        scope=("Fresh all-atom CDR:peptide chemistry and exact frozen marginal lineage; no score recalculation" if geometry_only else "Legacy all-TCR peptide energy; increments on CDR1/2/3 peptide contacts only; no pMHC or framework penalty"),
        surface_context=(dict(version="exposed-surface-reference-1", n_points=100, probe_A=1.4, sphere="native Fibonacci-y", opportunity="isolated CDR and intact-pMHC peptide exposed area; finite directed support", score="sum binder-allocated area * (V + logZ); A^2 times table scale", ownership=surface_ownership, ownership_definition=("minimum signed expanded-sphere gap; exact ties split equally" if surface_ownership == "deepest_sphere" else "nearest sampled isolated opposite boundary; exact nearest atom-owner ties split equally"), partitions="geometry-derived per structure, no labels or fitted weights", missing="omitted and explicitly reported area; no supported area gives empty sum 0") if surface_context else None),
        missing_matrix="Legacy unknown-entry omission retained; missing entries and partial totals explicitly counted",
        no_scoring_time_fitting=True, searched_indexes=lineage)
    (output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    if not summaries:
        raise ValueError("all input structures failed; see accounting.tsv")
    return ledger


def penalty_command(
    structures: list[str] = typer.Option(..., "-s", "--structures"),
    resource: Path = typer.Option(..., "--resource", exists=True, dir_okay=False),
    out: Path = typer.Option(..., "-o", "--out"),
    chunk: int = typer.Option(32, "--chunk", min=1),
    threads: int = typer.Option(1, "-t", "--threads"),
    calibration: Path | None = typer.Option(None, "--calibration", exists=True, dir_okay=False),
    normalizer: Path | None = typer.Option(None, "--normalizer", exists=True, dir_okay=False),
    surface_context: bool = typer.Option(False, "--surface-context"),
    surface_ownership: str = typer.Option("deepest_sphere", "--surface-ownership"),
    atom_count_resource: Path | None = typer.Option(None, "--atom-count-resource", exists=True, dir_okay=False),
):
    """Apply the same frozen geometry weight to four declared penalty arms."""
    penalty_export(structures, resource, out, chunk=chunk, threads=threads, calibration_path=calibration, normalizer_path=normalizer, surface_context=surface_context, surface_ownership=surface_ownership, atom_count_path=atom_count_resource)


def geometry_audit_command(
    structures: list[str] = typer.Option(..., "-s", "--structures"),
    resource: Path = typer.Option(..., "--resource", exists=True, dir_okay=False),
    out: Path = typer.Option(..., "-o", "--out"),
    chunk: int = typer.Option(32, "--chunk", min=1),
    threads: int = typer.Option(1, "-t", "--threads"),
):
    """Trace worst-tail geometry and all-atom chemistry without recalculating scores."""
    penalty_export(structures, resource, out, chunk=chunk, threads=threads, geometry_only=True)


def register(app, panel):
    """Keep the existing CLI entrypoint small."""
    app.command("contact-prior-fit", rich_help_panel=panel)(fit_command)
    from .atom_counts import fit_command as atom_count_fit
    app.command("contact-atom-count-fit", rich_help_panel=panel)(atom_count_fit)
    app.command("contact-penalties", rich_help_panel=panel)(penalty_command)
    app.command("contact-compare", rich_help_panel=panel)(comparison_command)
    app.command("contact-geometry-audit", rich_help_panel=panel)(geometry_audit_command)
    from .mechanism_summary import mechanism_command
    app.command("contact-mechanism", rich_help_panel=panel)(mechanism_command)
    from .context_normalization import normalizer_command
    app.command("contact-normalizer-fit", rich_help_panel=panel)(normalizer_command)
    from .weighted import calibration_command
    app.command("contact-weight-calibrate", rich_help_panel=panel)(calibration_command)
    from .weighted_summary import weighted_summary_command
    app.command("contact-weight-summary", rich_help_panel=panel)(weighted_summary_command)
