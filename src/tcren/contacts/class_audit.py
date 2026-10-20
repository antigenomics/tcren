"""Fresh bounded raw-structure driver for label-free contact-class observations."""
from __future__ import annotations

import json
from pathlib import Path
import time

import polars as pl
import typer

from .class_features import features_for_structure, opportunity_rows
from .class_statistics import STRATUM, class_statistics
from ..runtime import parallel_map, resolve_threads, single_threaded


def _audit_chunk(structures):
    """Process a contiguous annotated chunk; labels never enter the feature kernel."""
    results = []
    for structure in structures:
        try:
            results.append((structure.pdb_id, features_for_structure(structure), ""))
        except Exception as exc:
            results.append((structure.pdb_id, None, f"features: {type(exc).__name__}: {exc}"))
    return results


def _manifest(path):
    required = ["id", "receptor_cluster", *STRATUM, "source", "sha256"]
    frame = pl.read_csv(path, separator="\t", schema_overrides={
        **{name: pl.String for name in required}, "label": pl.Int64})
    if not set([*required, "label"]) <= set(frame.columns) or not frame.height:
        raise ValueError("nonempty manifest requires id, label, receptor_cluster, stratum, source and sha256")
    if frame.select(pl.any_horizontal(pl.col([*required, "label"]).is_null()).any()).item():
        raise ValueError("manifest required metadata must be nonnull")
    if frame.select(pl.any_horizontal(pl.col(required).str.len_chars() == 0).any()).item():
        raise ValueError("manifest required metadata must be nonempty")
    if frame["id"].n_unique() != frame.height or not set(frame["label"].unique()) <= {0, 1}:
        raise ValueError("manifest ids must be globally unique and labels binary 0/1")
    if frame.filter(~pl.col("sha256").str.contains(r"^[0-9a-f]{64}$")).height:
        raise ValueError("manifest sha256 must be lowercase hexadecimal SHA-256")
    if set(frame["mhc_class"].unique()) - {"I", "II"}:
        raise ValueError("manifest mhc_class must be I or II")
    if frame.group_by(STRATUM).agg(pl.col("label").n_unique().alias("classes")).filter(pl.col("classes") != 2).height:
        raise ValueError("each exact case/stratum requires both binder and nonbinder labels")
    return frame


def _finish_provenance(output, ledger, provenance, start):
    pl.DataFrame(ledger, infer_schema_length=None).write_csv(output / "accounting.tsv", separator="\t")
    provenance.update(elapsed_s=time.monotonic() - start,
                      processed_count=sum(row["status"] == "ok" for row in ledger),
                      failed_count=sum(row["status"] == "error" for row in ledger),
                      pending_count=sum(row["status"] == "pending" for row in ledger))
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, default=str) + "\n")


@single_threaded
def class_audit(manifest, output, *, chunk=128, threads=0):
    """Rebuild complete class-audit evidence from pinned raw inputs, without reuse.

    Exact strata are processed separately. Only one raw batch is retained at a
    time, while compact counts/marginals survive until the global chemistry union
    is known. Raw inventory/pair partitions and final C/O/dictionary/statistics
    are requested audit evidence, never inputs to this run. Any source, parsing,
    annotation or feature failure prevents all statistics and leaves a ledger.
    """
    from ..annotation.batch import annotate_batch
    from ..annotation import classify_chains
    from ..mhc import annotate_mhc_batch
    from ..structure import parse_structure
    from ..exports import _sha256
    from ..provenance import annotation_indexes

    if isinstance(chunk, bool) or not isinstance(chunk, int) or chunk < 1:
        raise ValueError("positive integer batch size required")
    threads = resolve_threads(threads)
    manifest = Path(manifest)
    metadata = _manifest(manifest)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    metadata.write_parquet(output / "metadata.parquet", compression="zstd")
    start = time.monotonic()
    ledger = [dict(id=row["id"], source=row["source"], expected_sha256=row["sha256"], sha256="",
                   status="pending", stage="pending", error="", n_inventory_atoms=0,
                   n_physical_pairs=0) for row in metadata.iter_rows(named=True)]
    by_id = {row["id"]: row for row in ledger}
    provenance = dict(complete=False, manifest_sha256=_sha256(manifest), source_count=metadata.height,
                      chunk_size=chunk, threads=threads, cutoff_A=5.0,
                      producer_sha256={name: _sha256(Path(__file__).with_name(name + ".py"))
                                       for name in ("class_audit", "class_features", "class_statistics", "observations")},
                      searched_indexes={}, coordinate_policy="selected observed heavy atoms; no repair",
                      contact_scope="all TCR regions to peptide and annotated MHC helices; full inventory retained",
                      views=["site", "chemistry"], levels=["atom", "residue"],
                      opportunities="full selected endpoint availability; C0 retained; O0 omitted explicitly",
                      no_potentials=True, no_scores=True, no_fitted_coefficients=True)
    cases, chemistry_counts, binder_marginals, target_marginals, chemistry_dictionaries = [], [], [], [], []
    try:
        with (output / "annotation.jsonl").open("w") as annotation:
            for case_number, case_meta in enumerate(metadata.sort(STRATUM + ["id"]).partition_by(STRATUM, maintain_order=True), 1):
                folder = output / f"case{case_number:03d}"
                (folder / "inventory").mkdir(parents=True)
                (folder / "pairs").mkdir()
                case_meta.write_parquet(folder / "metadata.parquet", compression="zstd")
                motifs, bm, tm, dictionaries = [], [], [], []
                for batch_number, offset in enumerate(range(0, case_meta.height, chunk)):
                    parsed = []
                    for source in case_meta.slice(offset, chunk).iter_rows(named=True):
                        row = by_id[source["id"]]
                        row.update(stage="source", status="error")
                        try:
                            row["sha256"] = _sha256(Path(source["source"]))
                            if row["sha256"] != source["sha256"]:
                                raise ValueError("raw source SHA-256 mismatch")
                            row["stage"] = "parse"
                            parsed.append((parse_structure(source["source"], pdb_id=source["id"], keep_modified=True), row))
                        except Exception as exc:
                            row["error"] = f"{row['stage']}: {type(exc).__name__}: {exc}"
                    typed = []
                    if parsed:
                        try:
                            records = annotate_batch([s for s, _ in parsed], threads=threads)
                            for (structure, row), record in zip(parsed, records, strict=True):
                                row["stage"] = "typing"
                                try:
                                    classify_chains(structure, organism="human", autodetect_species=True, precomputed_records=record)
                                    typed.append((structure, row))
                                except Exception as exc:
                                    row["error"] = f"typing: {type(exc).__name__}: {exc}"
                            annotate_mhc_batch([s for s, _ in typed], threads=threads)
                        except Exception as exc:
                            for _, row in parsed:
                                if not row["error"]:
                                    row.update(stage="annotation", error=f"annotation: {type(exc).__name__}: {exc}")
                            typed = []
                    ready = []
                    for structure, row in typed:
                        try:
                            actual_class = "II" if any(c.chain_supertype == "MHCII" for c in structure.chains) else "I"
                            if actual_class != case_meta["mhc_class"][0]:
                                raise ValueError("annotated coordinate MHC class differs from manifest")
                            provenance["searched_indexes"] = annotation_indexes([structure], previous=provenance["searched_indexes"])
                            annotation.write(json.dumps(dict(id=structure.pdb_id, searched_indexes=annotation_indexes([structure])), default=str) + "\n")
                            ready.append((structure, row))
                        except Exception as exc:
                            row.update(stage="annotation", error=f"annotation: {type(exc).__name__}: {exc}")
                    count = min(threads, len(ready))
                    work = [[s for s, _ in ready[i * len(ready) // count:(i + 1) * len(ready) // count]] for i in range(count)]
                    inventory_tables, pair_tables = [], []
                    try:
                        groups = parallel_map(_audit_chunk, work, threads)
                        for (_, row), (pid, tables, error) in zip(
                                ready, (result for group in groups for result in group), strict=True):
                            if pid != row["id"]:
                                raise ValueError("worker result identity differs from submitted model")
                            row["stage"] = "features"
                            if error:
                                row["error"] = error
                                continue
                            inventory, pairs, c, b, t, dictionary = tables
                            inventory_tables.append(inventory)
                            pair_tables.append(pairs)
                            motifs.append(c); bm.append(b); tm.append(t); dictionaries.append(dictionary)
                            row.update(status="ok", stage="complete", n_inventory_atoms=inventory.height,
                                       n_physical_pairs=pairs.height)
                    except Exception as exc:
                        for _, row in ready:
                            row.update(status="error", stage="workers", error=f"workers: {type(exc).__name__}: {exc}")
                    if inventory_tables:
                        pl.concat(inventory_tables).write_parquet(folder / "inventory" / f"batch{batch_number:03d}.parquet", compression="zstd")
                        pl.concat(pair_tables).write_parquet(folder / "pairs" / f"batch{batch_number:03d}.parquet", compression="zstd")
                    del parsed, typed, ready, work, inventory_tables, pair_tables
                    if "groups" in locals():
                        del groups
                    if "structure" in locals():
                        del structure
                    if "_" in locals():
                        del _
                    if "records" in locals():
                        del records
                    if "tables" in locals():
                        del tables
                    if "inventory" in locals():
                        del inventory, pairs
                    print(json.dumps(dict(case_number=case_number, processed=sum(r["status"] == "ok" for r in ledger),
                                          total=metadata.height, elapsed_s=time.monotonic() - start)), flush=True)
                if motifs:
                    c, b, t = (pl.concat(items) for items in (motifs, bm, tm))
                    dictionary = pl.concat(dictionaries).unique()
                    chemistry = dictionary.filter(pl.col("view") == "chemistry")
                    chemistry_counts.append(c.join(chemistry.select("feature"), on="feature", how="inner", validate="m:1"))
                    binder_marginals.append(b); target_marginals.append(t); chemistry_dictionaries.append(chemistry)
                    cases.append((folder, case_meta, c, b, t, dictionary.filter(pl.col("view") == "site")))
                del motifs, bm, tm, dictionaries
        if any(row["status"] != "ok" for row in ledger):
            raise ValueError("one or more input models failed; statistics aborted; see accounting.tsv")
        chemistry_dictionary = pl.concat(chemistry_dictionaries).unique().sort("feature")
        global_rows = opportunity_rows(pl.concat(chemistry_counts), pl.concat(binder_marginals),
                                       pl.concat(target_marginals), chemistry_dictionary,
                                       metadata["id"]).rename({"pdb.id": "id"})
        chemistry_per, overall = class_statistics(global_rows, metadata)
        overall.join(chemistry_dictionary, on="feature", how="left", validate="m:1").write_parquet(output / "overall.parquet", compression="zstd")
        chemistry_dictionary.write_parquet(output / "chemistry_dictionary.parquet", compression="zstd")
        global_rows.write_parquet(output / "chemistry_model_features.parquet", compression="zstd")
        for folder, case_meta, c, b, t, site_dictionary in cases:
            site_counts = c.join(site_dictionary.select("feature"), on="feature", how="inner", validate="m:1")
            sites = opportunity_rows(site_counts, b, t, site_dictionary, case_meta["id"]).rename({"pdb.id": "id"})
            chemistry_rows = global_rows.join(case_meta.select("id"), on="id", how="inner", validate="m:1")
            rows = pl.concat([sites, chemistry_rows]).sort("id", "feature")
            dictionary = pl.concat([site_dictionary, chemistry_dictionary]).sort("feature")
            per = chemistry_per.join(case_meta.select(STRATUM).unique(), on=STRATUM, how="inner", validate="m:1")
            if sites.height:
                site_per, _ = class_statistics(sites, case_meta)
                per = pl.concat([per, site_per])
            rows.write_parquet(folder / "model_features.parquet", compression="zstd")
            dictionary.write_parquet(folder / "feature_dictionary.parquet", compression="zstd")
            per.join(dictionary, on="feature", how="left", validate="m:1").write_parquet(folder / "statistics.parquet", compression="zstd")
        provenance.update(complete=True, case_count=len(cases), chemistry_feature_count=chemistry_dictionary.height,
                          n_inventory_atoms=sum(r["n_inventory_atoms"] for r in ledger),
                          n_physical_pairs=sum(r["n_physical_pairs"] for r in ledger))
    except Exception as exc:
        provenance["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        _finish_provenance(output, ledger, provenance, start)
    return ledger


def audit_command(
    manifest: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    out: Path = typer.Option(..., "-o", "--out"),
    chunk: int = typer.Option(128, "--chunk", min=1),
    threads: int = typer.Option(0, "-t", "--threads", help="0: all allocated CPUs; 1: serial; N: explicit CPU budget"),
):
    """Fresh contact-class counts and opportunity-aware clustered contrasts."""
    class_audit(manifest, out, chunk=chunk, threads=threads)
