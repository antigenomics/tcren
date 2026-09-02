"""Annotating a whole structure set in one pass, instead of one call per structure.

``classify_chains`` on its own spawns one ``arda`` / ``mmseqs easy-search`` per structure, each
building a temporary database from a handful of query sequences. Process startup dominates, and the
cost over a set the size of Native2026 is roughly an order of magnitude. Everything here sends every
chain of every structure in a **single call per organism** and slices the records back out
afterwards, which is why any path resolving to more than one structure goes through it.

This used to live in ``tcren.paper.helpers``, which is the paper's table-building module and sits at
the top of the stack. Five modules that are nowhere near the top -- ``footprint``, ``shuffle``,
``recent``, ``orient.superimpose`` and the descriptor dispatch -- were reaching up into it for the
batching, each with a function-local import to keep the cycle legal. The batching is infrastructure,
not paper code, so it lives here now and those imports point downwards. ``tcren.paper.helpers``
re-exports every name, so callers written against the old location keep working.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..structure.model import Structure

def iter_annotated_set(struct_dir: str | Path, on_error: str = "skip"):
    """Parse and chain-type every structure in a folder, yielding annotated structures.

    Annotation is **batched**: one arda call per organism covers every chain of every
    structure, because the per-call process startup dominates and per-structure annotation
    is an order of magnitude slower over a set the size of ``Native2026``. Global chain ids
    (``"<struct_idx>|<chain_id>"``) keep chains distinct across structures, and the records
    are sliced back per structure for ``classify_chains``.

    Args:
        struct_dir: Folder of PDB/mmCIF structures.
        on_error: ``"skip"`` (default) drops a structure that fails to parse or annotate;
            ``"raise"`` propagates.

    Yields:
        Chain-typed :class:`~tcren.structure.model.Structure` objects.
    """
    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..structure import parse_structure, structure_id_from_path, structure_paths

    structures: list[Structure] = []
    for path in structure_paths(Path(struct_dir)):
        # id resolved from the filename (handles "<id>.pdb(.gz)" and "<id>_renumbered.cif").
        try:
            structures.append(parse_structure(path, pdb_id=structure_id_from_path(path)))
        except Exception:
            if on_error == "raise":
                raise

    records_by_struct = _batch_annotate(structures, _import_arda())
    for idx, s in enumerate(structures):
        try:
            classify_chains(s, organism="human", autodetect_species=True,
                            precomputed_records=records_by_struct[idx])
        except Exception:
            if on_error == "raise":
                raise
            continue
        yield s
def annotate_structure_set(
    struct_dir: str | Path, on_error: str = "skip", count_atoms: bool = False
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Run the tcren pipeline over a folder of PDBs → ``(contacts, markup)`` tables.

    Replaces the legacy mir batch annotation. ``contacts`` is the stacked TCR↔peptide
    :func:`contact_table`; ``markup`` is one row per structure with the CDR3α/CDR3β/peptide
    sequences + species (the inputs to non-redundancy clustering and the benchmarks).
    Species is auto-detected per structure by alignment score (human vs mouse). All chains
    across the whole folder are annotated in a single mmseqs call per organism (the
    per-call process overhead dominates, so dataset-level batching is far faster than
    per-structure annotation).

    When ``count_atoms`` is set, each contact row carries an ``n_atom_contacts``
    heavy-atom-pair count (needed for atomic-weighted scoring).
    """
    # `contact_table` still lives in `paper.helpers`, which re-exports this module, so importing it
    # at the top would close the cycle this module was split out to open. Local, as everywhere else.
    from ..paper.helpers import contact_table  # noqa: PLC0415

    contacts, markup = [], []
    for s in iter_annotated_set(struct_dir, on_error=on_error):
        pdb_id = s.pdb_id
        try:
            ct = contact_table(s, count_atoms=count_atoms)
            if ct.height:
                contacts.append(ct)

            def _region_seq(chain_type, region):
                for c in s.chains:
                    if c.chain_type == chain_type:
                        for r in c.regions:
                            if r.region_type == region:
                                return r.sequence
                return None

            peptide = next((c.sequence() for c in s.chains if c.chain_type == "PEPTIDE"), None)
            markup.append({
                "pdb.id": pdb_id,
                "cdr3a": _region_seq("TRA", "CDR3"),
                "cdr3b": _region_seq("TRB", "CDR3"),
                "peptide": peptide,
                "species": s.complex_species,
            })
        except Exception:
            if on_error == "raise":
                raise
    contacts_df = pl.concat(contacts) if contacts else pl.DataFrame()
    markup_df = pl.DataFrame(markup) if markup else pl.DataFrame()
    return contacts_df, markup_df
def mhc_annotation(
    struct_dir, ids=None, organism: str = "human", on_error: str = "skip"
) -> pl.DataFrame:
    """Per-structure MHC allele + class for a folder (tcren mapper) — fully batched.

    Replaces the legacy ``PDB_MHC_annotation`` table. ``ids`` restricts to those PDB ids.
    Every chain is TCR-typed in one batched arda call (so MHC candidates can be found),
    then every candidate MHC chain across the whole folder is searched against the MHC
    reference in a **single** mmseqs ``easy_search`` (mmseqs parallelises internally — no
    Python process/thread pool, which would either deadlock on fork or re-pay the fixed
    mmseqs startup cost per structure). Returns ``pdb.id``, ``mhc.class``, ``mhc.allele``,
    ``status``.
    """
    import tempfile

    import arda.mmseqs as mmseqs

    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..mhc import reference
    from ..mhc.mapper import MhcCall, _best_hits, _candidate_chains, _reconcile_class
    from ..structure import parse_structure, structure_id_from_path, structure_paths

    struct_dir = Path(struct_dir)
    paths = structure_paths(struct_dir)
    if ids is not None:
        keep = set(ids)
        paths = [p for p in paths if structure_id_from_path(p) in keep]

    structures: list[Structure] = []
    for path in paths:
        pdb_id = structure_id_from_path(path)
        try:
            structures.append(parse_structure(path, pdb_id=pdb_id))
        except Exception:
            if on_error == "raise":
                raise

    # 1) Batched TCR chain-typing so the non-receptor / non-peptide MHC candidates are known.
    records_by_struct = _batch_annotate(structures, _import_arda())
    for idx, s in enumerate(structures):
        try:
            classify_chains(s, organism=organism, autodetect_species=True,
                            precomputed_records=records_by_struct[idx])
        except Exception:
            if on_error == "raise":
                raise

    # 2) One mmseqs search over every candidate MHC chain across all structures. Global ids
    #    "<struct_idx>|<chain_id>" keep chains unique; hits are sliced back per structure.
    flat = [
        (idx, c.chain_id, c.sequence())
        for idx, s in enumerate(structures)
        for c in _candidate_chains(s)
        if c.sequence()
    ]
    best: dict[str, dict] = {}
    if flat:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_fa = tmp / "query.fasta"
            with query_fa.open("w") as fh:
                for idx, cid, seq in flat:
                    fh.write(f">{idx}|{cid}\n{seq}\n")
            out_tsv = tmp / "hits.tsv"
            mmseqs.easy_search(query_fa, reference.reference_fasta(), out_tsv,
                               tmp / "mmseqs_tmp", search_type=1, sensitivity=5.7, max_seqs=50)
            best = _best_hits(out_tsv)

    # 3) Build MhcCalls per structure from the sliced hits, reconcile class, summarise.
    rows = []
    for idx, s in enumerate(structures):
        calls: list[MhcCall] = []
        for c in _candidate_chains(s):
            hit = best.get(f"{idx}|{c.chain_id}")
            if hit is None:
                continue
            meta = reference.parse_header(hit["target"])
            calls.append(MhcCall(
                chain_id=c.chain_id, chain_role=meta["chain_role"],
                mhc_class=meta["mhc_class"], allele=meta["allele"], locus=meta["locus"],
                species=meta["species"], identity=float(hit["pident"]), bits=float(hit["bits"]),
                qstart=int(hit["qstart"]), qend=int(hit["qend"]),
                tstart=int(hit["tstart"]), tend=int(hit["tend"]), cigar=hit["cigar"],
            ))
        _reconcile_class(calls)
        mhca = next((c for c in calls if c.chain_role == "MHCa"), None)
        if any(c.chain_role == "MHCb" for c in calls):
            mhc_class = "MHCII"
        elif mhca:
            mhc_class = "MHCI"
        else:
            mhc_class = None
        rows.append({
            "pdb.id": s.pdb_id, "mhc.class": mhc_class,
            "mhc.allele": mhca.allele if mhca else None, "status": "ok",
        })
    return pl.DataFrame(rows)
def annotate_batch(
    structures, arda=None, organisms=("human", "mouse"), threads: int = 0
) -> list[dict[str, dict[str, dict]]]:
    """Annotate every chain of every structure with one mmseqs call per organism.

    Public since 2.3.0. Batching matters: arda/mmseqs costs seconds per call, so annotating a
    1,000-structure cohort one at a time is minutes of pure index rebuild. Every benchmark that
    scores a cohort needs this, which is why four downstream scripts were reaching into the private
    name.

    ``threads`` caps the mmseqs thread count for THIS call. It matters whenever batches are
    annotated concurrently: ``arda.annotate_sequences`` does not forward a thread count, so
    ``annotate_records`` falls back to its ``threads=0`` default, which mmseqs reads as *all cores*.
    Twelve concurrent batches on a 16-core machine then ask for 192 threads and spend their time
    context-switching instead of searching. Pass ``max(1, cpu_count // concurrency)``. ``0`` keeps
    the all-cores default, which is right for a single non-concurrent call.

    Returns ``records[struct_idx][organism][chain_id]`` — the per-structure slices fed to
    :func:`~tcren.annotation.classify_chains` as ``precomputed_records``.

    ``arda`` defaults to the lazily imported backend; pass an instance to reuse one mmseqs handle
    across a large batch.
    """
    out: list[dict[str, dict[str, dict]]] = [
        {org: {} for org in organisms} for _ in structures
    ]
    flat = [
        (idx, c.chain_id, c.sequence())
        for idx, s in enumerate(structures)
        for c in s.chains
        if c.sequence()
    ]
    if not flat:
        return out
    if arda is None:  # match the single-structure API: resolve the backend here
        from ..annotation.arda_adapter import _import_arda
        arda = _import_arda()
    pairs = [(f"{idx}|{cid}", seq) for idx, cid, seq in flat]
    for org in organisms:
        records = _annotate(arda, pairs, org, threads)
        for (idx, cid, _seq), rec in zip(flat, records):
            out[idx][org][cid] = rec
    return out
def _annotate(arda, pairs, organism: str, threads: int):
    """``arda.annotate_sequences`` with a thread cap, falling back when arda cannot take one.

    The public adapter drops ``threads`` on the floor, so reach one level down to
    ``arda.annotate.mapper.annotate_records`` when a cap is asked for. Guarded: an arda without that
    module still works, just uncapped.
    """
    if threads > 0:
        try:
            from arda.annotate.mapper import annotate_records
            return annotate_records(pairs, organism=organism, seqtype="aa", threads=threads)
        except Exception:  # noqa: BLE001 - any arda without the private path keeps the default
            pass
    return arda.annotate_sequences(pairs, seqtype="aa", organism=organism)
def iter_typed(structures: Path, organism: str = "human"):
    """Yield chain-typed structures, annotating a whole set in one batch.

    ``classify_chains`` per structure spawns one ``mmseqs easy-search`` per structure, each
    building a temporary database from a handful of query sequences; the process startup dominates
    and the cost is roughly an order of magnitude. :func:`iter_annotated_set` sends
    every chain of every structure in a single call per organism, which is why anything resolving
    to more than one structure -- a directory, a glob, a manifest -- goes through it. A single
    file has nothing to batch and takes the direct path.
    """
    from ..structure import iter_structures, parse_structure, structure_paths
    try:
        many = len(structure_paths(structures)) > 1
    except Exception:  # .tar.gz and anything structure_paths cannot enumerate
        many = False
    from ..annotation import classify_chains

    if many:
        yield from iter_annotated_set(structures)
        return
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        yield s


#: Historical private alias; `annotate_batch` is the name to use.
_batch_annotate = annotate_batch
