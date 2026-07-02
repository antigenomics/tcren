"""Helpers for the Nat Comput Sci 2022 reproduction notebooks.

``contact_table`` replaces the legacy mir ``extract_contact_map`` (it returns the same
TCR↔peptide contact columns the R analyses consume, computed through the tcren pipeline).
``compare`` is the small regression utility behind ``07_compare_legacy.ipynb``.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from ..contactmap import ContactMap
from ..structure.model import Structure

# ContactMap.tcr_peptide() column -> R analysis column name.
_CONTACT_RENAME = {
    "pdb.id": "pdb.id",
    "chain.type.from": "chain.type.from",
    "region.type.from": "region.type.from",
    "residue.index.from": "residue.index.from",
    "residue.index.to": "residue.index.to",
    "pos.from": "pos.from",
    "pos.to": "pos.to",
    "residue.aa.from": "residue.aa.from",
    "residue.aa.to": "residue.aa.to",
}


def contact_table(
    structure: Structure, cutoff: float = 5.0, count_atoms: bool = False
) -> pl.DataFrame:
    """TCR↔peptide contact table for an annotated structure (the mir-replacement).

    The structure must already be chain-typed (``classify_chains``) and MHC-annotated
    (``annotate_mhc``). Returns the columns the R benchmarks use:
    ``pdb.id, chain.type.from, region.type.from, residue.index.from, residue.index.to,
    pos.from, pos.to, residue.aa.from, residue.aa.to``.

    When ``count_atoms`` is set, an extra ``n_atom_contacts`` column (the heavy-atom-pair
    count per residue pair) is carried through for atomic-weighted scoring. Default
    ``False`` keeps the schema byte-identical to the legacy output.
    """
    tp = ContactMap.from_structure(
        structure, cutoff=cutoff, count_atoms=count_atoms
    ).tcr_peptide()
    cols = list(_CONTACT_RENAME)
    if count_atoms:
        cols.append("n_atom_contacts")
    return tp.select(cols).unique()


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
    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..structure import parse_structure, structure_id_from_path, structure_paths

    struct_dir = Path(struct_dir)
    paths = structure_paths(struct_dir)
    structures: list[Structure] = []
    for path in paths:
        # id resolved from the filename (handles "<id>.pdb(.gz)" and "<id>_renumbered.cif").
        pdb_id = structure_id_from_path(path)
        try:
            structures.append(parse_structure(path, pdb_id=pdb_id))
        except Exception:
            if on_error == "raise":
                raise

    # One arda call per organism over every chain of every structure. Global ids
    # (``"<struct_idx>|<chain_id>"``) keep chains unique across structures; the records
    # are sliced back per structure and fed to classify_chains (no per-chain mmseqs).
    records_by_struct = _batch_annotate(structures, _import_arda())

    contacts, markup = [], []
    for idx, s in enumerate(structures):
        pdb_id = s.pdb_id
        try:
            classify_chains(s, organism="human", autodetect_species=True,
                            precomputed_records=records_by_struct[idx])
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


def _batch_annotate(
    structures, arda, organisms=("human", "mouse")
) -> list[dict[str, dict[str, dict]]]:
    """Annotate every chain of every structure with one mmseqs call per organism.

    Returns ``records[struct_idx][organism][chain_id]`` — the per-structure slices fed to
    :func:`~tcren.annotation.classify_chains` as ``precomputed_records``.
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
    pairs = [(f"{idx}|{cid}", seq) for idx, cid, seq in flat]
    for org in organisms:
        records = arda.annotate_sequences(pairs, seqtype="aa", organism=org)
        for (idx, cid, _seq), rec in zip(flat, records):
            out[idx][org][cid] = rec
    return out


def _read_any(path: str | Path) -> pl.DataFrame:
    """Read a CSV/TSV, transparently handling ``.gz`` and tab vs comma."""
    path = Path(path)
    name = path.name[:-3] if path.suffix == ".gz" else path.name
    sep = "\t" if name.endswith((".tsv", ".txt")) else ","
    return pl.read_csv(path, separator=sep, infer_schema_length=2000)


def compare(
    old_path: str | Path,
    new_path: str | Path,
    keys: list[str],
    value_cols: list[str] | None = None,
    tol: float = 1e-6,
) -> dict:
    """Compare two tables on ``keys`` and report row-set + max numeric differences.

    Returns ``{rows_old, rows_new, matched, only_old, only_new, max_abs_diff, status}``
    where ``status`` is ``"pass"`` when the key sets agree and every shared numeric column
    differs by ≤ ``tol``.
    """
    old, new = _read_any(old_path), _read_any(new_path)
    ko = set(map(tuple, old.select(keys).rows()))
    kn = set(map(tuple, new.select(keys).rows()))
    only_old, only_new = ko - kn, kn - ko

    max_abs = 0.0
    if value_cols is None:
        value_cols = [
            c for c in old.columns
            if c in new.columns and c not in keys and old[c].dtype.is_numeric()
        ]
    if value_cols and not only_old and not only_new:
        joined = old.join(new, on=keys, how="inner", suffix="__new")
        for c in value_cols:
            diff = (joined[c] - joined[f"{c}__new"]).abs().max()
            if diff is not None:
                max_abs = max(max_abs, float(diff))

    status = "pass" if not only_old and not only_new and max_abs <= tol else "FAIL"
    return {
        "rows_old": old.height, "rows_new": new.height,
        "matched": len(ko & kn), "only_old": len(only_old), "only_new": len(only_new),
        "max_abs_diff": max_abs, "status": status,
    }
