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


def contact_table(structure: Structure, cutoff: float = 5.0) -> pl.DataFrame:
    """TCR↔peptide contact table for an annotated structure (the mir-replacement).

    The structure must already be chain-typed (``classify_chains``) and MHC-annotated
    (``annotate_mhc``). Returns the columns the R benchmarks use:
    ``pdb.id, chain.type.from, region.type.from, residue.index.from, residue.index.to,
    pos.from, pos.to, residue.aa.from, residue.aa.to``.
    """
    tp = ContactMap.from_structure(structure, cutoff=cutoff).tcr_peptide()
    return tp.select(list(_CONTACT_RENAME)).unique()


def annotate_structure_set(
    struct_dir: str | Path, on_error: str = "skip"
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Run the tcren pipeline over a folder of PDBs → ``(contacts, markup)`` tables.

    Replaces the legacy mir batch annotation. ``contacts`` is the stacked TCR↔peptide
    :func:`contact_table`; ``markup`` is one row per structure with the CDR3α/CDR3β/peptide
    sequences + species (the inputs to non-redundancy clustering and the benchmarks).
    Species is auto-detected per structure by alignment score (human vs mouse). All chains
    across the whole folder are annotated in a single mmseqs call per organism (the
    per-call process overhead dominates, so dataset-level batching is far faster than
    per-structure annotation).
    """
    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..structure import parse_structure

    struct_dir = Path(struct_dir)
    paths = sorted(p for p in struct_dir.iterdir() if p.suffix in (".pdb", ".cif"))
    structures: list[Structure] = []
    for path in paths:
        # TCR3D CIFs are named "<pdb>_renumbered.cif"; normalise to the 4-char PDB id.
        pdb_id = path.stem.split("_")[0]
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
            ct = contact_table(s)
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
