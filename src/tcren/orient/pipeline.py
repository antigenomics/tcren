"""Orchestrate canonicalization of TCR-pMHC structures into the common MHC frame."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from .align import apply_transform
from ..structure.io import structure_id_from_path
from ..structure.model import Structure
from .chains import _has_multiple_copies, rename_chains, select_primary_complex
from .exceptions import detect_reverse_dock
from .frame import CanonResult, canonical_frame


def canonicalize_structure(
    structure: Structure,
    reference_id: str | None = None,
    force_pca: bool = False,
    select_primary: bool = True,
) -> tuple[Structure, CanonResult]:
    """Orient an (already chain-typed + MHC-annotated) structure into the canonical frame.

    Returns the oriented, A–E renamed structure and the populated :class:`CanonResult`
    (transform, frame, rmsd, reverse-dock flag, chain map). Coordinates are transformed; the
    chain roles drive the rename, so order matters (frame + reverse-dock are read before the
    transform clears region markup).
    """
    s = select_primary_complex(structure) if select_primary else structure
    result = canonical_frame(s, reference_id, force_pca)
    result.reversed_dock = detect_reverse_dock(s, result.rotation, result.translation)
    oriented = apply_transform(s, result)
    oriented.cell_type = s.cell_type
    oriented, chain_map = rename_chains(oriented)
    result.chain_map = chain_map
    return oriented, result


def align_to_canonical(
    structure: Structure,
    reference_id: str | None = None,
    organism: str = "human",
    force_pca: bool = False,
) -> tuple[Structure, CanonResult]:
    """Align a NEW (parsed) structure onto the Native2026 canonical frame.

    Runs chain typing + MHC annotation, then :func:`canonicalize_structure`. The stored
    per-class ``R_canon`` is reused, so the result is in the same frame as the dataset and the
    composed transform in the returned :class:`CanonResult` replays the placement exactly.
    """
    from ..annotation import classify_chains
    from ..mhc import annotate_mhc

    classify_chains(structure, organism=organism)
    annotate_mhc(structure)
    return canonicalize_structure(structure, reference_id=reference_id, force_pca=force_pca)


def check_oriented_complex(structure, max_peptide_len: int = 25, max_offset: float = 25.0,
                           max_tcr_gap: float = 15.0, max_orphan: float = 70.0):
    """Geometric sanity check on an oriented A–E complex; ``(ok, reason)``.

    Rejects structures whose canonical placement is inconsistent: missing / overlong peptide,
    peptide not at the groove centre (≈ origin), the TCR not engaging the peptide, or any chain
    stranded far from the complex (an orphan copy that survived primary-complex selection).
    """
    import numpy as np

    ca = {c.chain_id: np.asarray([r.ca for r in c.residues if r.ca is not None])
          for c in structure.chains}
    pep = ca.get("C")
    if pep is None or len(pep) < 2:
        return False, "no_peptide"
    if len(pep) > max_peptide_len:
        return False, "peptide_too_long"
    if np.linalg.norm(pep.mean(axis=0)) > max_offset:
        return False, "peptide_off_center"
    tcr = [ca[k] for k in ("A", "B") if k in ca and len(ca[k])]
    if tcr:
        t = np.vstack(tcr)
        if float(np.min(np.linalg.norm(t[:, None, :] - pep[None, :, :], axis=2))) > max_tcr_gap:
            return False, "tcr_not_engaged"
    for cid, arr in ca.items():
        if len(arr) and np.linalg.norm(arr.mean(axis=0)) > max_orphan:
            return False, f"orphan_chain_{cid}"
    return True, "ok"


_ROW_KEYS =("pdb.id", "status", "mhc.class", "species", "tcr.type", "cell.type", "frame",
             "reference.id", "rmsd", "n.anchor.atoms", "reversed.dock", "n.copies",
             "chain.map", "transform")


def _orient_row(pdb_id, status="ok"):
    row = dict.fromkeys(_ROW_KEYS)
    row["pdb.id"], row["status"] = pdb_id, status
    return row


def _finish_orient(structure, pdb_id, out, reference_id, force_pca,
                   mmcif=False, compress=True) -> dict:
    """Canonicalize an (already classified + MHC-annotated) structure → metadata row."""
    from ..structure.io import structure_output_path, write_structure

    row = _orient_row(pdb_id)
    row["mhc.class"] = "MHCII" if any(c.chain_type == "MHCb" for c in structure.chains) else "MHCI"
    row["species"] = structure.complex_species
    loci = {c.chain_type for c in structure.chains if c.chain_type in ("TRA", "TRB", "TRD", "TRG")}
    row["tcr.type"] = ("ab" if {"TRA", "TRB"} <= loci and not (loci & {"TRD", "TRG"})
                       else "gd" if {"TRD", "TRG"} <= loci else "other")
    row["n.copies"] = 2 if _has_multiple_copies(structure) else 1
    oriented, res = canonicalize_structure(structure, reference_id=reference_id,
                                           force_pca=force_pca)
    ok, reason = check_oriented_complex(oriented)
    row.update({"frame": res.frame, "reference.id": res.reference_id, "rmsd": res.rmsd,
                "n.anchor.atoms": res.n_anchor_atoms, "reversed.dock": res.reversed_dock,
                "cell.type": oriented.cell_type, "chain.map": json.dumps(res.chain_map),
                "transform": json.dumps({"rotation": res.rotation.tolist(),
                                         "translation": res.translation.tolist()})})
    if ok:
        write_structure(oriented, structure_output_path(out, pdb_id, mmcif, compress))
    else:
        row["status"] = f"rejected: {reason}"
    return row


def run_folder(
    structures: str | Path,
    out: str | Path,
    metadata: str | Path | None = None,
    organism: str = "human",
    reference_id: str | None = None,
    force_pca: bool = False,
    threads: int | None = None,
    mmcif: bool = False,
    compress: bool = False,
) -> pl.DataFrame:
    """Canonicalize a file or folder of structures; write oriented structures + a metadata table.

    Output format follows ``mmcif`` (``.cif`` vs ``.pdb``) and ``compress`` (trailing ``.gz``);
    plain PDB by default (pass ``compress=True`` to rebuild the gzipped Canonical2026 set).

    Annotation is BATCHED — one mmseqs search for all TCR chains (per organism) and one for all
    MHC chains across the whole set (mmseqs parallelises internally; never per-structure, never
    Python-threaded). Only the embarrassingly-parallel, mmseqs-free stages — parsing and the
    structural alignment + write — use a thread pool (``threads`` worker threads, default
    ``os.cpu_count()``).
    """
    import os
    from concurrent.futures import ThreadPoolExecutor

    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..mhc import annotate_mhc_batch
    from ..paper.helpers import _batch_annotate
    from ..structure.io import import_structure

    from ..structure.io import structure_paths

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    files = structure_paths(Path(structures))
    threads = threads or (os.cpu_count() or 1)

    # 1. Parse (I/O + gunzip — thread-friendly).
    def _parse(fp):
        try:
            return structure_id_from_path(fp), import_structure(fp, pdb_id=structure_id_from_path(fp),
                                                                keep_c_gene=True), None
        except Exception as exc:  # noqa: BLE001
            return structure_id_from_path(fp), None, f"error: {type(exc).__name__}: {str(exc)[:80]}"

    with ThreadPoolExecutor(max_workers=threads) as ex:
        parsed_all = list(ex.map(_parse, files))
    parsed = [(pid, s) for pid, s, err in parsed_all if s is not None]
    rows = [_orient_row(pid, err) for pid, s, err in parsed_all if s is None]

    # 2. Batched annotation: one mmseqs pass for TCR chains, one for MHC chains (no threads).
    records = _batch_annotate([s for _, s in parsed], _import_arda())
    for (_pid, s), recs in zip(parsed, records):
        classify_chains(s, organism=organism, precomputed_records=recs)
    annotate_mhc_batch([s for _, s in parsed])

    # 3. Structural alignment + write (CPU/SVD + I/O — thread-friendly, mmseqs-free).
    def _orient(item):
        pid, s = item
        try:
            return _finish_orient(s, pid, out, reference_id, force_pca, mmcif, compress)
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            return _orient_row(pid, f"error: {type(exc).__name__}: {str(exc)[:80]}")

    with ThreadPoolExecutor(max_workers=threads) as ex:
        rows += list(ex.map(_orient, parsed))

    df = pl.DataFrame(rows)
    if metadata is not None:
        df.write_csv(metadata)
    ok = df.filter(pl.col("status") == "ok").height
    print(f"oriented {ok}/{df.height} structures -> {out}")
    return df


def run_superimpose(
    structures: str | Path,
    out: str | Path,
    db_dir: str | Path | None = None,
    organism: str = "human",
    mmcif: bool = False,
    compress: bool = False,
) -> pl.DataFrame:
    """Superimpose each input structure onto a canonical database; write oriented structures.

    Thin folder/file driver over :func:`tcren.orient.superimpose` (see its module docstring for
    the MHC-ensemble method). ``db_dir`` defaults to ``data/Canonical2026``. Output format
    follows ``mmcif``/``compress`` (plain PDB by default for user-facing output).
    """
    from ..orient.superimpose import superimpose
    from ..structure.io import iter_structures, parse_structure, structure_output_path, write_structure

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for pid, s in iter_structures(structures, importer=parse_structure):
        try:
            oriented, res = superimpose(s, db_dir=db_dir, organism=organism)
            ok, reason = check_oriented_complex(oriented)
            if ok:
                write_structure(oriented, structure_output_path(out, pid, mmcif, compress))
                rows.append({"pdb.id": pid, "status": "ok", "reference.id": res.reference_id,
                             "rmsd": res.rmsd, "n.references": res.n_anchor_atoms})
            else:
                rows.append({"pdb.id": pid, "status": f"rejected: {reason}",
                             "reference.id": res.reference_id, "rmsd": res.rmsd,
                             "n.references": res.n_anchor_atoms})
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            rows.append({"pdb.id": pid, "status": f"error: {type(exc).__name__}: {str(exc)[:80]}",
                         "reference.id": None, "rmsd": None, "n.references": None})
    df = pl.DataFrame(rows)
    ok = df.filter(pl.col("status") == "ok").height if df.height else 0
    print(f"superimposed {ok}/{df.height} structures -> {out}")
    return df
