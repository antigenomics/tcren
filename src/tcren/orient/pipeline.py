"""Orchestrate canonicalization of TCR-pMHC structures into the common MHC frame."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from ..native.align import apply_transform
from ..native.database import NativeDatabase
from ..structure.model import Structure
from .chains import _has_multiple_copies, rename_chains, select_primary_complex
from .exceptions import detect_reverse_dock
from .frame import CanonResult, canonical_frame


def canonicalize_structure(
    structure: Structure,
    db: NativeDatabase | None = None,
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
    result = canonical_frame(s, db, reference_id, force_pca)
    result.reversed_dock = detect_reverse_dock(s, result.rotation, result.translation)
    oriented = apply_transform(s, result)
    oriented.cell_type = s.cell_type
    oriented, chain_map = rename_chains(oriented)
    result.chain_map = chain_map
    return oriented, result


def align_to_canonical(
    structure: Structure,
    db: NativeDatabase | None = None,
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
    return canonicalize_structure(structure, db=db, reference_id=reference_id, force_pca=force_pca)


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


def _structure_files(path: Path):
    if path.is_file():
        return [path]
    return sorted(p for p in path.iterdir() if p.suffix.lower() in (".pdb", ".cif", ".ent"))


def run_folder(
    structures: str | Path,
    out: str | Path,
    metadata: str | Path | None = None,
    organism: str = "human",
    reference_id: str | None = None,
    force_pca: bool = False,
    db: NativeDatabase | None = None,
) -> pl.DataFrame:
    """Canonicalize a file or folder of structures; write oriented PDBs + a metadata table."""
    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..mhc import annotate_mhc
    from ..paper.helpers import _batch_annotate
    from ..structure.io import import_structure, write_pdb

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    db = db or NativeDatabase()

    # Parse everything up front, then annotate TCR chains across the whole batch in one
    # mmseqs call per organism (the per-call process overhead dominates per-structure calls).
    files = _structure_files(Path(structures))
    parsed: list[tuple[str, object]] = []
    parse_errors: dict[str, str] = {}
    for fp in files:
        pdb_id = fp.stem.split("_")[0]
        try:
            parsed.append((pdb_id, import_structure(fp, pdb_id=pdb_id, keep_c_gene=True)))
        except Exception as exc:  # noqa: BLE001
            parse_errors[pdb_id] = f"error: {type(exc).__name__}: {str(exc)[:80]}"
    records = _batch_annotate([s for _, s in parsed], _import_arda())

    rows = []
    for (pdb_id, s), recs in zip(parsed, records):
        row = {"pdb.id": pdb_id, "status": "ok", "mhc.class": None, "species": None,
               "tcr.type": None, "cell.type": None, "frame": None, "reference.id": None,
               "rmsd": None, "n.anchor.atoms": None, "reversed.dock": None, "n.copies": None,
               "chain.map": None, "transform": None}
        try:
            classify_chains(s, organism=organism, precomputed_records=recs)
            annotate_mhc(s)
            row["mhc.class"] = "MHCII" if any(c.chain_type == "MHCb" for c in s.chains) else "MHCI"
            row["species"] = s.complex_species
            loci = {c.chain_type for c in s.chains if c.chain_type in ("TRA", "TRB", "TRD", "TRG")}
            row["tcr.type"] = ("ab" if {"TRA", "TRB"} <= loci and not (loci & {"TRD", "TRG"})
                               else "gd" if {"TRD", "TRG"} <= loci else "other")
            row["n.copies"] = 2 if _has_multiple_copies(s) else 1
            oriented, res = canonicalize_structure(s, db=db, reference_id=reference_id,
                                                   force_pca=force_pca)
            ok, reason = check_oriented_complex(oriented)
            row.update({
                "frame": res.frame, "reference.id": res.reference_id,
                "rmsd": res.rmsd, "n.anchor.atoms": res.n_anchor_atoms,
                "reversed.dock": res.reversed_dock, "cell.type": oriented.cell_type,
                "chain.map": json.dumps(res.chain_map),
                "transform": json.dumps({"rotation": res.rotation.tolist(),
                                         "translation": res.translation.tolist()}),
            })
            if ok:
                write_pdb(oriented, out / f"{pdb_id}.pdb")
            else:
                row["status"] = f"rejected: {reason}"
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            row["status"] = f"error: {type(exc).__name__}: {str(exc)[:80]}"
        rows.append(row)
    for pdb_id, status in parse_errors.items():
        rows.append({"pdb.id": pdb_id, "status": status})
    df = pl.DataFrame(rows)
    if metadata is not None:
        df.write_csv(metadata)
    ok = df.filter(pl.col("status") == "ok").height
    print(f"oriented {ok}/{df.height} structures -> {out}")
    return df
