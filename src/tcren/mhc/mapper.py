"""Map a structure's MHC chains to allele / class / role via mmseqs.

Searches each not-yet-typed chain against the curated MHC reference and assigns the
best hit's class (MHCI/MHCII), chain role (MHCa/MHCb/B2M), locus and allele. Class is
reconciled across the complex (B2M ⇒ class I; a class-II beta chain ⇒ class II).
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ..structure.model import PEPTIDE_TYPE, RECEPTOR_TYPES, Structure
from . import reference

_COLUMNS = [
    "query", "target", "qstart", "qend", "tstart", "tend", "qlen", "tlen",
    "alnlen", "mismatch", "gapopen", "cigar", "qaln", "taln", "evalue", "bits", "pident",
]


@dataclass(slots=True)
class MhcCall:
    """Result of mapping one chain to the MHC reference."""

    chain_id: str
    chain_role: str  # MHCa | MHCb | B2M
    mhc_class: str  # MHCI | MHCII
    allele: str
    locus: str
    species: str
    identity: float
    bits: float
    # Alignment of query (structure chain) to the reference target, for region projection.
    qstart: int
    qend: int
    tstart: int
    tend: int
    cigar: str


def _candidate_chains(structure: Structure):
    """Chains eligible for MHC mapping: not a receptor and not the peptide."""
    return [
        c
        for c in structure.chains
        if c.chain_type not in RECEPTOR_TYPES and c.chain_type != PEPTIDE_TYPE
    ]


def _best_hits(tsv: Path) -> dict[str, dict]:
    """Parse the mmseqs TSV, keeping the highest-bitscore hit per query chain."""
    if not tsv.exists() or tsv.stat().st_size == 0:
        return {}
    df = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=_COLUMNS)
    best = {}
    for row in df.sort("bits", descending=True).iter_rows(named=True):
        best.setdefault(row["query"], row)
    return best


def map_mhc(structure: Structure, sensitivity: float = 5.7) -> list[MhcCall]:
    """Map the structure's MHC chains against the curated reference.

    Args:
        structure: A structure whose TCR/peptide chains are already typed.
        sensitivity: mmseqs search sensitivity.

    Returns:
        One :class:`MhcCall` per chain that produced a reference hit.
    """
    import arda.mmseqs as mmseqs  # reuse arda's mmseqs wrapper

    candidates = _candidate_chains(structure)
    if not candidates:
        return []
    ref_db = reference.reference_db()  # prebuilt mmseqs DB (createdb cached), not the FASTA

    calls: list[MhcCall] = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        query_fa = tmp / "query.fasta"
        with query_fa.open("w") as fh:
            for chain in candidates:
                fh.write(f">{chain.chain_id}\n{chain.sequence()}\n")
        out_tsv = tmp / "hits.tsv"
        mmseqs.easy_search(query_fa, ref_db, out_tsv, tmp / "mmseqs_tmp",
                           search_type=1, sensitivity=sensitivity, max_seqs=50)
        best = _best_hits(out_tsv)

    calls = calls_from_hits(candidates, best)
    return calls


def calls_from_hits(candidates, best: dict[str, dict], key=None) -> list[MhcCall]:
    """Build reconciled :class:`MhcCall`s for ``candidates`` from precomputed mmseqs hits.

    ``key(chain) -> str`` maps a candidate chain to its key in ``best`` (default the chain id;
    a batched search uses ``"<struct_idx>|<chain_id>"``). Lets one mmseqs search over many
    structures' chains be sliced back per structure — no per-structure mmseqs call.
    """
    key = key or (lambda c: c.chain_id)
    calls: list[MhcCall] = []
    for chain in candidates:
        hit = best.get(key(chain))
        if hit is None:
            continue
        meta = reference.parse_header(hit["target"])
        calls.append(
            MhcCall(
                chain_id=chain.chain_id,
                chain_role=meta["chain_role"],
                mhc_class=meta["mhc_class"],
                allele=meta["allele"],
                locus=meta["locus"],
                species=meta["species"],
                identity=float(hit["pident"]),
                bits=float(hit["bits"]),
                qstart=int(hit["qstart"]),
                qend=int(hit["qend"]),
                tstart=int(hit["tstart"]),
                tend=int(hit["tend"]),
                cigar=hit["cigar"],
            )
        )
    _reconcile_class(calls)
    return calls


def _reconcile_class(calls: list[MhcCall]) -> None:
    """Reconcile MHC class across the complex (B2M ⇒ I; class-II beta ⇒ II)."""
    roles = {c.chain_role for c in calls}
    if "MHCb" in roles:
        resolved = "MHCII"
    elif "B2M" in roles:
        resolved = "MHCI"
    else:
        return
    for call in calls:
        if call.chain_role != "B2M":
            call.mhc_class = resolved


def apply_mhc_calls(structure: Structure, calls: list[MhcCall]) -> None:
    """Write MHC calls onto the structure's chains in place."""
    by_id = {c.chain_id: c for c in calls}
    for chain in structure.chains:
        call = by_id.get(chain.chain_id)
        if call is not None:
            chain.chain_type = call.chain_role
            chain.chain_supertype = call.mhc_class
            chain.allele_info = call.allele
