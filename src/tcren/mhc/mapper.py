"""Map a structure's MHC chains to allele / class / role via mmseqs.

Searches each not-yet-typed chain against the curated MHC reference and assigns the
best hit's class (MHCI/MHCII), chain role (MHCa/MHCb/B2M), locus and allele. Class is
reconciled across the complex (B2M ⇒ class I; a class-II beta chain ⇒ class II).
"""

from __future__ import annotations

import os
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
    # Co-best alleles among returned search hits; not an exhaustive genotype call.
    allele_candidates: tuple[str, ...] = ()
    representative_allele: str | None = None
    allele_source: str = "sequence_representative"
    supplied_support: str | None = None
    query_coverage: float | None = None
    exact_observed_match: bool = False
    search_hit_limit_reached: bool = False


def _candidate_chains(structure: Structure):
    """Chains eligible for MHC mapping: not a receptor and not the peptide."""
    return [
        c
        for c in structure.chains
        if c.chain_type not in RECEPTOR_TYPES and c.chain_type != PEPTIDE_TYPE
    ]


def _best_hits(tsv: Path) -> dict[str, dict]:
    """Keep co-best allele evidence and a deterministic representative alignment.

    Truncated groove sequences commonly match several alleles exactly, including
    null alleles whose distinguishing mutation lies outside the observed chain.
    Search output order must not turn that ambiguity into a null-allele call.
    """
    if not tsv.exists() or tsv.stat().st_size == 0:
        return {}
    df = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=_COLUMNS,
                     schema_overrides={"query": pl.String, "target": pl.String})
    counts = dict(df.group_by("query").len().iter_rows())
    best = {}
    for row in df.sort(["bits", "target", "tstart", "qstart"],
                       descending=[True, False, False, False]).iter_rows(named=True):
        query = row["query"]
        allele = reference.parse_header(row["target"])["allele"]
        if query not in best:
            row["allele_candidates"] = [allele]
            row["search_hit_limit_reached"] = counts[query] >= 50
            best[query] = row
        elif row["bits"] == best[query]["bits"]:
            candidates = best[query]["allele_candidates"]
            if allele not in candidates:
                candidates.append(allele)
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
    from ..provenance import searched_index_identity
    identity = searched_index_identity(ref_db)

    calls: list[MhcCall] = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        query_fa = tmp / "query.fasta"
        with query_fa.open("w") as fh:
            for chain in candidates:
                fh.write(f">{chain.chain_id}\n{chain.sequence()}\n")
        out_tsv = tmp / "hits.tsv"
        mmseqs.easy_search(query_fa, ref_db, out_tsv, tmp / "mmseqs_tmp",
                           search_type=1, sensitivity=sensitivity, max_seqs=50,
                           threads=int(os.environ.get("TCREN_MMSEQS_THREADS", "1")))
        best = _best_hits(out_tsv)

    if searched_index_identity(ref_db) != identity:
        raise RuntimeError("MHC annotation index changed during the search")
    structure.annotation_provenance["mhc"] = identity

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
                allele_candidates=tuple(hit.get("allele_candidates", (meta["allele"],))),
                representative_allele=meta["allele"],
                query_coverage=(int(hit["qend"]) - int(hit["qstart"]) + 1) / int(hit["qlen"])
                if hit.get("qlen") else None,
                exact_observed_match=(float(hit["pident"]) == 100.0
                    and int(hit.get("mismatch", 1)) == 0 and int(hit.get("gapopen", 1)) == 0
                    and int(hit["qstart"]) == 1 and int(hit["qend"]) == int(hit.get("qlen", -1))),
                search_hit_limit_reached=hit.get("search_hit_limit_reached", False),
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
    previous = {c.chain_id: c for c in structure.mhc_calls or []}
    for chain in structure.chains:
        call = by_id.get(chain.chain_id)
        if call is not None:
            old = previous.get(chain.chain_id)
            supplied = chain.allele_info
            if supplied and (old is None or old.allele_source == "supplied"
                             or supplied != old.allele):
                from .metadata import set_supplied_allele
                set_supplied_allele(call, supplied, chain.sequence())
            chain.chain_type = call.chain_role
            chain.chain_supertype = call.mhc_class
            chain.allele_info = call.allele
