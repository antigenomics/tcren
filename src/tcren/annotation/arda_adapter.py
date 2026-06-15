"""TCR chain annotation via the ``arda`` library.

Extracts each chain's amino-acid sequence, runs arda's AIRR annotation, and projects the
returned region coordinates (1-based, end-inclusive into the input sequence) back onto
structure residues as :class:`~tcren.structure.model.RegionMarkup`. Region names are
mapped to the legacy mir vocabulary (``FR1``/``CDR1``/…/``FR4``).
"""

from __future__ import annotations

# arda region key -> (mir region type). Ordered N→C.
_REGION_MAP: tuple[tuple[str, str], ...] = (
    ("fwr1", "FR1"),
    ("cdr1", "CDR1"),
    ("fwr2", "FR2"),
    ("cdr2", "CDR2"),
    ("fwr3", "FR3"),
    ("cdr3", "CDR3"),
    ("fwr4", "FR4"),
)

# Antigen-receptor loci arda may call. TCR-mimic antibodies (IG*) appear in the dataset
# as the receptor and were treated as TCR by the legacy mir, so we annotate them too.
_TCR_LOCI = ("TRA", "TRB", "TRD", "TRG")
_BCR_LOCI = ("IGH", "IGK", "IGL")
_RECEPTOR_LOCI = _TCR_LOCI + _BCR_LOCI


def _import_arda():
    try:
        import arda  # noqa: PLC0415
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ImportError(
            "the 'arda' package is required for TCR annotation; install it into the "
            "tcren environment (pip install -e ../arda)"
        ) from exc
    return arda


def _apply_record(chain, record: dict) -> dict:
    """Project one arda record onto ``chain`` in place; return the record.

    Sets ``chain_type``/``chain_supertype``/``allele_info`` and the region markup when the
    record is a receptor locus; a non-receptor record is returned untouched.
    """
    from ..structure.model import RegionMarkup

    locus = record.get("locus")
    if locus not in _RECEPTOR_LOCI:
        return record

    chain.chain_type = locus
    chain.chain_supertype = "TRAB" if locus in _TCR_LOCI else "IG"
    v_call, j_call = record.get("v_call"), record.get("j_call")
    if v_call or j_call:
        chain.allele_info = f"{v_call or ''}:{j_call or ''}"

    regions = []
    for arda_name, mir_name in _REGION_MAP:
        start = record.get(f"{arda_name}_start")
        end = record.get(f"{arda_name}_end")
        if start is None or end is None:
            continue
        try:
            start, end = int(start), int(end)
        except (TypeError, ValueError):
            continue  # arda left this region's coordinates unset for this record
        residues = chain.residues[start - 1 : end]  # 1-based inclusive -> slice
        if not residues:
            continue
        regions.append(
            RegionMarkup(
                region_type=mir_name,
                start_seq_index=residues[0].seq_index,
                end_seq_index=residues[-1].seq_index,
                sequence="".join(r.aa for r in residues),
                residues=residues,
            )
        )
    chain.regions = regions
    return record


def annotate_chain(chain, organism: str) -> dict | None:
    """Annotate one chain with arda; return the AIRR record if it is a TCR chain.

    Mutates ``chain`` in place when arda recognises it as TRA/TRB: sets ``chain_type``,
    ``chain_supertype`` (``"TRAB"``), ``allele_info`` and ``regions``. Returns the arda
    record (for any locus) or ``None`` if arda produced no locus.
    """
    seq = chain.sequence()
    if not seq:
        return None
    record = _import_arda().annotate_sequences(
        [(chain.chain_id, seq)], seqtype="aa", organism=organism
    )[0]
    return _apply_record(chain, record)


def apply_records(chains, by_id: dict[str, dict]) -> None:
    """Project a cached ``{chain_id: record}`` map onto chains in place (no arda call)."""
    for chain in chains:
        rec = by_id.get(chain.chain_id)
        if rec is not None:
            _apply_record(chain, rec)


def score_records(chains, by_id: dict[str, dict]) -> tuple[list[str], float]:
    """``(receptor_ids, summed mmseqs2_score)`` from already-computed records."""
    receptor_ids: list[str] = []
    total_score = 0.0
    for chain in chains:
        rec = by_id.get(chain.chain_id)
        if rec is not None and rec.get("locus") in _RECEPTOR_LOCI:
            receptor_ids.append(chain.chain_id)
            score = rec.get("mmseqs2_score")
            if isinstance(score, (int, float)):
                total_score += float(score)
    return receptor_ids, total_score


def annotate_chains(chains, organism: str) -> dict[str, dict]:
    """Annotate a batch of chains in a single arda call; apply records in place.

    One mmseqs invocation for all chains (the per-call process/DB overhead dominates the
    actual alignment, so batching is ~hundreds× faster than per-chain calls). Returns a
    ``{chain_id: record}`` map for chains that had a sequence.
    """
    items = [(c.chain_id, c.sequence()) for c in chains if c.sequence()]
    if not items:
        return {}
    records = _import_arda().annotate_sequences(items, seqtype="aa", organism=organism)
    by_id = {cid: rec for (cid, _), rec in zip(items, records)}
    apply_records(chains, by_id)
    return by_id


def annotate_tcr_chains(structure, organism: str = "human") -> list[str]:
    """Annotate all chains; return ids recognised as antigen-receptor (TCR/BCR) chains."""
    return annotate_tcr_chains_scored(structure, organism)[0]


def annotate_tcr_chains_scored(structure, organism: str = "human") -> tuple[list[str], float]:
    """Annotate all chains; return ``(receptor_ids, summed mmseqs2_score)``.

    The summed mmseqs2 alignment score over the receptor chains measures how well the
    structure's TCR/BCR chains match this organism's germline reference — the signal used
    to pick the correct species when annotating against human vs mouse.
    """
    return score_records(structure.chains, annotate_chains(structure.chains, organism))
