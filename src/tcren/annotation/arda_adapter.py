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


def annotate_chain(chain, organism: str) -> dict | None:
    """Annotate one chain with arda; return the AIRR record if it is a TCR chain.

    Mutates ``chain`` in place when arda recognises it as TRA/TRB: sets ``chain_type``,
    ``chain_supertype`` (``"TRAB"``), ``allele_info`` and ``regions``. Returns the arda
    record (for any locus) or ``None`` if arda produced no locus.
    """
    from ..structure.model import RegionMarkup

    arda = _import_arda()
    seq = chain.sequence()
    if not seq:
        return None
    record = arda.annotate_sequences([(chain.chain_id, seq)], seqtype="aa", organism=organism)[0]
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


def annotate_tcr_chains(structure, organism: str = "human") -> list[str]:
    """Annotate all chains; return ids recognised as antigen-receptor (TCR/BCR) chains."""
    receptor_ids = []
    for chain in structure.chains:
        record = annotate_chain(chain, organism)
        if record is not None and record.get("locus") in _RECEPTOR_LOCI:
            receptor_ids.append(chain.chain_id)
    return receptor_ids
