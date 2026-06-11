"""Compute tcren's own chain/complex annotation for native structures.

Used to validate the tcren pipeline against the TCR3D reference tables: for each native
CIF, tcren parses, types the chains (arda for TCR, the MHC mapper for MHC) and reports
the V/J genes, CDR3, MHC class/allele and epitope in a normalised form comparable to
``tcr_chain_data.tsv`` / ``tcr_complexes_data.tsv``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..annotation import classify_chains
from ..mhc import map_mhc
from ..structure import parse_structure
from .database import NativeDatabase


@dataclass(slots=True)
class ChainAnnotation:
    tcr_type: str  # "Alpha" | "Beta"
    v_gene: str  # allele stripped, e.g. "TRAV12-2"
    j_gene: str
    cdr3: str  # arda CDR3 (without the conserved C…F/W anchors)


@dataclass(slots=True)
class ComplexAnnotation:
    pdb_id: str
    mhc_class: str | None  # "CLASSI" | "CLASSII"
    mhc_allele: str | None  # MHCa allele, e.g. "HLA-A*02:01"
    epitope: str | None
    chains: list[ChainAnnotation] = field(default_factory=list)

    def chain(self, tcr_type: str) -> ChainAnnotation | None:
        return next((c for c in self.chains if c.tcr_type == tcr_type), None)


def _strip_allele(gene: str | None) -> str:
    return gene.split("*")[0] if gene else ""


def annotate_complex(
    db: NativeDatabase, pdb_id: str, organism: str = "human"
) -> ComplexAnnotation:
    """Annotate a single native complex with the tcren pipeline."""
    import arda

    s = parse_structure(db.cif_for(pdb_id), pdb_id=pdb_id)
    classify_chains(s, organism=organism)

    chains: list[ChainAnnotation] = []
    for chain in s.chains:
        if chain.chain_type in ("TRA", "TRB"):
            rec = arda.annotate_sequences(
                [(chain.chain_id, chain.sequence())], seqtype="aa", organism=organism
            )[0]
            chains.append(
                ChainAnnotation(
                    tcr_type="Alpha" if chain.chain_type == "TRA" else "Beta",
                    v_gene=_strip_allele(rec.get("v_call")),
                    j_gene=_strip_allele(rec.get("j_call")),
                    cdr3=rec.get("cdr3_aa") or "",
                )
            )

    calls = map_mhc(s)
    mhc_class = None
    if calls:
        mhc_class = "CLASSII" if any(c.chain_role == "MHCb" for c in calls) else "CLASSI"
    mhca = next((c for c in calls if c.chain_role == "MHCa"), None)
    epitope = next(
        (c.sequence() for c in s.chains if c.chain_type == "PEPTIDE"), None
    )
    return ComplexAnnotation(
        pdb_id=pdb_id,
        mhc_class=mhc_class,
        mhc_allele=mhca.allele if mhca else None,
        epitope=epitope,
        chains=chains,
    )


def cdr3_core(tcr3d_cdr3: str | None) -> str | None:
    """TCR3D CDR3 stripped of its conserved C (N-term) and F/W (C-term) anchors."""
    if not tcr3d_cdr3:
        return None
    return tcr3d_cdr3[1:-1] if len(tcr3d_cdr3) >= 2 else tcr3d_cdr3


def mhc_locus(allele: str | None) -> str:
    """First-field MHC locus token (``HLA-A*02:01`` → ``HLA-A*02``; ``I-Ak`` → ``I-Ak``)."""
    if not allele:
        return ""
    return allele.split(":")[0]
