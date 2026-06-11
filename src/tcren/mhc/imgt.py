"""Download and parse MHC allele references (IMGT/HLA + UniProt mouse H-2 + B2M).

Produces :class:`MhcAllele` records labelled with species, MHC class and chain role
(``MHCa`` = class-I heavy or class-II alpha; ``MHCb`` = class-II beta; ``B2M``). Human
alleles come from IMGT/HLA (``hla_prot.fasta``); mouse H-2 and beta-2-microglobulin come
from reviewed UniProt entries.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import requests

HLA_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta"
HUMAN_B2M_URL = "https://rest.uniprot.org/uniprotkb/P61769.fasta"
MOUSE_MHC_URL = (
    "https://rest.uniprot.org/uniprotkb/search?query="
    "%28organism_id%3A10090%29+AND+reviewed%3Atrue+AND+%28"
    "gene%3AH2-K1+OR+gene%3AH2-D1+OR+gene%3AH2-L+OR+gene%3AH2-Q1+OR+"
    "gene%3AH2-Aa+OR+gene%3AH2-Ab1+OR+gene%3AH2-Ea+OR+gene%3AH2-Eb1+OR+"
    "gene%3AB2m%29&format=fasta&size=200"
)

# Classical, peptide-presenting HLA loci -> (mhc_class, chain_role).
_HLA_LOCUS_ROLE: dict[str, tuple[str, str]] = {
    **{loc: ("MHCI", "MHCa") for loc in ("A", "B", "C", "E", "F", "G")},
    **{loc: ("MHCII", "MHCa") for loc in ("DRA", "DQA1", "DQA2", "DPA1")},
    **{
        loc: ("MHCII", "MHCb")
        for loc in ("DRB1", "DRB3", "DRB4", "DRB5", "DQB1", "DQB2", "DPB1")
    },
}


@dataclass(frozen=True, slots=True)
class MhcAllele:
    """A reference MHC allele sequence with its functional labels."""

    allele: str  # e.g. "HLA-A*02:01", "H2-Kb", "B2M"
    locus: str  # e.g. "A", "DRB1", "H2-K1"
    mhc_class: str  # "MHCI" | "MHCII"
    chain_role: str  # "MHCa" | "MHCb" | "B2M"
    species: str  # "human" | "mouse"
    sequence: str


def _read_fasta(path: Path):
    """Yield ``(header, sequence)`` pairs from a FASTA file."""
    header, seq = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line[1:], []
        else:
            seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def _download(url: str, dest: Path, force: bool) -> Path:
    if dest.exists() and not force:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    dest.write_text(resp.text)
    return dest


def download_human(cache_dir: Path, force: bool = False) -> Path:
    """Download IMGT/HLA ``hla_prot.fasta`` into ``cache_dir``."""
    return _download(HLA_URL, cache_dir / "hla_prot.fasta", force)


def download_mouse(cache_dir: Path, force: bool = False) -> tuple[Path, Path]:
    """Download mouse H-2 / B2m and human B2M into ``cache_dir``."""
    mouse = _download(MOUSE_MHC_URL, cache_dir / "mouse_mhc.fasta", force)
    human_b2m = _download(HUMAN_B2M_URL, cache_dir / "human_b2m.fasta", force)
    return mouse, human_b2m


def _two_field(allele: str) -> str:
    """Collapse an HLA allele designation to two fields (``A*02:01``)."""
    fields = allele.split(":")
    return ":".join(fields[:2]) if len(fields) >= 2 else allele


def parse_human(path: Path) -> list[MhcAllele]:
    """Parse IMGT/HLA, keeping classical loci collapsed to two-field resolution."""
    seen: set[tuple[str, str]] = set()  # (allele, sequence) dedup
    out: list[MhcAllele] = []
    for header, seq in _read_fasta(path):
        # header: "HLA:HLA00001 A*01:01:01:01 365 bp"
        parts = header.split()
        if len(parts) < 2 or not seq:
            continue
        designation = parts[1]  # "A*01:01:01:01"
        locus = re.split(r"[*]", designation)[0]
        role = _HLA_LOCUS_ROLE.get(locus)
        if role is None:
            continue  # non-classical / non-presenting locus (MICA, TAP, DM, DO, ...)
        mhc_class, chain_role = role
        allele = f"HLA-{_two_field(designation)}"
        key = (allele, seq)
        if key in seen:
            continue
        seen.add(key)
        out.append(
            MhcAllele(allele, locus, mhc_class, chain_role, "human", seq)
        )
    return out


def parse_mouse(mouse_path: Path, human_b2m_path: Path) -> list[MhcAllele]:
    """Parse reviewed mouse H-2 / B2m and human B2M from UniProt FASTA headers."""
    out: list[MhcAllele] = []
    gene_re = re.compile(r"\bGN=(\S+)")
    for header, seq in _read_fasta(mouse_path):
        if not seq:
            continue
        gene_m = gene_re.search(header)
        gene = gene_m.group(1) if gene_m else "?"
        name = header.lower()
        if gene == "B2m" or "beta-2-microglobulin" in name:
            mhc_class, role = "MHCI", "B2M"
        elif "class ii" in name:
            mhc_class = "MHCII"
            role = "MHCb" if "beta chain" in name else "MHCa"
        elif "class i" in name:
            mhc_class, role = "MHCI", "MHCa"
        else:
            continue
        accession = header.split("|")[1] if "|" in header else header.split()[0]
        out.append(MhcAllele(f"{gene}:{accession}", gene, mhc_class, role, "mouse", seq))

    for header, seq in _read_fasta(human_b2m_path):
        if seq:
            out.append(MhcAllele("B2M", "B2M", "MHCI", "B2M", "human", seq))
    return out
