"""End-to-end TCRen pipeline: structure → annotation → orientation → contacts → score.

One call takes a TCR-pMHC structure all the way through the tcren workflow:

1. **import** the structure (C-gene trimmed);
2. **annotate** chains — TCR loci/CDRs via arda, MHC allele/class/role + groove regions;
3. **superimpose** onto the canonical database (canonical Cα frame; optional);
4. **markup + contacts** — the per-residue region table and the 5 Å contact map;
5. **score** each interface with its potential: TCRen for TCR↔peptide, MJ for TCR↔MHC and
   peptide↔MHC, plus the total.

The interface energy is the sum of the residue-pair potential over the observed contacts of
that interface (the closest-atom contact per residue pair, as everywhere in tcren).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from .annotation import classify_chains
from .contactmap import ContactMap
from .contacts.table import residue_annotation
from .mhc import MhcCall, annotate_mhc
from .potential import Potential, keskin, mj, tcren
from .structure.io import import_structure
from .structure.model import Structure

# Interface → potential family (TCRen for the TCR↔peptide contact map; MJ elsewhere).
_INTERFACE_POTENTIAL = {"tcr_peptide": "tcren", "tcr_mhc": "mj", "peptide_mhc": "mj"}

# Bundled potential loaders, keyed by the name accepted in the ``potentials`` spec.
_BUNDLED_POTENTIALS = {"tcren": tcren, "mj": mj, "keskin": keskin}


def _resolve_potentials(
    spec: dict[str, str | Potential | None] | None,
) -> dict[str, Potential]:
    """Resolve a per-interface potential spec to ``{interface: Potential}``.

    Args:
        spec: Maps an interface name (``"tcr_peptide"``, ``"tcr_mhc"``, ``"peptide_mhc"``)
            to a :class:`Potential`, a bundled name (``"tcren"``/``"mj"``/``"keskin"``),
            a CSV path, or ``None``. A missing or ``None`` entry falls back to the default
            :data:`_INTERFACE_POTENTIAL` family for that interface.

    Returns:
        One resolved :class:`Potential` per interface in :data:`_INTERFACE_POTENTIAL`.
    """
    spec = spec or {}
    cache: dict[str, Potential] = {}

    def _load(value: str | Potential) -> Potential:
        if isinstance(value, Potential):
            return value
        if value in _BUNDLED_POTENTIALS:
            if value not in cache:
                cache[value] = _BUNDLED_POTENTIALS[value]()
            return cache[value]
        return Potential.from_csv(value)

    resolved: dict[str, Potential] = {}
    for iface, default_fam in _INTERFACE_POTENTIAL.items():
        value = spec.get(iface)
        resolved[iface] = _load(default_fam if value is None else value)
    return resolved


@dataclass(slots=True)
class PipelineResult:
    """Everything the pipeline produces for one structure."""

    pdb_id: str
    mhc_calls: list[MhcCall]
    markup: pl.DataFrame
    contacts: pl.DataFrame
    scores: dict[str, float]
    oriented: Structure | None = None
    rmsd: float | None = None
    extra: dict = field(default_factory=dict)


def _interface_energy(
    contacts: pl.DataFrame, potential: Potential, contact_weight: str = "residue"
) -> float:
    """Sum the residue-pair ``potential`` over an interface's contacts (unknown residues skipped).

    With ``contact_weight="residue"`` (default, legacy) each contacting residue pair adds
    ``potential[a, b]``. With ``contact_weight="atomic"`` each pair is multiplied by its
    ``n_atom_contacts`` heavy-atom-pair count (which the contacts table must carry).
    """
    if contacts.is_empty():
        return 0.0
    if contact_weight == "atomic":
        if "n_atom_contacts" not in contacts.columns:
            raise ValueError(
                "contact_weight='atomic' needs the n_atom_contacts column; build the "
                "contact map with count_atoms=True"
            )
        weights = contacts["n_atom_contacts"].to_list()
    else:
        weights = [1] * contacts.height
    total = 0.0
    for a, b, w in zip(
        contacts["residue.aa.from"], contacts["residue.aa.to"], weights
    ):
        try:
            total += potential.value(a, b) * w
        except (KeyError, IndexError):  # X / non-standard residue not in the potential
            continue
    return total


def run(
    structure: str | Path | Structure,
    organism: str = "human",
    superimpose: bool = True,
    db_dir: str | Path | None = None,
    cutoff: float = 5.0,
    potentials: dict[str, str | Potential | None] | None = None,
    tcr_regions: str = "all",
    contact_weight: str = "residue",
) -> PipelineResult:
    """Run the full pipeline on one structure (path or parsed :class:`Structure`).

    Args:
        structure: a structure file (any tcren-readable format) or an already-parsed structure.
        organism: organism for TCR annotation.
        superimpose: also orient onto the canonical database (sets ``oriented`` + ``rmsd``).
        db_dir: canonical database for ``superimpose`` (default ``data/Canonical2026``).
        cutoff: contact distance threshold (Å).
        potentials: optional per-interface potential override mapping an interface name
            (``"tcr_peptide"``, ``"tcr_mhc"``, ``"peptide_mhc"``) to a :class:`Potential`,
            a bundled name (``"tcren"``/``"mj"``/``"keskin"``), a CSV path, or ``None``.
            ``None`` (or a missing entry) keeps the default family for that interface, so
            the default output is unchanged.
        tcr_regions: which TCR regions to keep on the TCR side of the TCR-containing
            interfaces (``"all"`` default = no filter = legacy behaviour; ``"cdr"`` or
            ``"cdr+fr"`` to restrict).
        contact_weight: ``"residue"`` (default, legacy) weights each contacting residue
            pair by 1 on **all three** interfaces; ``"atomic"`` weights each pair by its
            ``n_atom_contacts`` heavy-atom-pair count (the contact map is then built with
            ``count_atoms=True``). Applies to ``tcr_peptide``, ``tcr_mhc`` and
            ``peptide_mhc`` alike.

    Returns:
        A :class:`PipelineResult` with the markup, contacts, per-interface scores and (if
        requested) the canonical-frame oriented structure.
    """
    if contact_weight not in ("residue", "atomic"):
        raise ValueError(f"contact_weight must be 'residue' or 'atomic', got {contact_weight!r}")
    s = structure if isinstance(structure, Structure) else import_structure(structure)
    classify_chains(s, organism=organism)
    calls = annotate_mhc(s)

    oriented = rmsd = None
    if superimpose:
        from .orient import superimpose as _superimpose

        oriented, info = _superimpose(s, db_dir=db_dir, organism=organism)
        rmsd = info.rmsd

    cm = ContactMap.from_structure(
        s, cutoff=cutoff, count_atoms=(contact_weight == "atomic")
    )
    resolved = _resolve_potentials(potentials)
    scores = {
        iface: _interface_energy(
            cm.interface(iface, tcr_regions=tcr_regions),
            resolved[iface],
            contact_weight=contact_weight,
        )
        for iface in _INTERFACE_POTENTIAL
    }
    scores["total"] = sum(scores.values())

    return PipelineResult(
        pdb_id=s.pdb_id, mhc_calls=calls, markup=residue_annotation(s),
        contacts=cm.contacts, scores=scores, oriented=oriented, rmsd=rmsd,
    )


def score_row(result: PipelineResult) -> dict:
    """Flatten a :class:`PipelineResult` to a one-row scores dict (for a CSV table)."""
    mhc = next((c for c in result.mhc_calls if c.chain_role == "MHCa"), None)
    return {
        "pdb.id": result.pdb_id,
        "mhc.class": mhc.mhc_class if mhc else None,
        "allele": mhc.allele if mhc else None,
        "rmsd": result.rmsd,
        "tcr_peptide.tcren": result.scores["tcr_peptide"],
        "tcr_mhc.mj": result.scores["tcr_mhc"],
        "peptide_mhc.mj": result.scores["peptide_mhc"],
        "total": result.scores["total"],
    }
