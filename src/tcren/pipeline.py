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

import numpy as np
import polars as pl

from .annotation import classify_chains
from .contactmap import ContactMap
from .contacts.table import residue_annotation
from .mhc import MhcCall, annotate_mhc
from .potential import Potential, keskin, mj, tcren
from .structure.io import import_structure
from .structure.model import PEPTIDE_TYPE, Structure

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
    """Everything the pipeline produces for one structure.

    ``extra`` carries the interface-sanity flag when the complex was oriented
    (``superimpose=True``): ``real_interface`` (``bool`` — ``False`` marks assay noise /
    a failed dock; see :func:`tcren.binder.is_real_interface`), and the raw descriptors it
    was computed from (``n_contacts``, ``scanning_angle``, ``pitch_angle``). With
    ``superimpose=False`` the docking angles are unavailable, so ``real_interface`` is
    ``None`` (``scanning_angle``/``pitch_angle`` ``None``) while ``n_contacts`` is still set.
    """

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
        weights = np.asarray(contacts["n_atom_contacts"].to_list(), dtype=np.float64)
    else:
        weights = np.ones(contacts.height, dtype=np.float64)
    # Vectorized gather off the dense matrix instead of a per-row polars filter
    # (Potential.value): O(contacts) lookups, not O(contacts × potential_rows). Pairs whose
    # residue is outside the alphabet, or absent from the matrix (nan), are dropped — exactly
    # as the per-row path skipped KeyError/IndexError.
    matrix, index = potential.as_matrix()
    rows_idx = np.array([index.get(a, -1) for a in contacts["residue.aa.from"].to_list()],
                        dtype=np.int64)
    cols_idx = np.array([index.get(b, -1) for b in contacts["residue.aa.to"].to_list()],
                        dtype=np.int64)
    valid = (rows_idx >= 0) & (cols_idx >= 0)
    vals = matrix[rows_idx[valid], cols_idx[valid]] * weights[valid]
    return float(np.nansum(vals))


def run(
    structure: str | Path | Structure,
    organism: str = "human",
    superimpose: bool = True,
    db_dir: str | Path | None = None,
    cutoff: float = 5.0,
    potentials: dict[str, str | Potential | None] | None = None,
    tcr_regions: str = "all",
    contact_weight: str = "residue",
    reference_aa: str | None = None,
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
        reference_aa: if set (typically ``"A"``), also report the reference-normalised
            energies ``delta_<interface>`` and ``delta_total`` --- each interface's
            :func:`tcren.ddg.reference_delta`, i.e. its energy minus the energy of a
            poly-``reference_aa`` peptide threaded onto the same contact map. Off by
            default, so the default ``scores`` dict is unchanged.

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

    if reference_aa is not None:
        # ΔF = F(peptide) − F(poly-reference peptide) per interface, on THIS structure's own
        # contact map. On a fixed map it is F minus a constant and changes no ranking; it only
        # bites across candidates that each carry their own pose (see ddg.reference_delta).
        # ΔF_tcr_mhc is identically 0 — the peptide is not in that interface.
        from .ddg import reference_delta

        peptide = next((c.sequence() for c in s.chains if c.chain_type == PEPTIDE_TYPE), None)
        if peptide is None:
            raise ValueError(f"{s.pdb_id}: no peptide chain, cannot compute a ΔF reference")
        for iface in _INTERFACE_POTENTIAL:
            scores[f"delta_{iface}"] = reference_delta(
                cm, peptide, resolved[iface], interface=iface, reference_aa=reference_aa,
                tcr_regions=tcr_regions, contact_weight=contact_weight,
            )
        scores["delta_total"] = sum(scores[f"delta_{i}"] for i in _INTERFACE_POTENTIAL)

    # Interface-sanity (assay-noise) flag: a cheap pre-energy check that the TCR:peptide
    # interface is a plausible dock (enough contacts + in-range docking geometry). The docking
    # angles only exist once the complex is oriented, so this is a no-op (real_interface=None)
    # when superimpose=False or the geometry is degenerate — we flag, never drop or raise.
    n_contacts = cm.interface("tcr_peptide", tcr_regions=tcr_regions).height
    scanning_angle = pitch_angle = None
    real_interface = None
    if superimpose:
        from .binder.noise import is_real_interface
        from .orient.docking import docking_angles

        try:
            angles = docking_angles(s)
            scanning_angle, pitch_angle = angles.crossing_angle, angles.incident_angle
        except ValueError:  # missing receptor pair / degenerate frame -> geometry unknown
            pass
        real_interface = is_real_interface(n_contacts, scanning_angle, pitch_angle)

    return PipelineResult(
        pdb_id=s.pdb_id, mhc_calls=calls, markup=residue_annotation(s),
        contacts=cm.contacts, scores=scores, oriented=oriented, rmsd=rmsd,
        extra={
            "real_interface": real_interface,
            "n_contacts": n_contacts,
            "scanning_angle": scanning_angle,
            "pitch_angle": pitch_angle,
        },
    )


def score_row(result: PipelineResult) -> dict:
    """Flatten a :class:`PipelineResult` to a one-row scores dict (for a CSV table).

    The ``d_*`` reference-normalised columns are present only when the pipeline was run with
    ``reference_aa`` set.
    """
    mhc = next((c for c in result.mhc_calls if c.chain_role == "MHCa"), None)
    row = {
        "pdb.id": result.pdb_id,
        "mhc.class": mhc.mhc_class if mhc else None,
        "allele": mhc.allele if mhc else None,
        "rmsd": result.rmsd,
        "tcr_peptide.tcren": result.scores["tcr_peptide"],
        "tcr_mhc.mj": result.scores["tcr_mhc"],
        "peptide_mhc.mj": result.scores["peptide_mhc"],
        "total": result.scores["total"],
    }
    if "delta_total" in result.scores:
        row.update({
            "d_tcr_peptide.tcren": result.scores["delta_tcr_peptide"],
            "d_tcr_mhc.mj": result.scores["delta_tcr_mhc"],
            "d_peptide_mhc.mj": result.scores["delta_peptide_mhc"],
            "d_total": result.scores["delta_total"],
        })
    return row
