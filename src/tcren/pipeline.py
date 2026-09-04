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
from .potential import Potential, keskin, mj, tcren, tcren2
from .energetics.scoring import _interface_energy, _phi_scale
from .structure.io import import_structure
from .structure.model import PEPTIDE_TYPE, Structure

# Interface → potential family (TCRen for the TCR↔peptide contact map; MJ elsewhere).
_INTERFACE_POTENTIAL = {"tcr_peptide": "tcren2", "tcr_mhc": "mj", "peptide_mhc": "mj"}

# Bundled potential loaders, keyed by the name accepted in the ``potentials`` spec. This is the
# ONE map -- the CLI imports it rather than keeping a second one, because when the two drifted
# apart `tcren2` resolved in `tcren score` and fell through to a path lookup in `tcren scoring`,
# which then reported "FileNotFoundError: tcren2" once per structure.
#
# `tcren2` is the current matrix; the 2022 one is `karnaukhov2022` after its paper. The bare name
# `tcren` is deliberately absent: it read as "the TCRen potential" while meaning the 2022 matrix.
_BUNDLED_POTENTIALS = {
    "tcren2": tcren2,
    "karnaukhov2022": tcren,
    "mj": mj,
    "keskin": keskin,
}


def _resolve_potentials(
    spec: dict[str, str | Potential | None] | None,
) -> dict[str, Potential]:
    """Resolve a per-interface potential spec to ``{interface: Potential}``.

    Args:
        spec: Maps an interface name (``"tcr_peptide"``, ``"tcr_mhc"``, ``"peptide_mhc"``, or
            ``"peptide_internal"`` for the intra-peptide term) to a :class:`Potential`, a bundled
            name (``"tcren2"``/``"karnaukhov2022"``/``"mj"``/``"keskin"``), a CSV path, or ``None``.
            A missing or ``None``
            entry falls back to the default family for that interface.

    Returns:
        One resolved :class:`Potential` per key of :data:`_INTERFACE_POTENTIAL` plus
        ``"peptide_internal"``.
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
    # The intra-peptide term is resolvable but is not an interface: it stays out of
    # _INTERFACE_POTENTIAL, which drives the three scores and the total. MJ, because TCRen is
    # derived from TCR↔peptide contacts and says nothing about a chain's contacts with itself.
    for iface, default_fam in {**_INTERFACE_POTENTIAL, "peptide_internal": "mj"}.items():
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
    intra_weight: float = 0.0,
    typed: bool = False,
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
            a bundled name (``"tcren2"``/``"karnaukhov2022"``/``"mj"``/``"keskin"``), a CSV
            path, or ``None``.
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
        typed: the structure's chains are already typed, so skip :func:`classify_chains`.
            Chain typing costs one ``mmseqs easy-search`` per structure and a caller scoring a
            whole cohort should have typed it in one batched call
            (:func:`tcren.paper.iter_annotated_set`) rather than paying that per structure.
            Only meaningful when ``structure`` is an already-parsed :class:`Structure`.
        intra_weight: weight of the intra-peptide term. Non-zero adds
            ``scores["peptide_internal"]`` — the peptide's contact energy with **itself**
            (:func:`tcren.intra_peptide_energy`), which every interface sum omits — and folds
            ``intra_weight *`` that energy into ``scores["total"]``. ``0.0`` (default) computes
            nothing and leaves ``scores`` unchanged. Its potential is MJ unless
            ``potentials["peptide_internal"]`` overrides it.

    Returns:
        A :class:`PipelineResult` with the markup, contacts, per-interface scores and (if
        requested) the canonical-frame oriented structure.
    """
    if contact_weight not in ("residue", "atomic"):
        raise ValueError(f"contact_weight must be 'residue' or 'atomic', got {contact_weight!r}")
    s = structure if isinstance(structure, Structure) else import_structure(structure)
    # `typed` means chain-typed AND MHC-annotated already, by the batched path
    # (`tcren.annotation.batch.iter_typed` + `tcren.mhc.annotate_mhc_batch`). Both spawn mmseqs,
    # so repeating either here is one process per structure -- see CLAUDE.md 0-mmseqs.
    if not typed:
        classify_chains(s, organism=organism)
    calls = getattr(s, "mhc_calls", None) if typed else None
    if calls is None:
        calls = annotate_mhc(s)

    oriented = rmsd = None
    if superimpose:
        from .docking import superimpose as _superimpose

        # annotate=False: `s` is chain-typed and MHC-annotated above. Leaving it True made
        # superimpose run classify_chains AND annotate_mhc a SECOND time per structure -- two more
        # mmseqs processes each, and 8.0 s of the 10.1 s a four-structure run spent.
        oriented, info = _superimpose(s, db_dir=db_dir, organism=organism, annotate=False)
        rmsd = info.rmsd

    cm = ContactMap.from_structure(
        s, cutoff=cutoff, count_atoms=(contact_weight == "atomic"),
        peptide_internal=bool(intra_weight),
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
    # Phi_TCRpMHC = c_TP Phi_TP + c_TM Phi_TM + c_PM Phi_PM. The three interfaces are scored with
    # DIFFERENT potentials -- TCRen2 for recognition, MJ or Keskin for presentation -- and those are
    # Boltzmann-inverted from different contact statistics, so their matrices are not on a common
    # scale (sd 0.4880 / 0.3270 / 1.3181). An unweighted sum is therefore a sum of three quantities
    # in three units. Each coefficient is 1 / the sd of that interface energy over the 374
    # Native2026 crystals -- the scale the quantity actually realises on real complexes, which
    # absorbs both the matrix scale and the interface's size. No label and no fit enter it, the
    # same standing as the reliability block moments. See _phi_scale.
    weights = {i: 1.0 / _phi_scale(i, resolved[i]) for i in _INTERFACE_POTENTIAL}
    scores["total"] = sum(weights[i] * scores[i] for i in _INTERFACE_POTENTIAL)

    if intra_weight:
        # The peptide's contacts with itself: reported raw, folded into the total at its weight,
        # so the term and the weight given to it stay separable in the output.
        from .energetics.scoring import intra_peptide_energy

        scores["peptide_internal"] = intra_peptide_energy(
            cm, resolved["peptide_internal"], contact_weight=contact_weight
        )
        scores["total"] += intra_weight * scores["peptide_internal"] / resolved["peptide_internal"].scale()

    if reference_aa is not None:
        # ΔF = F(peptide) − F(poly-reference peptide) per interface, on THIS structure's own
        # contact map. On a fixed map it is F minus a constant and changes no ranking; it only
        # bites across candidates that each carry their own pose (see ddg.reference_delta).
        # ΔF_tcr_mhc is identically 0 — the peptide is not in that interface.
        from .energetics.mutation import reference_delta

        peptide = next((c.sequence() for c in s.chains if c.chain_type == PEPTIDE_TYPE), None)
        if peptide is None:
            raise ValueError(f"{s.pdb_id}: no peptide chain, cannot compute a ΔF reference")
        for iface in _INTERFACE_POTENTIAL:
            scores[f"delta_{iface}"] = reference_delta(
                cm, peptide, resolved[iface], interface=iface, reference_aa=reference_aa,
                tcr_regions=tcr_regions, contact_weight=contact_weight,
            )
        scores["delta_total"] = sum(weights[i] * scores[f"delta_{i}"] for i in _INTERFACE_POTENTIAL)

    # Interface-sanity (assay-noise) flag: a cheap pre-energy check that the TCR:peptide
    # interface is a plausible dock (enough contacts + in-range docking geometry). The docking
    # angles only exist once the complex is oriented, so this is a no-op (real_interface=None)
    # when superimpose=False or the geometry is degenerate — we flag, never drop or raise.
    n_contacts = cm.interface("tcr_peptide", tcr_regions=tcr_regions).height
    scanning_angle = pitch_angle = None
    real_interface = None
    if superimpose:
        from .binder.noise import is_real_interface
        from .docking.angles import docking_angles

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
    ``reference_aa`` set, and ``Phi_pep_int`` only when it was run with a non-zero ``intra_weight``.
    """
    mhc = next((c for c in result.mhc_calls if c.chain_role == "MHCa"), None)
    row = {
        "pdb.id": result.pdb_id,
        "mhc.class": mhc.mhc_class if mhc else None,
        "allele": mhc.allele if mhc else None,
        "rmsd": result.rmsd,
        # Same names as tcren.recognition.RECOGNITION_FEATURES, so the two tables join and the
        # project has ONE vocabulary for these quantities.
        "Phi_tcr_pep": result.scores["tcr_peptide"],
        "Phi_tcr_mhc": result.scores["tcr_mhc"],
        "Phi_pep_mhc": result.scores["peptide_mhc"],
        "Phi_total": result.scores["total"],
    }
    if "peptide_internal" in result.scores:
        row["Phi_pep_int"] = result.scores["peptide_internal"]
    if "delta_total" in result.scores:
        row.update({
            "dPhi_tcr_pep": result.scores["delta_tcr_peptide"],
            "dPhi_tcr_mhc": result.scores["delta_tcr_mhc"],
            "dPhi_pep_mhc": result.scores["delta_peptide_mhc"],
            "dPhi_total": result.scores["delta_total"],
        })
    return row
