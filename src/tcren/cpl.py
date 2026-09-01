"""Predict a combinatorial-peptide-library response matrix from one template TCR:pMHC structure.

A positional-scanning combinatorial peptide library (CPL) measures a T-cell clone's peptide
preference one position at a time. For a peptide of length ``L``, each of the ``L x 20``
sublibraries fixes position ``i`` to amino acid ``a`` and leaves **every other position an equimolar
1/20 mixture**, so the measured cell is an ensemble mean,

    R[i, a] = E[ response | x_i = a ] .

This module predicts that matrix from a single deposited or modelled complex. Position ``i`` has a
fixed set of contact partners in the template, so threading each of the twenty residues through the
same contact map and re-reading the potential costs one batched call per interface. Nothing is
re-docked, nothing moves, and nothing is fitted to any assay.

WHAT A CELL IS SCORED WITH. The assay reads **activation**, which needs the peptide presented *and*
the receptor engaged: a substitution that abolishes MHC binding abolishes the response whatever the
receptor would have done. Every cell therefore carries the sum of both peptide-bearing interfaces,

    Phi = Phi(TCR:peptide) + Phi(peptide:MHC),

TCRen over the first and Miyazawa--Jernigan over the second by default. The two channels are
statistically uncorrelated over these cells, so they add rather than duplicate. A position the
receptor never touches is an *anchor*: its TCR term is identically zero, so scoring by the sum
degrades gracefully to presentation alone, and the older "TCRen at receptor-facing positions, MJ at
anchors" partition is this rule's special case rather than a separate mode.

TWO REFERENCE STATES, BOTH USEFUL. A raw ``Phi`` carries a large per-position offset that says only
how many contacts the position makes, so a cell is meaningful only *relative to the other residues
that could sit there*. Two references are offered and both are exact differences of whole-peptide
energies, so the "rest of the peptide" cancels:

``"wild_type"``
    ``Phi(x_{i->wt}) - Phi(x_{i->a})`` -- a **mutation scan** off the residue the template carries.
    This is the neoantigen / epitope-design question: is this substitution better than what is
    there? It is :func:`tcren.ddg.ddg` resolved per position.

``"equimolar"`` (default)
    ``mean_b Phi(x_{i->b}) - Phi(x_{i->a})`` -- referenced to the **1/20 mixture**, which is the
    assay's own null: the CPL background at position ``i`` is exactly that mixture. Use this to
    compare against measured CPL cells. It is also the only one of the two under which the
    template's own residue is an ordinary measurement rather than a forced zero.

The two differ by a per-position constant, and that constant is not noise: it is how far the
template's residue sits above its column's mean. Referencing to the wild type folds that
between-position quantity into every cell of the column, which is the right thing for a mutation
scan and the wrong thing for a comparison against an assay whose background is the mixture.

SIGN. Lower energy is a better binder throughout ``tcren``, and both references are written as
``reference - candidate``, so **positive means favourable**: a positive ``wild_type`` value says the
substitution improves on the template residue, and a positive ``equimolar`` value says the residue is
better than the average residue at that position.

Example:
    >>> from tcren import ContactMap, parse_structure, response_matrix, position_scan
    >>> from tcren.annotation import classify_chains
    >>> from tcren.mhc import annotate_mhc
    >>> s = parse_structure("3HG1.pdb", pdb_id="3HG1")
    >>> classify_chains(s, organism="human"); annotate_mhc(s)
    >>> rm = response_matrix(ContactMap.from_structure(s, cutoff=5.0))
    >>> rm.peptide
    'ELAGIGILTV'
    >>> position_scan(rm, 5).head(3)          # every residue at position 5
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import polars as pl

from .contactmap import ContactMap
from .potential import Potential, mj
from .potential import tcren as _tcren_potential  # noqa: F401  (name clashes with the package)
from .energetics.scoring import _PEPTIDE_SIDE, score_peptides

#: Column order of every matrix this module returns.
AA20: tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWY")

#: The two peptide-bearing interfaces, and the potential family each is scored with by default.
#: ``tcr_mhc`` is absent on purpose: a peptide substitution cannot change it, so its contribution to
#: every cell of the matrix is identically zero.
INTERFACES: tuple[str, ...] = ("tcr_peptide", "peptide_mhc")

#: Accepted values of the ``reference`` argument.
REFERENCES: tuple[str, ...] = ("equimolar", "wild_type")


@dataclass(frozen=True)
class ResponseMatrix:
    """A predicted CPL response matrix: peptide positions x twenty amino acids.

    ``phi[i, a]`` is the **total** interface energy of the template with amino acid ``aa[a]``
    threaded at ``positions[i]``, everything else held at the template sequence and geometry. It is
    a whole-complex energy, not a per-position contribution, so differences *within a row* are exact
    and differences *between rows* are meaningless -- which is why every accessor references a row
    against something in that same row.

    Use :meth:`referenced` (or :func:`mutation_effect` / :func:`position_scan` /
    :func:`equimolar_effect`) rather than reading ``phi`` directly.
    """

    peptide: str                        #: the template's own peptide sequence
    positions: tuple[int, ...]          #: 1-based peptide positions that contact either interface
    interface_class: tuple[str, ...]    #: ``"receptor"`` or ``"anchor"``, one per row
    aa: tuple[str, ...]                 #: the twenty column labels, :data:`AA20`
    phi: object                         #: ``(n_positions, 20)`` total energy, TCR:peptide + peptide:MHC
    phi_tcr: object                     #: the TCR:peptide component alone
    phi_mhc: object                     #: the peptide:MHC component alone
    structure_id: str                   #: the template's ``pdb_id``

    # ---------------------------------------------------------------- lookups
    def row_of(self, position: int) -> int:
        """Row index of a 1-based peptide ``position``, or raise if it contacts nothing.

        Raising is deliberate. A position that touches neither interface has an entirely flat row,
        and returning that silently would report "this position tolerates everything" for what is
        really "this template says nothing about this position".
        """
        try:
            return self.positions.index(position)
        except ValueError:
            raise KeyError(
                f"peptide position {position} makes no contact in {self.structure_id}; "
                f"contacting positions are {list(self.positions)}"
            ) from None

    def column_of(self, aa: str) -> int:
        """Column index of a one-letter amino acid code."""
        try:
            return self.aa.index(aa.upper())
        except ValueError:
            raise KeyError(f"{aa!r} is not one of the twenty amino acids") from None

    def wild_type_at(self, position: int) -> str:
        """The residue the template carries at a 1-based peptide ``position``."""
        return self.peptide[position - 1]

    # ---------------------------------------------------------------- referencing
    def referenced(self, reference: str = "equimolar"):
        """The matrix as effects, ``(n_positions, 20)``, positive = favourable.

        Args:
            reference: ``"equimolar"`` references each row to the mean over its twenty residues --
                the 1/20 mixture the assay holds the other positions at. ``"wild_type"`` references
                it to the residue the template carries, making the matrix a mutation scan whose
                wild-type column is identically zero.

        Returns:
            ``reference_energy - phi``, so a positive entry is a favourable residue.
        """
        if reference not in REFERENCES:
            raise ValueError(f"reference must be one of {REFERENCES}, got {reference!r}")
        phi = np.asarray(self.phi, dtype=float)
        if reference == "equimolar":
            base = np.nanmean(phi, axis=1, keepdims=True)
        else:
            wt = np.array([self.column_of(self.wild_type_at(p)) for p in self.positions])
            base = phi[np.arange(phi.shape[0]), wt][:, None]
        return base - phi

    def to_frame(self, reference: str | None = None) -> pl.DataFrame:
        """Long form: one row per (position, amino acid) cell.

        Args:
            reference: emit only this reference's ``effect`` column. ``None`` (default) emits both
                as ``effect_equimolar`` and ``effect_wild_type``, which is what a side-by-side
                comparison against a measured matrix wants.

        Returns:
            Columns ``structure_id``, ``pos``, ``wt_aa``, ``aa``, ``is_wt``, ``interface_class``,
            ``phi``, ``phi_tcr``, ``phi_mhc``, and the requested effect column(s).
        """
        refs = REFERENCES if reference is None else (reference,)
        eff = {f"effect_{r}": np.asarray(self.referenced(r)).ravel() for r in refs}
        n = len(self.positions)
        pos = np.repeat(np.asarray(self.positions), 20)
        aa = np.tile(np.asarray(self.aa), n)
        wt = np.array([self.wild_type_at(p) for p in self.positions]).repeat(20)
        return pl.DataFrame({
            "structure_id": [self.structure_id] * (n * 20),
            "pos": pos, "wt_aa": wt, "aa": aa, "is_wt": aa == wt,
            "interface_class": np.repeat(np.asarray(self.interface_class), 20),
            "phi": np.asarray(self.phi, dtype=float).ravel(),
            "phi_tcr": np.asarray(self.phi_tcr, dtype=float).ravel(),
            "phi_mhc": np.asarray(self.phi_mhc, dtype=float).ravel(),
            **eff,
        })


# ---------------------------------------------------------------- construction
def _engaged(contact_map: ContactMap, interface: str, tcr_regions: str) -> set[int]:
    """1-based peptide positions that actually touch ``interface``, read off the contact map.

    Read from the contact map and never inferred from an alanine scan: a position whose template
    residue already *is* alanine scores a zero substitution because the mutation is a no-op, not
    because the position makes no contact, and inferring absence from that silently drops it.
    """
    side = _PEPTIDE_SIDE[interface]
    kw = {"tcr_regions": tcr_regions} if interface == "tcr_peptide" else {}
    return {int(p) + 1 for p in contact_map.interface(interface, **kw)[f"pos.{side}"].to_list()}


def _threaded_energies(
    contact_map: ContactMap, peptide: str, potential: Potential, interface: str,
    positions: Sequence[int], tcr_regions: str, contact_weight: str,
    weights: "np.ndarray | None" = None,
):
    """``(len(positions), 20)`` of ``Phi(x_{i->a})`` for one interface, in ONE batched call."""
    variants, keys = [], []
    for i in positions:
        for a in AA20:
            variants.append(peptide[: i - 1] + a + peptide[i:])
            keys.append((i, a))
    kw = {"tcr_regions": tcr_regions} if interface == "tcr_peptide" else {}
    scored = score_peptides(contact_map, variants, potential, interface=interface,
                            contact_weight=contact_weight, weights=weights, **kw)
    by_seq = dict(zip(scored["peptide"].to_list(), scored["score"].to_list()))
    out = np.full((len(positions), 20), np.nan)
    row = {p: k for k, p in enumerate(positions)}
    for (i, a), seq in zip(keys, variants):
        if seq in by_seq:
            out[row[i], AA20.index(a)] = float(by_seq[seq])
    return out


def response_matrix(
    contact_map: ContactMap,
    peptide: str | None = None,
    *,
    tcr_potential: Potential | None = None,
    mhc_potential: Potential | None = None,
    tcr_regions: str = "all",
    contact_weight: str = "residue",
    tcr_weights: "np.ndarray | None" = None,
) -> ResponseMatrix:
    """Predict the CPL response matrix of a template TCR:pMHC complex.

    Every peptide position that contacts either peptide-bearing interface gets a row; each of the
    twenty residues is threaded there in turn on the template's own contact map and scored with both
    potentials, and the two interface energies are summed. See the module docstring for why both
    interfaces enter and what the two reference states mean.

    Args:
        contact_map: the template's contact map, e.g.
            ``ContactMap.from_structure(s, cutoff=5.0)``. The structure must have been chain-typed
            (``tcren.annotation.classify_chains``) and MHC-annotated
            (``tcren.mhc.annotate_mhc``) first, or the peptide:MHC interface comes out empty.
        peptide: sequence to thread substitutions off. Defaults to the template's own peptide, taken
            from the contact map. Must match the template's peptide length.
        tcr_potential: potential for TCR:peptide. Default: the bundled TCRen.
        mhc_potential: potential for peptide:MHC. Default: Miyazawa--Jernigan, which is what the
            presentation interface is scored with throughout ``tcren`` -- TCRen is a *directed*
            TCR-to-peptide recognition potential and does not describe the groove.
        tcr_regions: which TCR regions contribute on the TCR side (``"all"``/``"cdr"``/``"cdr+fr"``).
        contact_weight: ``"residue"`` (default) or ``"atomic"``; passed through to
            :func:`tcren.scoring.score_peptides`.
        tcr_weights: an explicit per-contact multiplier for the **TCR:peptide** interface only, one
            value per row of ``contact_map.interface("tcr_peptide")`` and in its row order. Its use
            is the same as :func:`tcren.ddg.ddg`'s ``weights``: replace the map's hard 0/1 contact
            indicator with a contact **probability** -- :func:`tcren.potts.contact_probabilities`'
            ``p_model`` -- so every threaded substitution is scored against how often each pair
            actually touches rather than against one frozen snapshot of whether it did. The
            presentation interface is left alone, because the shipped Potts model is fitted on
            TCR:peptide. ``None`` (default) leaves the result byte-identical.

    Returns:
        A :class:`ResponseMatrix`.

    Raises:
        ValueError: if the complex has no peptide contacts at all, or if ``peptide`` has the wrong
            length. Both are silent-empty failures otherwise: an empty matrix is indistinguishable
            from a genuinely flat one.
    """
    peptide = peptide or _template_peptide(contact_map)
    if contact_map.peptide_length and len(peptide) != contact_map.peptide_length:
        raise ValueError(
            f"peptide {peptide!r} has length {len(peptide)}, but {contact_map.pdb_id} carries a "
            f"{contact_map.peptide_length}-mer"
        )
    engaged = {name: _engaged(contact_map, name, tcr_regions) for name in INTERFACES}
    positions = tuple(sorted(engaged["tcr_peptide"] | engaged["peptide_mhc"]))
    if not positions:
        raise ValueError(
            f"{contact_map.pdb_id} has no peptide contacts on either interface -- nothing to "
            "predict. Was the structure chain-typed and MHC-annotated before the contact map?"
        )

    phi_tcr = _threaded_energies(contact_map, peptide, tcr_potential or _tcren_potential(),
                                 "tcr_peptide", positions, tcr_regions, contact_weight,
                                 tcr_weights)
    phi_mhc = _threaded_energies(contact_map, peptide, mhc_potential or mj(),
                                 "peptide_mhc", positions, tcr_regions, contact_weight)
    # A position engaged on only one interface contributes nothing on the other, and a NaN there
    # would poison the sum. Absent contacts are a flat zero, not a missing measurement.
    phi_tcr = np.nan_to_num(phi_tcr, nan=0.0)
    phi_mhc = np.nan_to_num(phi_mhc, nan=0.0)
    return ResponseMatrix(
        peptide=peptide, positions=positions, aa=AA20,
        interface_class=tuple("receptor" if p in engaged["tcr_peptide"] else "anchor"
                              for p in positions),
        phi=phi_tcr + phi_mhc, phi_tcr=phi_tcr, phi_mhc=phi_mhc,
        structure_id=contact_map.pdb_id or "",
    )


def _template_peptide(contact_map: ContactMap) -> str:
    """The template's peptide, reconstructed from the contact map's peptide-side residues."""
    seen: dict[int, str] = {}
    for name in INTERFACES:
        side = _PEPTIDE_SIDE[name]
        iface = contact_map.interface(name)
        for p, a in zip(iface[f"pos.{side}"].to_list(), iface[f"residue.aa.{side}"].to_list()):
            seen[int(p)] = a
    if not seen:
        raise ValueError(f"{contact_map.pdb_id} has no peptide contacts; pass `peptide` explicitly")
    length = contact_map.peptide_length or (max(seen) + 1)
    missing = [i + 1 for i in range(length) if i not in seen]
    if missing:
        raise ValueError(
            f"{contact_map.pdb_id}: positions {missing} contact nothing, so the template peptide "
            "cannot be read off the contact map -- pass `peptide` explicitly"
        )
    return "".join(seen[i] for i in range(length))


# ---------------------------------------------------------------- the three queries
def mutation_effect(rm: ResponseMatrix, position: int, aa: str, *,
                    reference: str = "equimolar") -> float:
    """Effect of ONE substitution at ONE position, positive = favourable.

    Args:
        rm: a predicted :class:`ResponseMatrix`.
        position: 1-based peptide position.
        aa: the one-letter residue to put there.
        reference: ``"equimolar"`` (default) scores it against the 1/20 mixture at that position;
            ``"wild_type"`` scores it against the residue the template carries, i.e. the classical
            mutation-scan ``ddG``.

    Returns:
        A single number. Under ``"wild_type"`` the template's own residue returns exactly ``0.0``.
    """
    return float(rm.referenced(reference)[rm.row_of(position), rm.column_of(aa)])


def position_scan(rm: ResponseMatrix, position: int, *,
                  reference: str = "equimolar") -> pl.DataFrame:
    """Effect of EVERY substitution at one position -- one column of the response matrix.

    Args:
        rm: a predicted :class:`ResponseMatrix`.
        position: 1-based peptide position.
        reference: as :func:`mutation_effect`.

    Returns:
        Twenty rows -- ``pos``, ``aa``, ``wt_aa``, ``is_wt``, ``interface_class``, ``phi``,
        ``effect`` -- sorted best residue first, so the head of the frame is what the receptor
        prefers at that position.
    """
    row = rm.row_of(position)
    wt = rm.wild_type_at(position)
    return pl.DataFrame({
        "pos": [position] * 20,
        "aa": list(rm.aa),
        "wt_aa": [wt] * 20,
        "is_wt": [a == wt for a in rm.aa],
        "interface_class": [rm.interface_class[row]] * 20,
        "phi": np.asarray(rm.phi, dtype=float)[row],
        "effect": np.asarray(rm.referenced(reference))[row],
    }).sort("effect", descending=True)


def equimolar_effect(rm: ResponseMatrix, position: int, aa: str | None = None) -> float:
    """Effect of replacing a residue with a random 1/20 mixture at that position.

    This is the sublibrary construction read backwards. The CPL background at position ``i`` is the
    equimolar mixture, so the cost of giving up a defined residue for that mixture is

        ``mean_b Phi(x_{i->b}) - Phi(x_{i->a})``

    -- exactly the ``"equimolar"`` reference. A **positive** result means ``aa`` is better than the
    average residue at that position, so scrambling it to the mixture loses that much; a negative
    result means the mixture is an improvement, i.e. the template residue is disfavoured there.

    Args:
        rm: a predicted :class:`ResponseMatrix`.
        position: 1-based peptide position.
        aa: the residue being given up. Defaults to the one the template carries, which is the usual
            question: how much does this position's identity matter?

    Returns:
        A single number in the potential's units.
    """
    return mutation_effect(rm, position, aa or rm.wild_type_at(position), reference="equimolar")
