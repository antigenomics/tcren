"""Candidate-peptide scoring by amino-acid substitution.

Ports the second half of ``run_TCRen.R``: for each candidate peptide, substitute its
amino acids at the contacted peptide positions of a structure's contact map and sum the
pairwise potential over all contacts. Lower scores indicate more favourable interactions.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import polars as pl

from .contactmap import ContactMap, Interface
from .potential import Potential

# Which side of each interface carries the (substituted) peptide.
_PEPTIDE_SIDE: dict[str, str] = {
    "tcr_peptide": "to",
    "tcr_mhc": "to",  # substitutes the MHC side; peptide is fixed
    "peptide_mhc": "from",
}


# --- The intra-peptide term -------------------------------------------------------------------
# Every interface energy above sums over contacts between two *different* chains, so a peptide held
# in its bound conformation by its own side chains costs the same as one that is not. This term is
# that omission, made explicit and optional: the same pairwise potential summed over the contacts
# the peptide makes with itself (tcren.peptide_internal_contacts).


def _pair_sum(matrix, index, aa_from, aa_to, weights) -> float:
    """Sum ``matrix[aa_from, aa_to] * weights`` over pairs, skipping residues the potential lacks.

    Mirrors :func:`score_peptides`: an amino acid outside the potential's alphabet, or a pair it
    leaves undefined (``NaN``, e.g. cysteine in TCRen), drops out of the sum rather than poisoning it.
    """
    i = np.array([index.get(a, -1) for a in aa_from], dtype=np.int64)
    j = np.array([index.get(b, -1) for b in aa_to], dtype=np.int64)
    valid = (i >= 0) & (j >= 0)
    return float(np.nansum(matrix[i[valid], j[valid]] * weights[valid]))


def _intra_arrays(contact_map: ContactMap, potential: Potential, contact_weight: str):
    """``(symmetrised matrix, aa index, pos.from, pos.to, weights)`` for the peptide's own contacts.

    The potential is symmetrised — ``(F + Fᵀ) / 2`` — because an intra-chain pair has no ``from``/``to``
    orientation to respect. It matters: TCRen is directed (rows read as the TCR side, columns as the
    peptide side), and the canonical residue order that labels a pair's sides is an artefact of the
    contact table, not chemistry. Symmetric potentials such as MJ are unchanged by this.
    """
    pairs = contact_map.peptide_internal
    if pairs is None:
        raise ValueError(
            "the intra-peptide term needs the peptide's own contacts; build the contact map "
            "with ContactMap.from_structure(..., peptide_internal=True)"
        )
    if contact_weight == "atomic":
        if "n_atom_contacts" not in pairs.columns:
            raise ValueError(
                "contact_weight='atomic' needs the n_atom_contacts column; build the contact "
                "map with ContactMap.from_structure(..., count_atoms=True, peptide_internal=True)"
            )
        weights = np.asarray(pairs["n_atom_contacts"].to_list(), dtype=np.float64)
    else:
        weights = np.ones(pairs.height, dtype=np.float64)
    matrix, index = potential.as_matrix()
    return (
        0.5 * (matrix + matrix.T),
        index,
        np.asarray(pairs["pos.from"].to_list(), dtype=np.int64),
        np.asarray(pairs["pos.to"].to_list(), dtype=np.int64),
        weights,
    )


def intra_peptide_energy(
    contact_map: ContactMap,
    potential: Potential,
    peptide: str | None = None,
    contact_weight: str = "residue",
) -> float:
    """The peptide's contact energy **with itself**, the term the interface energies omit.

    Summed over :func:`tcren.peptide_internal_contacts` (4 Å, sequence separation ≥ 3, so sequence
    neighbours — in contact because they are bonded, not because the peptide folded — are excluded),
    under the symmetrised potential. Lower is more favourable, as everywhere in tcren.

    On an extended class-I 9-mer this is a small, sparse term: most such peptides make one or two
    internal contacts, so it moves a score only where a peptide is genuinely bulged or packed
    against itself. That is the point — it is the sequence-dependence the interface sum cannot see.

    Args:
        contact_map: a map built with ``ContactMap.from_structure(..., peptide_internal=True)``.
        potential: pairwise potential (MJ is the sensible default here — TCRen is derived from
            TCR↔peptide contacts, not from a chain's contacts with itself).
        peptide: candidate sequence threaded onto the structure's peptide positions. ``None``
            (default) scores the structure's own residues.
        contact_weight: ``"residue"`` (default, one per contacting pair) or ``"atomic"``
            (weight each pair by its ``n_atom_contacts`` heavy-atom-pair count).

    Returns:
        The summed energy; ``0.0`` when the peptide contacts nothing of itself.
    """
    if contact_weight not in ("residue", "atomic"):
        raise ValueError(f"contact_weight must be 'residue' or 'atomic', got {contact_weight!r}")
    matrix, index, pos_from, pos_to, weights = _intra_arrays(
        contact_map, potential, contact_weight
    )
    if len(weights) == 0:
        return 0.0
    if peptide is None:
        pairs = contact_map.peptide_internal
        aa_from = pairs["residue.aa.from"].to_list()
        aa_to = pairs["residue.aa.to"].to_list()
    else:
        aa_from = [peptide[p] if 0 <= p < len(peptide) else "" for p in pos_from]
        aa_to = [peptide[p] if 0 <= p < len(peptide) else "" for p in pos_to]
    return _pair_sum(matrix, index, aa_from, aa_to, weights)


def score_peptides(
    contact_map: ContactMap,
    candidates: Iterable[str],
    potential: Potential,
    interface: Interface = "tcr_peptide",
    require_same_length: bool = True,
    substituted_side: str | None = None,
    tcr_regions: str = "all",
    contact_weight: str = "residue",
    intra_weight: float = 0.0,
    intra_potential: Potential | None = None,
) -> pl.DataFrame:
    """Score candidate peptides against a structure's contact map.

    Args:
        contact_map: The structure's contact map.
        candidates: Candidate peptide sequences (one-letter).
        potential: Pairwise potential to score with.
        interface: Which interface to score over (default ``"tcr_peptide"``).
        require_same_length: Only score candidates whose length matches the structure's
            peptide length (mirrors the legacy length join). Ignored when the contact
            map has no recorded peptide length.
        substituted_side: ``"to"`` or ``"from"`` — which contact side the candidate is
            threaded onto. Defaults to the peptide side of ``interface``.
        tcr_regions: which TCR regions to keep on the TCR side (``"all"`` default = no
            filter = legacy behaviour; ``"cdr"`` or ``"cdr+fr"`` to restrict).
        contact_weight: ``"residue"`` (default, legacy) gives every contacting residue
            pair unit weight; ``"atomic"`` weights each residue pair by its
            ``n_atom_contacts`` heavy-atom-pair count, so the energy tracks the LJ+Coulomb
            atom-pair sum more closely. ``"atomic"`` requires the contact map to have been
            built with ``count_atoms=True``.
        intra_weight: weight ``w`` of the intra-peptide term, added as
            ``score = E_interface + w * E_intra`` (:func:`intra_peptide_energy`, threaded with the
            same candidate). ``0.0`` (default) leaves the score byte-identical to the interface
            sum. A non-zero weight requires the contact map to have been built with
            ``peptide_internal=True``. The term is on the same energy scale as the interface sum,
            so ``w=1`` treats an internal contact as worth an interface contact.
        intra_potential: potential for the intra-peptide term; defaults to ``potential``. Pass MJ
            when scoring the interface with TCRen, which is derived from TCR↔peptide contacts and
            says nothing about a chain's contacts with itself.

    Returns:
        Columns ``complex.id``, ``peptide``, ``potential``, ``score`` sorted by
        ``complex.id`` then ascending ``score``.
    """
    if contact_weight not in ("residue", "atomic"):
        raise ValueError(f"contact_weight must be 'residue' or 'atomic', got {contact_weight!r}")
    side = substituted_side or _PEPTIDE_SIDE[interface]
    if side not in ("to", "from"):
        raise ValueError(f"substituted_side must be 'to' or 'from', got {side!r}")
    fixed = "from" if side == "to" else "to"

    iface = contact_map.interface(interface, tcr_regions=tcr_regions)
    matrix, index = potential.as_matrix()

    pos_col = iface[f"pos.{side}"]
    if pos_col.null_count():
        raise ValueError(
            f"{pos_col.null_count()} of {len(pos_col)} contacts have no 'pos.{side}': the "
            f"'{side}' chain carries no region markup. Run tcren.annotation.classify_chains("
            "structure) (and tcren.mhc.annotate_mhc for an MHC side) before ContactMap.from_structure."
        )
    pos = np.asarray(pos_col.to_list(), dtype=np.int64)
    fixed_aa = iface[f"residue.aa.{fixed}"].to_list()
    fixed_idx = np.array([index.get(a, -1) for a in fixed_aa], dtype=np.int64)

    if contact_weight == "atomic":
        if "n_atom_contacts" not in iface.columns:
            raise ValueError(
                "contact_weight='atomic' needs the n_atom_contacts column; build the "
                "contact map with ContactMap.from_structure(..., count_atoms=True)"
            )
        weights = np.asarray(iface["n_atom_contacts"].to_list(), dtype=np.float64)
    else:
        weights = np.ones(len(pos), dtype=np.float64)

    # The intra-peptide arrays do not depend on the candidate, so they are gathered once.
    intra = (
        _intra_arrays(contact_map, intra_potential or potential, contact_weight)
        if intra_weight
        else None
    )

    candidates = list(candidates)
    rows = []
    for peptide in candidates:
        if require_same_length and contact_map.peptide_length is not None:
            if len(peptide) != contact_map.peptide_length:
                continue
        # Gather the substituted amino acid for each contact from the candidate.
        subst_idx = np.array(
            [index.get(peptide[p], -1) if 0 <= p < len(peptide) else -1 for p in pos],
            dtype=np.int64,
        )
        if side == "to":
            rows_idx, cols_idx = fixed_idx, subst_idx
        else:
            rows_idx, cols_idx = subst_idx, fixed_idx
        valid = (rows_idx >= 0) & (cols_idx >= 0)
        vals = matrix[rows_idx[valid], cols_idx[valid]] * weights[valid]
        # Pairs absent from the potential (e.g. Cys on the 'from' axis) are dropped,
        # exactly as the inner join in run_TCRen.R drops unmatched rows.
        score = float(np.nansum(vals))
        if intra is not None:
            i_matrix, i_index, i_from, i_to, i_weights = intra
            score += intra_weight * _pair_sum(
                i_matrix, i_index,
                [peptide[p] if 0 <= p < len(peptide) else "" for p in i_from],
                [peptide[p] if 0 <= p < len(peptide) else "" for p in i_to],
                i_weights,
            )
        rows.append({"complex.id": contact_map.pdb_id, "peptide": peptide, "score": score})

    out = pl.DataFrame(
        rows,
        schema={"complex.id": pl.Utf8, "peptide": pl.Utf8, "score": pl.Float64},
    ).with_columns(pl.lit(potential.name).alias("potential"))
    return out.select("complex.id", "peptide", "potential", "score").sort(
        "complex.id", "score"
    )


def score_structures(
    contact_maps: Iterable[ContactMap],
    candidates: Iterable[str],
    potential: Potential,
    **kwargs,
) -> pl.DataFrame:
    """Score candidates against several structures and stack the results."""
    candidates = list(candidates)
    frames = [score_peptides(cm, candidates, potential, **kwargs) for cm in contact_maps]
    return pl.concat(frames) if frames else pl.DataFrame()


# --- Recognition matrix: the per-position amino-acid preference landscape from F ---------------
# Generalises the CPL positional-scan matrix to either side of an interface. Substituting the peptide
# side gives the CPL-matrix analog (position x AA over the epitope); substituting the TCR side gives
# the motif analog (position x AA over the CDR3) -- a physics-derived recognition matrix.

from dataclasses import dataclass as _dataclass  # noqa: E402

_AA20 = tuple("ACDEFGHIKLMNPQRSTVWY")


@_dataclass(slots=True)
class RecognitionMatrix:
    """Per-position × amino-acid substitution-energy landscape (see :func:`recognition_matrix`).

    ``energy[i, a]`` is the summed pairwise potential over position ``positions[i]``'s contacts when
    amino acid ``aa[a]`` sits there and the other side is held fixed. **Lower = more favourable**, so a
    per-position preference is ``-energy`` (higher = preferred). Positions with no contact are omitted.
    Entries are ``NaN`` for amino acids the potential leaves undefined (e.g. cysteine pairs in TCRen),
    exactly as :func:`score_peptides` drops those contacts — so reduce columns with ``np.nan*`` ops.
    """

    positions: list          #: one ``(chain_type, region, pos, native_aa)`` tuple per row of ``energy``
    aa: tuple                #: the 20 amino-acid column order
    energy: object           #: ``(n_positions, 20)`` float ndarray of substitution energies
    side: str                #: which side was scanned (``"from"`` = TCR, ``"to"`` = peptide)
    interface: str


def recognition_matrix(
    contact_map: ContactMap,
    potential: Potential,
    *,
    interface: Interface = "tcr_peptide",
    side: str | None = None,
    tcr_regions: str = "all",
) -> RecognitionMatrix:
    """The per-position × 20-AA substitution-energy matrix for one interface side.

    For ``interface="tcr_peptide"``, ``side="from"`` scans the **TCR/CDR3** (the motif-matrix analog)
    and ``side="to"`` scans the **peptide** (the CPL-matrix analog). Each entry is the F energy summed
    over that position's contacts with the given amino acid substituted in, the opposite side fixed —
    the same virtual-substitution path as :func:`score_peptides`, resolved per position rather than
    summed over the whole sequence.

    Args:
        contact_map: the structure's contact map.
        potential: pairwise potential (TCRen for TCR:peptide).
        interface: which interface to score over.
        side: ``"from"`` or ``"to"``; defaults to the **non**-peptide side for ``tcr_peptide``
            (i.e. the TCR), and to the peptide side for the presentation interfaces.
        tcr_regions: TCR-region filter (``"all"``/``"cdr"``/``"cdr+fr"``) — use ``"cdr"`` to restrict
            a TCR-side scan to the CDRs.

    Returns:
        A :class:`RecognitionMatrix`. Rows are the contacted positions in ``(chain, region, pos)``
        order; columns are :data:`RecognitionMatrix.aa`.
    """
    if side is None:
        side = "from" if interface == "tcr_peptide" else _PEPTIDE_SIDE[interface]
    if side not in ("from", "to"):
        raise ValueError(f"side must be 'from' or 'to', got {side!r}")
    fixed = "to" if side == "from" else "from"
    iface = contact_map.interface(interface, tcr_regions=tcr_regions)
    matrix, index = potential.as_matrix()
    aa_idx = np.array([index.get(a, -1) for a in _AA20], dtype=np.int64)

    fixed_aa = iface[f"residue.aa.{fixed}"].to_list()
    fixed_idx = np.array([index.get(a, -1) for a in fixed_aa], dtype=np.int64)
    keys = list(zip(iface[f"chain.type.{side}"].to_list(), iface[f"region.type.{side}"].to_list(),
                    iface[f"pos.{side}"].to_list(), iface[f"residue.aa.{side}"].to_list()))

    order: list = []
    rows: list = []
    seen: dict = {}
    for row_i, key in enumerate(keys):
        if key not in seen:
            seen[key] = len(order)
            order.append(key)
            rows.append(np.zeros(20))
        fj = fixed_idx[row_i]
        if fj < 0:
            continue
        for a in range(20):
            ai = aa_idx[a]
            if ai < 0:
                continue
            # potential is directed [from, to]; substitute on `side`, keep `fixed`.
            rows[seen[key]][a] += matrix[ai, fj] if side == "from" else matrix[fj, ai]
    return RecognitionMatrix(positions=order, aa=_AA20,
                             energy=np.vstack(rows) if rows else np.zeros((0, 20)),
                             side=side, interface=interface)
