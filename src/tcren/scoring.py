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


def score_peptides(
    contact_map: ContactMap,
    candidates: Iterable[str],
    potential: Potential,
    interface: Interface = "tcr_peptide",
    require_same_length: bool = True,
    substituted_side: str | None = None,
    tcr_regions: str = "all",
    contact_weight: str = "residue",
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

    pos = np.asarray(iface[f"pos.{side}"].to_list(), dtype=np.int64)
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
