"""Candidate-peptide scoring by amino-acid substitution.

Ports the second half of ``run_TCRen.R``: for each candidate peptide, substitute its
amino acids at the contacted peptide positions of a structure's contact map and sum the
pairwise potential over all contacts. Lower scores indicate more favourable interactions.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import polars as pl

from ..contactmap import ContactMap, Interface
from ..potential import Potential, mj

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


def intra_peptide_energy(
    contact_map: ContactMap,
    potential: Potential,
    peptide: str | None = None,
    contact_weight: str = "residue",
) -> float:
    """The peptide's contact energy **with itself**, the term the interface energies omit.

    Summed over :func:`tcren.peptide_internal_contacts` (5 Å, sequence separation ≥ 3, so sequence
    neighbours — in contact because they are bonded, not because the peptide folded — are excluded),
    under the **symmetrised** potential ``(F + Fᵀ) / 2``: an intra-chain pair has no ``from``/``to``
    orientation to respect, and the canonical residue order labelling a pair's sides is an artefact
    of the contact table, not chemistry. That matters for a directed potential such as TCRen and is
    a no-op for a symmetric one such as MJ. Lower is more favourable, as everywhere in tcren.

    On an extended class-I 9-mer this is a small, sparse term: such peptides make one or two internal
    contacts, so it moves a score only where a peptide is genuinely bulged or packed against itself.
    That is the point — it is the sequence-dependence the interface sum cannot see.

    Args:
        contact_map: a map built with ``ContactMap.from_structure(..., peptide_internal=True)``.
        potential: pairwise potential (MJ is the sensible choice — TCRen is derived from
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
    pairs = contact_map.peptide_internal
    if pairs is None:
        raise ValueError(
            "the intra-peptide term needs the peptide's own contacts; build the contact map "
            "with ContactMap.from_structure(..., peptide_internal=True)"
        )
    if contact_weight == "atomic" and "n_atom_contacts" not in pairs.columns:
        raise ValueError(
            "contact_weight='atomic' needs the n_atom_contacts column; build the contact "
            "map with ContactMap.from_structure(..., count_atoms=True, peptide_internal=True)"
        )
    if pairs.is_empty():
        return 0.0
    weights = (np.asarray(pairs["n_atom_contacts"].to_list(), dtype=np.float64)
               if contact_weight == "atomic" else np.ones(pairs.height))
    if peptide is None:
        aa_from, aa_to = pairs["residue.aa.from"].to_list(), pairs["residue.aa.to"].to_list()
    else:
        aa_from, aa_to = ([peptide[p] if 0 <= p < len(peptide) else "" for p in side]
                          for side in (pairs["pos.from"], pairs["pos.to"]))
    matrix, index = potential.as_matrix()
    # Pairs whose residue is outside the alphabet, or that the potential leaves undefined (NaN,
    # e.g. cysteine in TCRen), drop out of the sum — exactly as score_peptides drops them.
    i = np.array([index.get(a, -1) for a in aa_from], dtype=np.int64)
    j = np.array([index.get(b, -1) for b in aa_to], dtype=np.int64)
    valid = (i >= 0) & (j >= 0)
    sym = 0.5 * (matrix + matrix.T)
    return float(np.nansum(sym[i[valid], j[valid]] * weights[valid]))


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
    weights: "np.ndarray | None" = None,
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
        intra_potential: potential for the intra-peptide term; defaults to **MJ**, not to
            ``potential`` — TCRen is derived from TCR↔peptide contacts and says nothing about a
            chain's contacts with itself.
        weights: an explicit per-contact multiplier applied on top of ``contact_weight``, one value
            per row of the selected interface and in its row order. This is how a rotamer-averaged
            contact probability (:func:`tcren.rotamers.contact_probabilities`), a per-position
            weight (:func:`position_weights`) or a contact-type filter enters the sum. ``None``
            (default) leaves the score byte-identical.

    Returns:
        Columns ``complex.id``, ``peptide``, ``potential``, ``score`` sorted by
        ``complex.id`` then ascending ``score``.
    """
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


    weights = _contact_weights(iface, contact_weight, weights)

    if intra_weight:
        intra_potential = intra_potential or mj()

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
        if intra_weight:
            score += intra_weight * intra_peptide_energy(
                contact_map, intra_potential, peptide, contact_weight
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
    and ``side="to"`` scans the **peptide** (the CPL-matrix analog). Each entry is the Φ energy summed
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


# =========================================================================================
# peptide position: role annotation, weighting, and the per-position energy profile
# =========================================================================================
#: Named per-position weighting schemes for :func:`position_weights`.
POSITION_SCHEMES = ("uniform", "central", "tcr_facing")


def peptide_positions(contact_map: ContactMap, structure=None, interface: Interface = "tcr_peptide",
                      tcr_regions: str = "all") -> pl.DataFrame:
    """Annotate an interface's contacts with the peptide position and role they involve.

    The position was always there — ``pos.to`` on the ``tcr_peptide`` interface is the 0-based
    peptide index, because the peptide chain carries one full-length region starting at 0 — and
    :mod:`tcren.refine.anchors` has always predicted anchors. The two were never joined, so nothing
    downstream could ask whether a contact sits on an anchor or in the TCR-facing bulge.

    Args:
        contact_map: the structure's contact map.
        structure: the source structure. Passed to :func:`tcren.refine.predict_anchors`, which then
            uses the real MHC-class call rather than the peptide-length heuristic. Recommended for
            class II, where a 12-20mer would otherwise be misread as class I.
        interface: which interface (must have a peptide side).
        tcr_regions: TCR-region filter, passed through to :meth:`ContactMap.interface`.

    Returns:
        The interface frame plus ``peptide.pos`` (1-based P-number), ``peptide.aa`` and
        ``peptide.role`` (``"anchor"`` or ``"tcr_facing"``).

    Raises:
        ValueError: if the peptide side carries no region markup (null positions).
    """
    from ..refine.anchors import predict_anchors

    side = _PEPTIDE_SIDE[interface]
    iface = contact_map.interface(interface, tcr_regions=tcr_regions)
    if iface.height == 0:
        return iface.with_columns(pl.lit(None, dtype=pl.Int64).alias("peptide.pos"),
                                  pl.lit(None, dtype=pl.Utf8).alias("peptide.aa"),
                                  pl.lit(None, dtype=pl.Utf8).alias("peptide.role"))
    pos_col = iface[f"pos.{side}"]
    if pos_col.null_count():
        raise ValueError(f"the peptide side carries no region markup; {pos_col.null_count()} of "
                         f"{len(pos_col)} contacts have a null 'pos.{side}'")

    # Take the sequence from the structure when we have it. Reassembling it from the contacts
    # leaves 'X' wherever the TCR touches nothing, and the class-II register heuristic slides a
    # 9-mer window over the sequence — on a peptide that is mostly X it picks a register from
    # almost no evidence.
    peptide = None
    if structure is not None:
        from ..refine.anchors import native_peptide
        try:
            peptide = native_peptide(structure)
        except (ValueError, KeyError):
            peptide = None
    if peptide is None:
        aa_by_pos = dict(zip(pos_col.to_list(), iface[f"residue.aa.{side}"].to_list()))
        length = contact_map.peptide_length or (max(aa_by_pos) + 1)
        peptide = "".join(aa_by_pos.get(i, "X") for i in range(length))
    anchors = set(predict_anchors(peptide, structure).anchors)

    pos = np.asarray(pos_col.to_list(), dtype=np.int64)
    return iface.with_columns(
        pl.Series("peptide.pos", pos + 1),
        pl.Series("peptide.aa", [peptide[p] if 0 <= p < len(peptide) else "X" for p in pos]),
        pl.Series("peptide.role", ["anchor" if p in anchors else "tcr_facing" for p in pos]),
    )


def position_weights(annotated: pl.DataFrame, scheme: str = "uniform",
                     length: int | None = None) -> np.ndarray:
    """Per-contact weights from where along the peptide each contact sits.

    A contact potential sums every contact alike, so a clash at an anchor — which the groove
    tolerates and a TCR never touches — costs the same as one under the CDR3 loops. These schemes
    let the sum say otherwise; feed the result to ``score_peptides(..., weights=...)``.

    Args:
        annotated: the frame :func:`peptide_positions` returns.
        scheme: ``"uniform"`` (all ones — the default everywhere, so nothing moves unless asked),
            ``"central"`` (triangular in ``peptide.pos``, peaking at the middle of the peptide and
            falling to 0 at either terminus), or ``"tcr_facing"`` (1 off the anchors, 0 on them).
        length: peptide length for the ``"central"`` ramp; taken from the annotation when omitted.

    Returns:
        One float per row of ``annotated``, in its row order.

    Raises:
        ValueError: for an unknown ``scheme``.
    """
    if scheme not in POSITION_SCHEMES:
        raise ValueError(f"scheme must be one of {POSITION_SCHEMES}, got {scheme!r}")
    n = annotated.height
    if scheme == "uniform" or n == 0:
        return np.ones(n, dtype=np.float64)
    if scheme == "tcr_facing":
        return (np.asarray(annotated["peptide.role"].to_list()) == "tcr_facing").astype(np.float64)

    pos = np.asarray(annotated["peptide.pos"].to_list(), dtype=np.float64)
    length = float(length or annotated["peptide.pos"].max())
    centre = (length + 1.0) / 2.0
    half = max(centre - 1.0, 1.0)
    return np.clip(1.0 - np.abs(pos - centre) / half, 0.0, 1.0)


def position_profile(contact_map: ContactMap, potential: Potential, structure=None,
                     interface: Interface = "tcr_peptide", tcr_regions: str = "all",
                     contact_weight: str = "residue") -> pl.DataFrame:
    """Per-peptide-position decomposition of the interface energy.

    The sum :func:`score_peptides` reports, resolved along the peptide instead of collapsed: which
    positions carry the interaction, and which carry strain. Summing ``phi`` reproduces the total.

    Args:
        contact_map: the structure's contact map.
        potential: the pairwise potential (TCRen for TCR:peptide).
        structure: source structure, for the anchor call (see :func:`peptide_positions`).
        interface: which interface.
        tcr_regions: TCR-region filter.
        contact_weight: ``"residue"`` or ``"atomic"``, as elsewhere.

    Returns:
        One row per contacted position: ``complex.id``, ``peptide.pos``, ``peptide.aa``,
        ``peptide.role``, ``n_contacts``, ``phi``.
    """

    ann = peptide_positions(contact_map, structure, interface, tcr_regions)
    if ann.height == 0:
        return pl.DataFrame(schema={"complex.id": pl.Utf8, "peptide.pos": pl.Int64,
                                    "peptide.aa": pl.Utf8, "peptide.role": pl.Utf8,
                                    "n_contacts": pl.UInt32, "phi": pl.Float64})
    w = _contact_weights(ann, contact_weight)
    matrix, index = potential.as_matrix()
    i = np.array([index.get(a, -1) for a in ann["residue.aa.from"].to_list()], dtype=np.int64)
    j = np.array([index.get(b, -1) for b in ann["residue.aa.to"].to_list()], dtype=np.int64)
    e = np.where((i >= 0) & (j >= 0), matrix[np.clip(i, 0, None), np.clip(j, 0, None)], np.nan) * w

    return (ann.with_columns(pl.Series("phi", np.nan_to_num(e, nan=0.0)))
            .group_by("peptide.pos", "peptide.aa", "peptide.role", maintain_order=True)
            .agg(pl.len().alias("n_contacts"), pl.col("phi").sum())
            .with_columns(pl.lit(contact_map.pdb_id).alias("complex.id"))
            .select("complex.id", "peptide.pos", "peptide.aa", "peptide.role", "n_contacts", "phi")
            .sort("peptide.pos"))


def central_strain(profile: pl.DataFrame, band: float = 1 / 3) -> float:
    """Interface energy carried by the peptide's central, TCR-facing band.

    The review's concern, made a number: a TCR has to clear the middle of the peptide to dock at
    all, so an unfavourable (positive) energy there is a viability question in a way that the same
    value at P1 or PΩ is not. Positive = the centre is repulsive.

    Args:
        profile: the frame :func:`position_profile` returns.
        band: fraction of the peptide's length counted as central, centred on the middle.

    Returns:
        Summed ``phi`` over the central band, or ``nan`` for an empty profile.
    """
    if profile.height == 0:
        return float("nan")
    pos = np.asarray(profile["peptide.pos"].to_list(), dtype=np.float64)
    lo, hi = pos.min(), pos.max()
    centre, half = (lo + hi) / 2.0, max((hi - lo + 1) * band / 2.0, 0.5)
    sel = np.abs(pos - centre) <= half
    return float(np.asarray(profile["phi"].to_list())[sel].sum())


# --- the interface energy sum, and the two things it needs -----------------------------------------
# These were in `tcren.pipeline`, which is the top-of-stack runner. Four modules below it in the
# stack -- `ddg`, `rotamers`, this one, and the descriptor dispatch -- each reached up for one of
# them with a function-local import, because summing a potential over a contact table is not a
# pipeline concern, it is the energy concern this module already owns. `pipeline` now imports them
# from here, which is the direction the dependency was always meant to run.

def _phi_scale(interface: str, potential: Potential) -> float:
    """The scale that makes one interface energy commensurate with the other two.

    The Native2026 standard deviation of that interface's energy under that potential, read from
    the frozen moments; :meth:`Potential.scale` (the sd of the potential's own matrix) when the
    pair is not tabulated, which is what an unbundled potential or a non-default assignment gets.
    """
    from ..reliability import moments

    key = f"{interface}|{potential.name}"
    entry = moments().get("phi", {}).get(key)
    return float(entry["sd"]) if entry else potential.scale()


def _contact_weights(contacts: pl.DataFrame, contact_weight: str = "residue",
                     weights: "np.ndarray | None" = None) -> np.ndarray:
    """Per-contact multiplier for an energy sum.

    Every score in the package is ``sum_ij w_ij * e(a_i, b_j)``; this is the only place ``w``
    comes from, so a rotamer-averaged contact probability
    (:func:`tcren.rotamers.contact_probabilities`), a position weight
    (:func:`tcren.scoring.position_weights`) or a contact-type filter all enter the same way.

    Args:
        contacts: the interface frame the energy is summed over.
        contact_weight: ``"residue"`` (unit weight per contacting residue pair) or ``"atomic"``
            (its ``n_atom_contacts`` heavy-atom-pair count).
        weights: an explicit per-row multiplier, applied **on top of** ``contact_weight``. Must
            be one value per row.

    Raises:
        ValueError: for an unknown ``contact_weight``, a missing ``n_atom_contacts`` column, or a
            ``weights`` array of the wrong length.
    """
    if contact_weight not in ("residue", "atomic"):
        raise ValueError(f"contact_weight must be 'residue' or 'atomic', got {contact_weight!r}")
    if contact_weight == "atomic":
        if "n_atom_contacts" not in contacts.columns:
            raise ValueError(
                "contact_weight='atomic' needs the n_atom_contacts column; build the "
                "contact map with count_atoms=True"
            )
        out = np.asarray(contacts["n_atom_contacts"].to_list(), dtype=np.float64)
    else:
        out = np.ones(contacts.height, dtype=np.float64)
    if weights is not None:
        weights = np.asarray(weights, dtype=np.float64)
        if weights.shape != (contacts.height,):
            raise ValueError(f"weights must have one value per contact "
                             f"({contacts.height}), got {weights.shape}")
        out = out * weights
    return out


def _interface_energy(
    contacts: pl.DataFrame, potential: Potential, contact_weight: str = "residue",
    weights: "np.ndarray | None" = None,
) -> float:
    """Sum the residue-pair ``potential`` over an interface's contacts (unknown residues skipped).

    With ``contact_weight="residue"`` (default, legacy) each contacting residue pair adds
    ``potential[a, b]``. With ``contact_weight="atomic"`` each pair is multiplied by its
    ``n_atom_contacts`` heavy-atom-pair count (which the contacts table must carry). ``weights``
    multiplies on top — see :func:`_contact_weights`.
    """
    if contacts.is_empty():
        return 0.0
    weights = _contact_weights(contacts, contact_weight, weights)
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
