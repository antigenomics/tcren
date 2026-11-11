"""Closed-form CDR3 Ω-loop reachability — the feasibility filter.

The CDR3 loop is a chain of ``n = L + 1`` virtual Cα–Cα bonds of length :data:`BOND_CA`, pinned at
its neck separation ``R`` (the ``cdr3{a,b}_ext`` frame descriptor). Pure geometry then bounds how far
the loop's apex can reach from the neck axis, independent of sequence or energy:

* :func:`reach_max` — the maximum apex reach for a loop of length ``L`` pinned at ``R``.
* :func:`reachability_floor` — the shortest loop that can span a target at distance ``d``.
* :func:`span_saturation` — the fraction of maximal reach a structure actually uses,
  ``reach / reach_max``; ``> 1`` is geometrically impossible and real repertoires sit well below 1
  (~0.3–0.5 at CDR3 length 10–20), i.e. with large conformational slack.

These are **feasibility / mode** descriptors, not binder classifiers — backbone Cα geometry does not
discriminate binders (every apparent signal was a provenance / epitope confound). Bill them as
``structure -> geometry`` featurisation only.

"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = ["BOND_CA", "CB_OFFSET", "enforce_bond", "graphon_weights", "orientation_cosine",
           "read_back", "reach_max", "reachability_floor", "resample_graphon", "span_saturation",
           "cdr3_internal_coords", "virtual_cb", "LoopInternalCoords"]

#: Graphon conventions :func:`graphon_weights` and :func:`read_back` understand. They are not
#: interchangeable; which pairing to use is measured, not deduced --- see :func:`graphon_weights`.
CONVENTIONS = ("interior", "endpoint")

#: Virtual Cα–Cα bond length (Å) used for the loop chain.
BOND_CA = 3.80

#: Cβ displacement from Cα in the local three-Cα frame of :func:`virtual_cb`, ``(u, v, w)`` in Å.
#: Least squares over 4,346 CDR residues of the ``Canonical2026`` crystals; the residual is 0.547 Å
#: mean and 0.412 Å median.
CB_OFFSET = (0.0803, -1.0313, -0.7252)


def virtual_cb(ca: np.ndarray) -> np.ndarray:
    """Cβ positions built from a Cα trace alone, ``(..., n, 3) -> (..., n, 3)``.

    A drawn geometry has no side chains, so the question *does this residue turn toward the
    peptide* cannot be asked of it --- unless Cβ is reconstructed from the backbone, which it can
    be: three consecutive Cα atoms fix a local frame and Cβ sits at a nearly constant offset in it.

    Fitted on the ``Canonical2026`` CDR residues, the position residual is **0.547 Å mean, 0.412 Å
    median**, and --- what actually matters --- the resulting :func:`orientation_cosine` correlates
    with the one from the real Cβ at **Pearson 0.81**, separating contacts within a fixed Cα-distance
    bin at ROC-AUC 0.68-0.79 against the real side chain's 0.70-0.84 (50,446 pairs, 6,603 contacts).

    Glycine is not a special case: this is where a Cβ *would* sit given the backbone, which is the
    quantity the orientation term wants, so there is no missing value to encode.

    Args:
        ca: ``(..., n, 3)`` Cα coordinates, ``n >= 3``.

    Returns:
        ``(..., n, 3)`` Cβ positions. The two termini reuse their neighbour's frame, having no
        three-Cα frame of their own; three collinear Cα fall back to an arbitrary frame orthogonal
        to the chain rather than returning ``NaN``.

    Raises:
        ValueError: if the trace has fewer than 3 residues or the last axis is not length 3.
    """
    arr = np.asarray(ca, float)
    if arr.shape[-1] != 3:
        raise ValueError(f"expected (..., n, 3) coordinates, got shape {arr.shape}")
    if arr.shape[-2] < 3:
        raise ValueError(f"a three-Cα frame needs at least 3 residues, got {arr.shape[-2]}")
    prev, mid, nxt = arr[..., :-2, :], arr[..., 1:-1, :], arr[..., 2:, :]
    u = nxt - prev
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    v = 0.5 * (prev + nxt) - mid
    v -= (v * u).sum(-1, keepdims=True) * u
    # Three collinear Cα leave no bisector. A real chain never is, but a synthetic or a heavily
    # extended one can be, and a NaN here would propagate silently into every cosine.
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    fallback = np.cross(u, np.where(np.abs(u[..., :1]) < 0.9, np.array([1.0, 0.0, 0.0]),
                                    np.array([0.0, 1.0, 0.0])))
    v = np.where(norm > 1e-8, v / np.where(norm > 1e-8, norm, 1.0),
                 fallback / np.linalg.norm(fallback, axis=-1, keepdims=True))
    w = np.cross(u, v)
    a, b, c = CB_OFFSET
    out = np.empty_like(arr)
    out[..., 1:-1, :] = mid + a * u + b * v + c * w
    out[..., 0, :] = arr[..., 0, :] + (out[..., 1, :] - arr[..., 1, :])
    out[..., -1, :] = arr[..., -1, :] + (out[..., -2, :] - arr[..., -2, :])
    return out


def orientation_cosine(ca: np.ndarray, partner_ca: np.ndarray,
                       cb: np.ndarray | None = None) -> np.ndarray:
    """Does the residue turn toward the partner? ``(..., n, m)`` cosines in ``[-1, 1]``.

    The cosine between a residue's Cα→Cβ direction and its Cα→partner-Cα direction: ``+1`` is a side
    chain aimed straight at the partner residue, ``-1`` pointed away. **This is where the energy
    is** --- the Cα geometry says a pair is close enough, and this says whether the residue is
    facing the right way to use it. Held at a fixed Cα distance it separates contact from
    non-contact at ROC-AUC 0.69-0.84 across the 4-12 Å range.

    Args:
        ca: ``(..., n, 3)`` receptor Cα.
        partner_ca: ``(..., m, 3)`` partner Cα.
        cb: ``(..., n, 3)`` Cβ. Defaults to :func:`virtual_cb` of ``ca``, which is what a drawn
            geometry has; pass the real Cβ when scoring a solved structure.

    Returns:
        ``(..., n, m)`` cosines. Zero where the Cα and Cβ coincide.
    """
    ca = np.asarray(ca, float)
    partner_ca = np.asarray(partner_ca, float)
    side = (virtual_cb(ca) if cb is None else np.asarray(cb, float)) - ca
    to_partner = partner_ca[..., None, :, :] - ca[..., :, None, :]
    side_norm = np.linalg.norm(side, axis=-1)
    partner_norm = np.linalg.norm(to_partner, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos = np.einsum("...ik,...ijk->...ij", side, to_partner) / (
            side_norm[..., :, None] * partner_norm)
    return np.where(side_norm[..., :, None] > 1e-6, np.nan_to_num(cos), 0.0)


def graphon_weights(n: int, grid: int, *, convention: str = "interior") -> np.ndarray:
    r"""``(grid, n)`` bilinear weights carrying ``n`` residues onto ``grid`` graphon cells.

Registration and read-back are the two directions of one map, and the geometry prior does
    **not** pair them adjointly: the archived geometry prior registers residue :math:`i` at
    grid coordinate :math:`i(G-1)/(L-1)` while its read-back reads it
    back from :math:`(i+1)(G-1)/(L+1)`. A round trip therefore contracts a loop's footprint toward
    its midpoint by :math:`(L-1)/(L+1)` --- 0.714 for a 6-mer against 0.889 for a 17-mer, zero error
    at the midpoint, 1.29 Å at the termini of a 6-mer, and **a fixed bias no grid refinement
    removes**.

    **That mismatch is deliberate to keep, and the reason is measured** (2026-09-07, 7,732 observed
    contacts over 370 leave-one-out ``Native2026`` crystals at grid 32; predicted minus measured
    Cα distance, in Å):

    ==========================  =========  =========  =========
    registration / read-back    CDR2       CDR3       all
    ==========================  =========  =========  =========
    endpoint / interior (used)  +5.76      +2.29      **+2.91**
    endpoint / endpoint         +5.57      +2.78      +3.26
    interior / interior         +5.63      +2.79      +3.28
    ==========================  =========  =========  =========

    An adjoint pairing buys 0.19 Å on CDR2 and costs 0.49 Å on CDR3. Two things make the mismatched
    pairing the better one, and neither is the round trip: ``"endpoint"`` registration lands residue
    0 on cell 0 and residue :math:`L-1` on cell :math:`G-1`, so loop **termini align exactly across
    structures of different length**, which is the operation an average over crystals actually
    performs; and the ``"interior"`` read-back samples away from the surface's extrapolated edge
    cells, a shrinkage toward the loop centre that helps because the average over-estimates
    distance everywhere. **So the CDR2 deficit is not a registration artefact** --- it survives
    every pairing at 5.6-5.8 Å against CDR3's 2.3-2.8 Å.

    ``convention`` names which map:

    ``"interior"``
        Cell ``g`` samples residue coordinate :math:`g(n+1)/(G-1) - 1`, the inverse of the
        :math:`(i+1)/(n+1)` read-back. Adjoint to ``read_back(..., convention="interior")``, whose
        round trip then converges to the identity as the grid is refined.
    ``"endpoint"``
        Cell ``g`` samples ``np.linspace(0, n - 1, grid)[g]``. What ``_resample`` and
        ``contactmap._register`` have always done, so :func:`tcren.contactmap.registered_map` and
        the archived geometry priors stay reproducible. Adjoint to
        ``read_back(..., convention="endpoint")``, with the same endpoint coordinate convention. Interpolation can still lose detail.

    Out-of-range coordinates clamp, matching ``scipy.ndimage.map_coordinates(mode="nearest")``; the
    matrix form agrees with it to 7e-15 Å and, unlike it, applies to a whole stack in one ``@``.

    Args:
        n: Residues in the chain being registered.
        grid: Graphon cells per axis, ``G``.
        convention: ``"interior"`` or ``"endpoint"``.

    Returns:
        ``(grid, n)`` array ``W`` with ``W @ f`` the resampled signal.

    Raises:
        ValueError: on an unknown ``convention``, ``n < 1``, or ``grid < 2`` under ``"interior"``
            (which needs two cells to define the coordinate).

    Example:
        >>> graphon_weights(4, 4, convention="endpoint").round(3)
        array([[1.   , 0.   , 0.   , 0.   ],
               [0.   , 1.   , 0.   , 0.   ],
               [0.   , 0.   , 1.   , 0.   ],
               [0.   , 0.   , 0.   , 1.   ]])
    """
    if convention not in CONVENTIONS:
        raise ValueError(f"convention must be one of {CONVENTIONS}, got {convention!r}")
    if n < 1:
        raise ValueError(f"need at least one residue, got n={n}")
    if convention == "endpoint":
        coord = np.linspace(0.0, n - 1, grid)
    else:
        if grid < 2:
            raise ValueError(f"the 'interior' convention needs grid >= 2, got {grid}")
        coord = np.arange(grid) * (n + 1) / (grid - 1) - 1.0
    return _bilinear(coord, n)


def _bilinear(coord: np.ndarray, size: int) -> np.ndarray:
    """``(len(coord), size)`` linear-interpolation weights, clamped like ``mode="nearest"``."""
    coord = np.clip(np.asarray(coord, float), 0.0, size - 1)
    lo = np.floor(coord).astype(int)
    hi = np.minimum(lo + 1, size - 1)
    frac = coord - lo
    w = np.zeros((len(coord), size))
    rows = np.arange(len(coord))
    w[rows, lo] = 1.0 - frac
    w[rows, hi] += frac
    return w


def _read_weights(n: int, grid: int, convention: str = "interior") -> np.ndarray:
    r"""``(n, grid)`` weights reading ``n`` residues out of a ``grid``-cell axis.

    ``"interior"`` samples :math:`(i+1)/(n+1)`; ``"endpoint"`` samples :math:`i/(n-1)`, putting the
    first and last residue exactly on the first and last cell. Each is adjoint to the
    :func:`graphon_weights` of the same name.
    """
    if convention == "endpoint":
        coord = np.linspace(0.0, grid - 1, n)
    elif convention == "interior":
        coord = (np.arange(1, n + 1) / (n + 1)) * (grid - 1)
    else:
        raise ValueError(f"convention must be one of {CONVENTIONS}, got {convention!r}")
    return _bilinear(coord, grid)


def resample_graphon(dmap: np.ndarray, grid: int, *, convention: str = "interior") -> np.ndarray:
    """Resample an ``(..., L, M)`` map onto ``(..., grid, grid)`` graphon cells.

    Args:
        dmap: Map to register; the last two axes are the receptor and partner indices, and any
            leading axes are carried through untouched.
        grid: Graphon cells per axis.
        convention: Passed to :func:`graphon_weights`.

    Returns:
        The registered ``(..., grid, grid)`` map.
    """
    a = graphon_weights(dmap.shape[-2], grid, convention=convention)
    b = graphon_weights(dmap.shape[-1], grid, convention=convention)
    return np.einsum("gl,...lm,hm->...gh", a, np.asarray(dmap, float), b, optimize=True)


def read_back(surface: np.ndarray, loop_len: int, pep_len: int, *,
              convention: str = "interior") -> np.ndarray:
    r"""Read a ``(..., G, G)`` graphon surface back at ``(loop_len, pep_len)``.

    The surface is sampled at :math:`(u, v) = ((i+1)/(L+1), (j+1)/(M+1))`, which is what carries a
    shape averaged over one set of lengths to a pair of lengths it never saw. Leading axes are
    carried through, so a whole member stack reads back in one call rather than one call per member.

    Args:
        surface: ``(..., G, G)`` registered map.
        loop_len: Receptor residues ``L`` to read back at.
        pep_len: Partner residues ``M`` to read back at.
        convention: Sampling positions; adjoint to the :func:`graphon_weights` of the same name.
            ``"interior"`` is what the geometry prior has always read back at.

    Returns:
        ``(..., loop_len, pep_len)``.
    """
    grid = surface.shape[-1]
    a = _read_weights(loop_len, grid, convention)        # (L, G)
    b = _read_weights(pep_len, grid, convention)         # (M, G)
    return np.einsum("lg,...gh,mh->...lm", a, np.asarray(surface, float), b, optimize=True)


def enforce_bond(ca: np.ndarray, bond: float = BOND_CA, *, tol: float = 1e-7,
                 max_sweeps: int = 200) -> np.ndarray:
    """Project a ``(..., n, 3)`` Cα trace onto ``|ca[i+1] - ca[i]| == bond``.

    A sampled or interpolated trace is not a chain: consecutive Cα atoms are a fixed virtual bond
    apart and nothing in a Gaussian draw or a bilinear read-back knows that. One SHAKE-style sweep
    corrects every even-index pair at once and then every odd-index pair --- within each parity the
    pairs are disjoint, so each half-sweep is an exact simultaneous projection, and alternating them
    converges.

    **The constant is measured, and its exception is not modelled.** Over 642 CDR1/CDR2/CDR3 Cα gaps
    in the 17 bundled crystals the separation is **3.7985 Å, s.d. 0.0668**, with 638 of 642
    (99.38 %) between 3.5 and 3.9 Å. The three that are not are **cis** peptide bonds at 2.854,
    2.892 and 2.996 Å (0.47 %), which this projection deliberately drives to ``bond`` like any
    other: the point is to make a *generated* trace physical, not to reproduce one crystal's cis
    proline. A caller reconstructing a specific solved loop should not use this.

    Args:
        ca: ``(..., n, 3)`` Cα coordinates; leading axes are independent traces.
        bond: Target Cα–Cα separation, Å.
        tol: Stop once **every** gap of every trace in the stack is within this of ``bond``.
            Sweeps needed vary with how far the input starts out, so this iterates to a residual
            rather than a fixed count. The default is 10 million times finer than the 0.5 Å bin the
            model's distance term reads, and leaves headroom for float accumulation over a big
            stack --- 1e-9 is not reliably reachable over a few hundred traces at once.
        max_sweeps: Give up after this many even/odd sweeps.

    Returns:
        A new array of the same shape whose consecutive gaps are ``bond`` to within ``tol``.

    Raises:
        ValueError: if the last axis is not length 3, if two consecutive Cα coincide exactly (the
            correction direction is undefined), or if ``max_sweeps`` is reached without converging.

    Example:
        >>> rng = np.random.default_rng(0)
        >>> trace = enforce_bond(rng.normal(0.0, 4.0, (12, 3)))
        >>> bool(np.allclose(np.linalg.norm(np.diff(trace, axis=0), axis=1), BOND_CA))
        True
    """
    out = np.array(ca, dtype=float, copy=True)
    if out.shape[-1] != 3:
        raise ValueError(f"expected (..., n, 3) coordinates, got shape {out.shape}")
    n = out.shape[-2]
    if n < 2:
        return out
    for _ in range(max_sweeps):
        for start in (0, 1):
            i = np.arange(start, n - 1, 2)
            if not len(i):
                continue
            d = out[..., i + 1, :] - out[..., i, :]
            length = np.linalg.norm(d, axis=-1)
            if np.any(length == 0.0):
                raise ValueError("coincident consecutive Cα atoms; bond direction is undefined")
            shift = (0.5 * (1.0 - bond / length))[..., None] * d
            out[..., i, :] += shift
            out[..., i + 1, :] -= shift
        gaps = np.linalg.norm(np.diff(out, axis=-2), axis=-1)
        if np.abs(gaps - bond).max() <= tol:
            return out
    raise ValueError(f"enforce_bond did not converge in {max_sweeps} sweeps; worst gap is off by "
                     f"{np.abs(gaps - bond).max():.3e} Å")



def reach_max(length: int, neck: float, bond: float = BOND_CA) -> float:
    """Maximum apex reach ``sqrt(((L+1)b/2)^2 - (R/2)^2)`` of a length-``length`` loop pinned at ``neck``.

    Args:
        length: CDR3 loop length ``L`` (number of residues).
        neck: neck separation ``R`` in Å (the ``ext`` frame descriptor).
        bond: virtual Cα–Cα bond length (default :data:`BOND_CA`).

    Returns:
        Maximum reach in Å, or ``0.0`` if the loop is too short to span the neck.
    """
    half = (length + 1) * bond / 2.0
    return math.sqrt(half * half - (neck / 2.0) ** 2) if half > neck / 2.0 else 0.0


def reachability_floor(distance: float, neck: float, bond: float = BOND_CA) -> int:
    """Shortest loop length that can reach a target at ``distance`` from a neck of separation ``neck``.

    Inverts :func:`reach_max`: the minimum ``L`` with ``reach_max(L, neck) >= distance``.
    """
    return int(math.ceil((2.0 / bond) * math.sqrt(distance * distance + (neck / 2.0) ** 2) - 1.0))


def span_saturation(structure) -> dict[str, float]:
    """Per-loop ``reach / reach_max`` for ``cdr3a`` and ``cdr3b`` of a chain-typed structure.

    Composes the existing CDR3 frame descriptors: ``reach`` (loop-centroid distance from the groove
    origin) and ``ext`` (neck separation ``R``), with ``L`` the resolved CDR3 length. ``> 1`` flags an
    infeasible pose (the loop cannot reach that far); real structures sit ~0.3–0.5. Returns ``NaN`` for
    a loop whose frame or CDR3 span is undefined.

    Note:
        The structure must be chain-typed (``classify_chains``) so its CDR3 regions are populated.
    """
    from .recognition import _cdr3_frame_features

    frame = _cdr3_frame_features(structure)
    lengths = _cdr3_lengths(structure)
    out: dict[str, float] = {}
    for loop in ("cdr3a", "cdr3b"):
        reach, neck, length = frame.get(f"{loop}_reach"), frame.get(f"{loop}_ext"), lengths.get(loop)
        rmax = reach_max(length, neck) if (length and neck is not None and math.isfinite(neck)) else 0.0
        out[loop] = reach / rmax if (rmax > 0 and reach is not None and math.isfinite(reach)) else math.nan
    return out


def _cdr3_lengths(structure) -> dict[str, int]:
    """Resolved CDR3 Cα count per loop (``cdr3a``/``cdr3b``) from the region markup."""
    out: dict[str, int] = {}
    for loop, ctype in (("cdr3a", "TRA"), ("cdr3b", "TRB")):
        for c in structure.chains:
            if getattr(c, "chain_type", None) != ctype:
                continue
            for reg in getattr(c, "regions", []) or []:
                if reg.region_type == "CDR3":
                    out[loop] = sum(1 for r in reg.residues if r.ca is not None)
                    break
    return out


@dataclass(slots=True)
class LoopInternalCoords:
    """Virtual-bond internal coordinates of a CDR3 Cα trace (the input to the Julia loop model).

    A Cα chain of ``N`` points has ``N-1`` bonds, ``N-2`` bond angles (interior vertices) and ``N-3``
    pseudo-torsions — the inverse of a discrete-Frenet reconstruction. ``neck`` is the end-to-end
    span ``|p_last - p_first|`` (the ``ext`` frame descriptor). Angles and torsions are in **radians**.
    """

    bonds: object            #: (N-1,) virtual Cα–Cα bond lengths, Å
    angles: object           #: (N-2,) Cα–Cα–Cα pseudo bond angles θ, radians
    torsions: object         #: (N-3,) Cα pseudo-dihedrals τ, radians in (−π, π]
    neck: float              #: end-to-end span ``R = |p_last − p_first|``, Å
    n_ca: int


def _bond_angle(p1, p2, p3) -> float:
    """Cα–Cα–Cα pseudo bond angle at ``p2`` (radians)."""
    a = p1 - p2
    b = p3 - p2
    c = float(a.dot(b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))
    return math.acos(max(-1.0, min(1.0, c)))


def cdr3_internal_coords(structure, loop: str = "cdr3b") -> LoopInternalCoords | None:
    """Virtual-bond internal coordinates (θ, τ, neck R) of a CDR3 loop's Cα trace.

    A thin ``structure -> coordinates`` extractor: the loop's *model* (torsion sampling, closure,
    ensemble entropy) lives in the algem monograph's Julia code and does not migrate — tcren only
    hands off the intrinsic internal coordinates of the observed loop. Reuses the CDR3 region markup
    (``classify_chains``) for the span and the shared ``_dihedral`` primitive for torsions.

    Args:
        structure: a chain-typed structure.
        loop: ``"cdr3a"`` (TRA) or ``"cdr3b"`` (TRB).

    Returns:
        A :class:`LoopInternalCoords`, or ``None`` if the loop has fewer than 4 resolved Cα (a
        torsion needs 4 points).
    """
    from .contactmap import _cdr3_ca
    from .docking.tcrdock_geometry import _dihedral

    ctype = {"cdr3a": "TRA", "cdr3b": "TRB"}.get(loop)
    if ctype is None:
        raise ValueError(f"loop must be 'cdr3a' or 'cdr3b', got {loop!r}")
    ca = _cdr3_ca(structure, ctype)
    if ca is None or len(ca) < 4:
        return None
    bonds = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    angles = np.array([_bond_angle(ca[k], ca[k + 1], ca[k + 2]) for k in range(len(ca) - 2)])
    torsions = np.array([_dihedral(ca[k], ca[k + 1], ca[k + 2], ca[k + 3])
                         for k in range(len(ca) - 3)])
    neck = float(np.linalg.norm(ca[-1] - ca[0]))
    return LoopInternalCoords(bonds=bonds, angles=angles, torsions=torsions, neck=neck, n_ca=len(ca))
