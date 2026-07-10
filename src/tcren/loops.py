"""Omega-loop geometry for immunoglobulin/TCR junctions.

A TCR (or antibody) junction is the segment from the conserved Cys104 to the Phe/Trp118 of
the J-region ``[FW]GXG`` motif -- what AIRR calls ``junction``, anchors included. Measured
over the 374 canonical TCR-pMHC complexes (731 junctions, 360 TRA + 371 TRB), it is a
textbook omega loop by the original Leszczynski & Rose criterion as quoted in
Fetrow 1995 (FASEB J 9:708-717, p.709):

    "the segment size was restricted to be between 6 and 16 residues, and the distance
    between segment termini was strictly limited to less than 10 A and less than two-thirds
    the longest alpha-carbon-alpha-carbon distance across the segment"

Three conditions, not two -- the third (compactness) is what separates a loop from an
extended segment, and it is the one that matters here:

============ ==== =========== ============= ========= ========== ======== ===========
region         n   median len  d_end (A)     d/d_max   (a) 6-16   (b) <10  all three
============ ==== =========== ============= ========= ========== ======== ===========
CDR1          277      5          10.97       1.00        44%       2%        0%
CDR2          271      6           6.35       0.67        90%      99%       49%
CDR3 (IMGT)   378     11           5.32       0.37       100%     100%       99%
junction      378     13           5.82       0.32        98%     100%       98%
============ ==== =========== ============= ========= ========== ======== ===========

**CDR1 is not an omega loop**: its termini are the most distant pair in the segment
(``d_end/d_max == 1.00``), so the chain never reverses direction. CDR2 sits exactly on the
two-thirds boundary. Only CDR3 and the junction qualify, and they qualify overwhelmingly.
Everything in this module is therefore scoped to the junction.

The neck is nearly invariant -- Ca(Cys)-Ca(Phe) is 5.81 A with a standard deviation of
0.48 A -- while the junction spans 8 to 19 residues. A fixed-geometry neck holding a
variable-length body is what makes a gap-position prior learnable and transferable across
loops. (The same "bottleneck" geometry, with a rigid N-terminal half and a flexible
C-terminal one, is described for the 16-residue omega loop of TEM-type beta-lactamases:
Egorov, Ulyashova & Rubtsova, Biomolecules 2019;9(12):854, doi:10.3390/biom9120854.)

Note: no atom pair sits at 4.5 A. Over 731 junctions the minimum Ca(C)-Ca(F/W) distance is
4.98 A and only 3.3% fall in 4.0-5.5 A, so a 4.5 A cutoff would select essentially nothing.
:data:`NECK_RANGE` records the measured band instead.
"""
from __future__ import annotations

import itertools
import re
from dataclasses import dataclass

import numpy as np

__all__ = [
    "JUNCTION_MOTIF", "NECK_RANGE", "Junction",
    "find_junctions", "omega_stats", "is_omega_loop",
    "frenet", "kabsch_rmsd", "block_layouts", "structural_block_position",
    "structural_align", "gap_runs", "is_single_block",
    "frenet_frame", "virtual_cb", "cb_orientation", "ramachandran",
]

#: The J-region motif that closes the junction: Phe or Trp, then Gly-X-Gly.
JUNCTION_MOTIF = re.compile(r"[FW]G.G")

#: Measured Ca(Cys) - Ca(Phe/Trp) neck distance band, in Angstrom (5th-95th percentile over
#: 731 junctions from data/Canonical2026). Use this, not the folklore 4.5 A.
NECK_RANGE = (5.0, 7.5)

# Original Leszczynski & Rose omega-loop criterion, as quoted by Fetrow 1995 p.709.
_MIN_LEN, _MAX_LEN = 6, 16
_MAX_END_DIST = 10.0
_MAX_END_RATIO = 2.0 / 3.0

_MIN_JUNCTION, _MAX_JUNCTION = 8, 26


@dataclass(frozen=True)
class Junction:
    """One ``C ... [FW]GXG`` junction located in a chain.

    Attributes:
        cys: Index of the conserved Cys in the chain's residue list.
        fw: Index of the Phe/Trp that closes the loop.
        seq: Junction sequence, anchors included (AIRR ``junction``).
        ca: ``(len(seq), 3)`` array of C-alpha coordinates, in order.
    """
    cys: int
    fw: int
    seq: str
    ca: np.ndarray

    @property
    def cdr3(self) -> str:
        """IMGT CDR3: the junction with both anchors stripped. Two residues shorter."""
        return self.seq[1:-1]

    @property
    def neck(self) -> float:
        """Ca(Cys) - Ca(Phe/Trp) distance in Angstrom."""
        return float(np.linalg.norm(self.ca[0] - self.ca[-1]))


def find_junctions(seq: str, ca: np.ndarray) -> list[Junction]:
    """Locate every ``C ... [FW]GXG`` junction in a chain.

    For each ``[FW]GXG`` occurrence, walks back to the nearest preceding Cys that yields a
    junction of plausible length (8-26 residues). Validated against
    ``notebooks/natcompsci2022/results_new/markup_2026.csv``: the recovered CDR3 matches the
    curated column exactly for the structures tested.

    Args:
        seq: One-letter chain sequence over resolved residues.
        ca: ``(len(seq), 3)`` C-alpha coordinates, aligned with ``seq``.

    Returns:
        Junctions in N-to-C order.
    """
    if len(seq) != len(ca):
        raise ValueError(f"sequence has {len(seq)} residues but {len(ca)} CA coordinates")
    out: list[Junction] = []
    for mo in JUNCTION_MOTIF.finditer(seq):
        f = mo.start()
        for c in range(f - 1, max(-1, f - _MAX_JUNCTION - 1), -1):
            if seq[c] != "C":
                continue
            if _MIN_JUNCTION <= f - c + 1 <= _MAX_JUNCTION:
                out.append(Junction(c, f, seq[c:f + 1], np.asarray(ca[c:f + 1], dtype=float)))
            break  # nearest preceding Cys only
    return out


def omega_stats(ca: np.ndarray) -> dict:
    """Geometry behind the three omega-loop conditions.

    Returns ``n_residues``, ``d_end`` (terminus separation), ``d_max`` (longest Ca-Ca
    distance across the segment), ``ratio`` = ``d_end / d_max``, and a bool per condition.
    """
    ca = np.asarray(ca, dtype=float)
    if len(ca) < 3:
        raise ValueError("need at least 3 C-alpha coordinates")
    d_end = float(np.linalg.norm(ca[0] - ca[-1]))
    d_max = float(np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=-1).max())
    return {
        "n_residues": len(ca),
        "d_end": d_end,
        "d_max": d_max,
        "ratio": d_end / d_max if d_max else float("inf"),
        "length_ok": _MIN_LEN <= len(ca) <= _MAX_LEN,
        "termini_close": d_end < _MAX_END_DIST,
        "compact": d_end < _MAX_END_RATIO * d_max,
    }


def is_omega_loop(ca: np.ndarray, relax_length: bool = False) -> bool:
    """All three Leszczynski-Rose conditions.

    ``relax_length=True`` drops the 16-residue ceiling, following Fetrow's note that the
    working definition was relaxed for longer segments. Needed for antibody CDR-H3, which
    reaches 30+ residues.
    """
    s = omega_stats(ca)
    length_ok = s["n_residues"] >= _MIN_LEN if relax_length else s["length_ok"]
    return bool(length_ok and s["termini_close"] and s["compact"])


def frenet(ca: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Discrete Frenet curvature and torsion along a C-alpha trace, in degrees.

    Following Hu, Lundgren & Niemi, *Phys. Rev. E* **83**, 061908 (2011),
    doi:10.1103/PhysRevE.83.061908 (arXiv:1102.5658), eqs. 14-22. With
    ``t_i = (r_{i+1} - r_i) / |.|`` the unit tangents and
    ``b_i = (t_{i-1} x t_i) / |.|`` the binormals:

    * ``cos kappa_i = t_i . t_{i+1}``  -- the virtual bond angle
    * ``cos tau_i   = b_i . b_{i+1}``  -- the virtual torsion angle, signed

    Returns ``(kappa, tau)`` of lengths ``n-2`` and ``n-3``. Both are invariant under rigid
    motion, so they describe loop *shape* without superposition -- the basis for a structural
    alphabet. Torsion is undefined at an inflection point (collinear tangents); those entries
    come back as ``nan``.
    """
    ca = np.asarray(ca, dtype=float)
    if len(ca) < 4:
        raise ValueError("need at least 4 C-alpha coordinates for curvature and torsion")
    t = np.diff(ca, axis=0)
    norms = np.linalg.norm(t, axis=1, keepdims=True)
    t = t / norms

    kappa = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", t[:-1], t[1:]), -1.0, 1.0)))

    cross = np.cross(t[:-1], t[1:])
    cn = np.linalg.norm(cross, axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        b = np.where(cn > 1e-8, cross / cn, np.nan)

    cos_tau = np.clip(np.einsum("ij,ij->i", b[:-1], b[1:]), -1.0, 1.0)
    sin_tau = np.einsum("ij,ij->i", np.cross(b[:-1], b[1:]), t[1:-1])
    tau = np.degrees(np.arctan2(sin_tau, cos_tau))
    return kappa, tau


def frenet_frame(ca: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Orthonormal ``(t, n, b)`` frame attached to each interior C-alpha.

    Built from the same tangents and binormals as :func:`frenet`: ``t_i`` is the unit bond
    ``i -> i+1``, ``b_i`` the unit binormal ``t_{i-1} x t_i``, and ``n_i = b_i x t_i``. Defined
    for residues ``1 .. len(ca) - 2``, so each array has ``len(ca) - 2`` rows.

    The frame is what lets a side-chain direction be expressed without superposition: see
    :func:`cb_orientation`. Rows where consecutive tangents are collinear come back as ``nan``.
    """
    ca = np.asarray(ca, dtype=float)
    if len(ca) < 4:
        raise ValueError("need at least 4 C-alpha coordinates for a Frenet frame")
    t = np.diff(ca, axis=0)
    t = t / np.linalg.norm(t, axis=1, keepdims=True)

    cross = np.cross(t[:-1], t[1:])
    cn = np.linalg.norm(cross, axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        b = np.where(cn > 1e-8, cross / cn, np.nan)
    tt = t[1:]                      # tangent at residue i, for i = 1 .. len(ca)-2
    n = np.cross(b, tt)
    return tt, n, b


def virtual_cb(nn: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Idealised C-beta position from the backbone N, CA and C of the same residue.

    Glycine has no C-beta, so its side-chain direction has to be constructed. Uses the standard
    tetrahedral placement; ``scripts/cb_contacts.py`` checks it against the observed C-beta of
    every non-Gly residue and refuses to proceed if the RMSD exceeds 0.15 A.
    """
    nn, ca, c = (np.asarray(x, dtype=float) for x in (nn, ca, c))
    b = ca - nn
    cc = c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def cb_orientation(ca: np.ndarray, cb: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Side-chain direction in the local Frenet frame, in degrees.

    For each interior residue, ``u = unit(cb - ca)`` is resolved against ``(t, n, b)``:

    * ``polar``   -- angle between ``u`` and the tangent, in ``[0, 180]``
    * ``azimuth`` -- ``atan2(u . b, u . n)`` in ``(-180, 180]``, the rotation about the tangent

    Both are invariant under rigid motion, so they describe *where the side chain points*
    relative to the loop's own geometry without any superposition. The azimuth flips sign under
    reflection, exactly as the Frenet torsion does, so it carries the chain's handedness.

    .. note::
       It does **not** help predict which junction residues touch the peptide. Over 3,883
       interior residues in 368 crystal junctions, the azimuth alone reaches ROC-AUC 0.728, but
       simple distance from the loop apex reaches 0.847, and adding the C-beta orientation on
       top of position gains **+0.000** (5-fold logistic, grouped by junction). Measured by
       ``scripts/cb_contacts.py``. Do not build a positional weight profile out of it.

    Args:
        ca: ``(n, 3)`` C-alpha coordinates.
        cb: ``(n, 3)`` C-beta coordinates, real or from :func:`virtual_cb`.

    Returns:
        ``(polar, azimuth)``, each of length ``n - 2``.
    """
    ca = np.asarray(ca, dtype=float)
    cb = np.asarray(cb, dtype=float)
    if ca.shape != cb.shape:
        raise ValueError(f"shape mismatch: {ca.shape} vs {cb.shape}")
    t, n, b = frenet_frame(ca)
    u = cb[1:-1] - ca[1:-1]
    u = u / np.linalg.norm(u, axis=1, keepdims=True)
    polar = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", u, t), -1.0, 1.0)))
    azimuth = np.degrees(np.arctan2(np.einsum("ij,ij->i", u, b), np.einsum("ij,ij->i", u, n)))
    return polar, azimuth


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1).sum(-1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(-1, keepdims=True) * b1
    x = (v * w).sum(-1)
    y = (np.cross(b1, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def ramachandran(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Backbone ``(phi, psi)`` in degrees, for residues ``1 .. len(ca) - 2``.

    Needs N and C, which :func:`frenet` does not. Whether it *adds* anything over the C-alpha
    curvature and torsion was measured, not assumed (``scripts/cb_contacts.py``, 3,883 residues,
    5-fold k-NN): ``psi`` is recoverable from ``(kappa, tau)`` with circular R^2 = 0.822 and is
    redundant; ``phi`` is not (R^2 = 0.442) and does carry independent information.

    That information is nonetheless worth nothing downstream: adding ``phi`` on top of loop
    position and shape gains **+0.000** ROC-AUC for peptide-contact prediction. Provided for
    completeness; nothing in this package consumes it.
    """
    n, ca, c = (np.asarray(x, dtype=float) for x in (n, ca, c))
    if not (len(n) == len(ca) == len(c)) or len(ca) < 3:
        raise ValueError("need matching N, CA, C arrays of at least 3 residues")
    i = np.arange(1, len(ca) - 1)
    phi = _dihedral(c[i - 1], n[i], ca[i], c[i])
    psi = _dihedral(n[i], ca[i], c[i], n[i + 1])
    return phi, psi


def kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Minimal RMSD between two equal-length point sets after optimal superposition."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    if len(a) == 0:
        raise ValueError("empty point set")
    a = a - a.mean(0)
    b = b - b.mean(0)
    v, _, wt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    return float(np.sqrt((((a @ rot) - b) ** 2).sum(1).mean()))


def block_layouts(m: int, n: int) -> list[list[tuple[int, int]]]:
    """Residue correspondences induced by every single-gap-block position.

    Layout ``i`` pairs the first ``i`` residues on the diagonal, then jumps the block of
    ``d = |m - n|`` residues in the longer sequence. Mirrors ``seqtree.gapblock``: ``i``
    ranges over ``[0, min(m, n)]`` inclusive at both ends.
    """
    d = abs(m - n)
    shorter = min(m, n)
    out = []
    for i in range(shorter + 1):
        if m >= n:
            out.append([(j, j) for j in range(i)] + [(j + d, j) for j in range(i, n)])
        else:
            out.append([(j, j) for j in range(i)] + [(j, j + d) for j in range(i, m)])
    return out


def structural_block_position(ca_q: np.ndarray, ca_r: np.ndarray) -> tuple[int, float, list[float]]:
    """Gap-block position best supported by the two loops' backbones.

    For every layout, superposes the residue pairs it induces and reports the resulting
    C-alpha RMSD; returns ``(best_position, best_rmsd, rmsd_per_layout)``.

    .. warning::
       This ranks positions *within* the single-gap-block family, because
       :func:`block_layouts` enumerates only that family. It therefore cannot test whether one
       block is enough -- it has no way to express any other answer. Use
       :func:`structural_align` for a correspondence that makes no such assumption.
    """
    rmsds = [kabsch_rmsd(np.asarray([ca_q[x] for x, _ in pairs]),
                         np.asarray([ca_r[y] for _, y in pairs]))
             for pairs in block_layouts(len(ca_q), len(ca_r))]
    best = int(np.argmin(rmsds))
    return best, rmsds[best], rmsds


# ---------------------------------------------------------------- model-independent alignment

def _kabsch_transform(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation and translation carrying ``a`` onto ``b``: ``a @ rot + shift``."""
    ca, cb = a.mean(0), b.mean(0)
    v, _, wt = np.linalg.svd((a - ca).T @ (b - cb))
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    return rot, cb - ca @ rot


def _affine_dp(score: np.ndarray, gap_open: float, gap_extend: float) -> str:
    """Unrestricted global affine alignment maximising ``score``. Returns an op string.

    Ops are ``M`` (a pair), ``D`` (a residue of the row sequence against a gap) and ``I`` (a
    residue of the column sequence against a gap). Any number of gap blocks, anywhere -- the
    whole point: nothing here presumes a single block.
    """
    m, n = score.shape
    neg = -np.inf
    M = np.full((m + 1, n + 1), neg)
    X = np.full((m + 1, n + 1), neg)   # gap in the column sequence, row residue consumed
    Y = np.full((m + 1, n + 1), neg)   # gap in the row sequence, column residue consumed
    M[0, 0] = 0.0
    for i in range(1, m + 1):
        X[i, 0] = -gap_open - (i - 1) * gap_extend
    for j in range(1, n + 1):
        Y[0, j] = -gap_open - (j - 1) * gap_extend
    bM = np.zeros((m + 1, n + 1), np.int8)
    bX = np.zeros((m + 1, n + 1), np.int8)
    bY = np.zeros((m + 1, n + 1), np.int8)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cand = (M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1])
            k = int(np.argmax(cand))
            M[i, j] = score[i - 1, j - 1] + cand[k]
            bM[i, j] = k

            cand = (M[i - 1, j] - gap_open, X[i - 1, j] - gap_extend, Y[i - 1, j] - gap_open)
            k = int(np.argmax(cand))
            X[i, j] = cand[k]
            bX[i, j] = k

            cand = (M[i, j - 1] - gap_open, X[i, j - 1] - gap_open, Y[i, j - 1] - gap_extend)
            k = int(np.argmax(cand))
            Y[i, j] = cand[k]
            bY[i, j] = k

    i, j = m, n
    state = int(np.argmax((M[m, n], X[m, n], Y[m, n])))
    ops: list[str] = []
    while i > 0 or j > 0:
        if state == 0:
            ops.append("M")
            state = int(bM[i, j])
            i, j = i - 1, j - 1
        elif state == 1:
            ops.append("D")
            state = int(bX[i, j])
            i -= 1
        else:
            ops.append("I")
            state = int(bY[i, j])
            j -= 1
    return "".join(reversed(ops))


def _pairs_from_ops(ops: str) -> list[tuple[int, int]]:
    pairs, i, j = [], 0, 0
    for op in ops:
        if op == "M":
            pairs.append((i, j))
            i, j = i + 1, j + 1
        elif op == "D":
            i += 1
        else:
            j += 1
    return pairs


def structural_align(
    ca_q: np.ndarray,
    ca_r: np.ndarray,
    d0: float = 3.0,
    gap_open: float = 0.6,
    gap_extend: float = 0.1,
    max_iter: int = 20,
) -> tuple[list[tuple[int, int]], float, str]:
    """Residue correspondence between two C-alpha traces, with **no** gap-model assumption.

    Seeds the superposition on the loop's own anchors -- the three N-terminal and three
    C-terminal residues, whose Ca(Cys)-Ca(Phe) separation is invariant to +/- 0.48 A across
    junctions -- then iterates *superpose, rescore, realign* to a fixed point. The realignment
    is an unrestricted affine DP, free to open any number of gap blocks anywhere.

    This is the oracle :func:`structural_block_position` cannot be: it can return a
    correspondence that is *not* a single block, so asking how often it does is a real
    question with a real answer.

    Args:
        ca_q: ``(m, 3)`` C-alpha coordinates.
        ca_r: ``(n, 3)`` C-alpha coordinates.
        d0: Distance scale of the residue similarity ``1 / (1 + (dist/d0)^2)``, in Angstrom.
        gap_open: Similarity forfeited to open a gap block.
        gap_extend: Similarity forfeited per additional gap column.
        max_iter: Fixed-point iteration cap.

    Returns:
        ``(pairs, rmsd, ops)`` -- the matched residue index pairs, their superposed C-alpha
        RMSD, and the alignment op string over ``M``/``D``/``I``.

    Raises:
        ValueError: If either trace has fewer than 6 residues (the anchor seed needs them).
    """
    q = np.asarray(ca_q, dtype=float)
    r = np.asarray(ca_r, dtype=float)
    if len(q) < 6 or len(r) < 6:
        raise ValueError("need at least 6 C-alpha coordinates per loop to seed on the anchors")

    pairs = [(0, 0), (1, 1), (2, 2),
             (len(q) - 3, len(r) - 3), (len(q) - 2, len(r) - 2), (len(q) - 1, len(r) - 1)]
    ops = ""
    for _ in range(max_iter):
        rot, shift = _kabsch_transform(np.array([q[i] for i, _ in pairs]),
                                       np.array([r[j] for _, j in pairs]))
        moved = q @ rot + shift
        dist = np.linalg.norm(moved[:, None, :] - r[None, :, :], axis=-1)
        ops_new = _affine_dp(1.0 / (1.0 + (dist / d0) ** 2), gap_open, gap_extend)
        pairs_new = _pairs_from_ops(ops_new)
        if len(pairs_new) < 3:
            break
        if pairs_new == pairs:
            ops = ops_new
            break
        pairs, ops = pairs_new, ops_new

    rot, shift = _kabsch_transform(np.array([q[i] for i, _ in pairs]),
                                   np.array([r[j] for _, j in pairs]))
    moved = np.array([q[i] for i, _ in pairs]) @ rot + shift
    ref = np.array([r[j] for _, j in pairs])
    rmsd = float(np.sqrt(((moved - ref) ** 2).sum(1).mean()))
    return pairs, rmsd, ops


def gap_runs(ops: str) -> list[tuple[str, int, int]]:
    """Maximal runs of gap ops as ``(op, start_column, length)``."""
    out, col = [], 0
    for op, grp in itertools.groupby(ops):
        n = len(list(grp))
        if op != "M":
            out.append((op, col, n))
        col += n
    return out


def is_single_block(ops: str) -> bool:
    """Does this alignment place all its gaps in one contiguous block of one sequence?

    ``True`` for ``MMMDDMMM`` and for a gapless ``MMMM``; ``False`` for ``MMDMMIMM``, which
    needs two blocks. This is the predicate that decides whether the single-gap-block model
    is a restriction the data actually pays for.
    """
    return len(gap_runs(ops)) <= 1
