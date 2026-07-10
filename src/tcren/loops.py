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

import re
from dataclasses import dataclass

import numpy as np

__all__ = [
    "JUNCTION_MOTIF", "NECK_RANGE", "Junction",
    "find_junctions", "omega_stats", "is_omega_loop",
    "frenet", "kabsch_rmsd", "block_layouts", "structural_block_position",
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
    C-alpha RMSD; returns ``(best_position, best_rmsd, rmsd_per_layout)``. This is the
    structural ground truth against which a sequence-only block choice is scored.
    """
    rmsds = [kabsch_rmsd(np.asarray([ca_q[x] for x, _ in pairs]),
                         np.asarray([ca_r[y] for _, y in pairs]))
             for pairs in block_layouts(len(ca_q), len(ca_r))]
    best = int(np.argmin(rmsds))
    return best, rmsds[best], rmsds
