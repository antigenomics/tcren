"""TCR:pMHC docking geometry — native reimplementation of TCRdock's rigid-body parameterisation.

Reimplemented from the Bradley lab's TCRdock (https://github.com/phbradley/TCRdock, MIT license), commit
``c5a7af42eeb0c2a4492a4d4fe803f1f9aafb6193`` (2024-03-04), specifically ``tcrdock/docking_geometry.py``,
``tcrdock/mhc_util.py``, ``tcrdock/tcr_util.py`` and ``tcrdock/superimpose.py``. No TCRdock code is imported;
the geometry is ported to tcren's own :class:`~tcren.structure.model.Structure` and annotation.

The docking geometry describes how the TCR sits on the peptide-MHC groove as a rigid-body transform between
two coordinate frames ("stubs"):

* the **MHC stub** — from the ~180° pseudo-symmetry of the class-I α1α2 (or class-II α1β1) β-sheet floor:
  x points toward the peptide, z from one half of the sheet to the other, origin at the two halves' midpoint;
* the **TCR stub** — from the ~180° pseudo-symmetry relating Vα and Vβ: x toward the CDR loops, z from Vα to
  Vβ, origin at the two domains' midpoint.

Six numbers (:class:`DockingGeometry`) fix the relative pose: ``d`` (frame separation), ``torsion`` (dihedral
about the MHC–TCR line), and ``tcr_unit_y/z`` + ``mhc_unit_y/z`` (each frame's direction to the other, in the
other's local axes). For interpretable in-plane / tilt angles use :func:`tcren.orient.docking.docking_angles`
(``crossing_angle`` = the groove-plane "scanning" angle; ``incident_angle`` = the tilt); this module adds the
full rigid-body pose that those two scalars do not capture.

MHC-I core positions are mapped by BLOSUM-aligning the α chain to TCRdock's class-I template; TCR core
positions are the conserved IMGT framework positions, located from tcren's arda region markup.

Provenance note (validated on 618 TCRvdb TCRmodel2 models, 2026-07-05): tcren's ``crossing_angle`` reproduces
the "scanning_angle" reported by upstream AF/TCRmodel2 annotation tables (r≈0.88), so that quantity is a
genuine, reproducible interface geometry. The upstream "pitch_angle", however, matches **no** clean geometric
angle (best correlate is ``d``, r≈0.42) and discriminates TCRvdb binders better (macro-PR≈0.72) than any
clean docking feature computed here (d≈0.64, torsion≈0.62, tilt≈0.58) — i.e. its extra signal is not
reproducible from coordinates and is likely AlphaFold-confidence contamination. Prefer these documented,
crystal-computable descriptors over the opaque upstream pitch.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

_TCR_TYPES = ("TRA", "TRB", "TRD", "TRG")
_AB_PAIR = ("TRA", "TRB")
_GD_PAIR = ("TRD", "TRG")
_MHC_A_TYPES = ("MHCa", "MHC")   # class-I heavy chain / class-II α (groove-bearing)
_MHC_B_TYPES = ("MHCb",)         # class-II β

# --- TCRdock class-I template (3pqyA) and its 12 β-sheet-floor core columns (1-indexed) ------------------
CLASS1_TEMPLATE_SEQ = (
    "PHSMRYFETAVSRPGLEEPRYISVGYVDNKEFVRFDSDAENPRYEPRAPWMEQEGPEYWERETQKAKGQEQWFRVSLRNLLGYYNQSAGGSHTLQQMSGC"
    "DLGSDWRLLRGYLQFAYEGRDYIALNEDLKTWTAADMAAQITRRKWEQSGAAEHYKAYLEGECVEWLHRYLKNGNATLLRTDSPKAHVTHHPRSKGEVTL"
    "RCWALGFYPADITLTWQLNGEELTQDMELVETRPAGDGTFQKWASVVVPLGKEQNYTCRVYHEGLPEPLTLRWEP"
)
CLASS1_CORE_POS_0X = [p - 1 for p in (4, 6, 8, 10, 23, 25, 94, 96, 98, 100, 113, 115)]

# --- TCR core: conserved IMGT framework positions located from region boundaries -------------------------
# TCRdock uses 13 IMGT positions [21,23,25, 39,41, 53,54,55, 78, 89, 102,103,104]; we use the 11 that anchor
# unambiguously off a region boundary (IMGT 78/89 sit inside FR3's variable gap zone and are dropped — the
# symmetry frame averages over the rest). Anchors: (region, offset) in tcren seq_index.
#   FR1 tail (21,23,25) counts back from CDR1 start; FR2 head (39,41) forward from FR2 start; FR2 tail
#   (53,54,55) back from CDR2 start; FR3 tail (102,103,104) back from CDR3 start (104 = conserved Cys).
_TCR_CORE_ANCHORS = [
    ("CDR1", -6), ("CDR1", -4), ("CDR1", -2),       # IMGT 21, 23, 25
    ("FR2", 0), ("FR2", 2),                          # IMGT 39, 41
    ("CDR2", -3), ("CDR2", -2), ("CDR2", -1),        # IMGT 53, 54, 55
    ("CDR3", -3), ("CDR3", -2), ("CDR3", -1),        # IMGT 102, 103, 104
]


# =======================================================================================================
# geometry helpers (ported from TCRdock geom_util / superimpose; pure numpy)
# =======================================================================================================
def _kabsch(fix: np.ndarray, mov: np.ndarray) -> np.ndarray:
    """Return ``mov`` rigidly superimposed onto ``fix`` (least-squares, Kabsch)."""
    fc, mc = fix.mean(0), mov.mean(0)
    H = (mov - mc).T @ (fix - fc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (R @ (mov - mc).T).T + fc


def _stub_from_three_points(a, b, c) -> tuple[np.ndarray, np.ndarray]:
    """Orthonormal axes (rows) from three points: x=a→ (a-a=0 → use a,b,c per TCRdock stub_from_four_points)."""
    ihat = (a - b)
    ihat = ihat / (np.linalg.norm(ihat) + 1e-12)
    jhat = (c - b) - ihat * ihat.dot(c - b)
    jhat = jhat / (np.linalg.norm(jhat) + 1e-12)
    khat = np.cross(ihat, jhat)
    return np.stack([ihat, jhat, khat]), a


def _rotation_axis(axes1: np.ndarray, axes2: np.ndarray) -> np.ndarray:
    """Unit axis of the rotation taking frame ``axes1`` to ``axes2`` (both rows), via matrix log."""
    R = axes2.T @ axes1
    # rotvec from R (Rodrigues): angle from trace, axis from skew part
    cos = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    ang = np.arccos(cos)
    if ang < 1e-8:
        return np.array([1.0, 0.0, 0.0])
    ax = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    n = np.linalg.norm(ax)
    if n < 1e-8:                       # ~180°: axis from the symmetric part (R + I)
        M = (R + np.eye(3)) / 2.0
        ax = M[:, int(np.argmax(np.diag(M)))]
        return ax / (np.linalg.norm(ax) + 1e-12)
    return ax / n


def _symmetry_stub(coords: np.ndarray, point_towards: np.ndarray) -> dict:
    """Stub from a set of core CAs whose two halves are related by ~180° (TCRdock get_symmetry_stub)."""
    m = len(coords) // 2
    acom, bcom = coords[:m].mean(0), coords[m:].mean(0)
    swapped = np.vstack([coords[m:], coords[:m]])
    axes1, _ = _stub_from_three_points(coords[0], coords[1], coords[2])
    moved = _kabsch(swapped, coords)                 # coords onto swapped
    axes2, _ = _stub_from_three_points(moved[0], moved[1], moved[2])
    x = _rotation_axis(axes1, axes2)                 # symmetry axis
    origin = 0.5 * (acom + bcom)
    if (point_towards - origin).dot(x) < 0:
        x = -x
    z = bcom - acom
    z = z - x * x.dot(z)
    z = z / (np.linalg.norm(z) + 1e-12)
    y = np.cross(z, x)
    return {"axes": np.stack([x, y, z]), "origin": origin}


def _global2local(stub: dict, v: np.ndarray) -> np.ndarray:
    return stub["axes"].dot(v - stub["origin"])


def _dihedral(p1, p2, p3, p4) -> float:
    a = p2 - p1; a /= np.linalg.norm(a)
    b = p3 - p2; b /= np.linalg.norm(b)
    c = p4 - p3; c /= np.linalg.norm(c)
    x = -a.dot(c) + a.dot(b) * b.dot(c)
    y = a.dot(np.cross(b, c))
    return float(np.arctan2(y, x))


# =======================================================================================================
# core-CA extraction from a chain-typed, annotated tcren Structure
# =======================================================================================================
def _region_start(chain, region_type: str) -> int | None:
    for r in chain.regions:
        if r.region_type == region_type:
            return r.start_seq_index
    return None


def _ca(chain, seq_index: int):
    r = chain.by_seq_index(seq_index)
    return None if r is None or r.ca is None else np.asarray(r.ca, float)


def _tcr_core_ca(structure, va_type: str, vb_type: str) -> np.ndarray | None:
    """13→11 conserved IMGT framework CAs for Vα then Vβ (corresponding order)."""
    out = []
    for ctype in (va_type, vb_type):
        chain = next((c for c in structure.chains if c.chain_type == ctype), None)
        if chain is None:
            return None
        starts = {reg: _region_start(chain, reg) for reg in ("CDR1", "FR2", "CDR2", "CDR3")}
        if any(v is None for v in starts.values()):
            return None
        for reg, off in _TCR_CORE_ANCHORS:
            ca = _ca(chain, starts[reg] + off)
            if ca is None:
                return None
            out.append(ca)
    return np.asarray(out)


def _blosum_map(template: str, query: str) -> dict[int, int]:
    """0-indexed template→query position map from a BLOSUM62 local alignment (TCRdock blosum_align)."""
    from Bio.Align import PairwiseAligner, substitution_matrices
    aln = PairwiseAligner()
    aln.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aln.mode = "local"
    aln.open_gap_score = -11
    aln.extend_gap_score = -1
    best = aln.align(template, query)[0]
    mp = {}
    for (t0, t1), (q0, q1) in zip(best.aligned[0], best.aligned[1]):
        for k in range(t1 - t0):
            mp[t0 + k] = q0 + k
    return mp


def _mhc_core_ca(structure) -> np.ndarray | None:
    """12 class-I β-sheet-floor CAs (α1 half then α2 half) via BLOSUM alignment to the TCRdock template."""
    chain = next((c for c in structure.chains if c.chain_type in _MHC_A_TYPES), None)
    if chain is None:
        return None
    seq = "".join(r.aa for r in chain.residues if getattr(r, "aa", None))
    if len(seq) < 120:
        return None
    tmap = _blosum_map(CLASS1_TEMPLATE_SEQ, seq)
    out = []
    for pos in CLASS1_CORE_POS_0X:
        q = tmap.get(pos)
        if q is None:
            return None
        ca = _ca(chain, q)
        if ca is None:
            return None
        out.append(ca)
    return np.asarray(out)


# =======================================================================================================
# public API
# =======================================================================================================
@dataclass(slots=True)
class DockingGeometry:
    """TCR:pMHC rigid-body docking geometry — TCRdock's 6-parameter form."""

    d: float                 # MHC-stub → TCR-stub origin separation (Å)
    torsion: float           # dihedral about the MHC–TCR line (radians, [0, 2π))
    tcr_unit_y: float        # TCR direction in MHC-local y
    tcr_unit_z: float        # TCR direction in MHC-local z
    mhc_unit_y: float        # MHC direction in TCR-local y
    mhc_unit_z: float        # MHC direction in TCR-local z

    def to_dict(self) -> dict:
        return asdict(self)


def docking_geometry(structure) -> DockingGeometry:
    """Compute the TCRdock docking geometry of a chain-typed, MHC-annotated TCR-pMHC structure.

    The structure must already be chain-typed (:func:`tcren.annotation.classify_chains`) and MHC-annotated
    (:func:`tcren.mhc.annotate_mhc`) with arda region markup on the TCR chains.

    Args:
        structure: a chain-typed, annotated TCR-pMHC :class:`~tcren.structure.model.Structure`.

    Returns:
        The :class:`DockingGeometry`.

    Raises:
        ValueError: if the MHC β-sheet core or a complete TCR Vα/Vβ core cannot be located.
    """
    mhc_ca = _mhc_core_ca(structure)
    if mhc_ca is None:
        raise ValueError("could not locate the class-I MHC β-sheet core (class-II not yet supported)")
    tcr_ca = None
    for va, vb in (_AB_PAIR, _GD_PAIR):
        tcr_ca = _tcr_core_ca(structure, va, vb)
        if tcr_ca is not None:
            break
    if tcr_ca is None:
        raise ValueError("could not locate a complete TCR Vα/Vβ framework core")

    # peptide + CDR-loop centroids disambiguate the two symmetry-axis directions
    pep = np.array([r.ca for c in structure.chains if c.chain_type == "PEPTIDE"
                    for r in c.residues if r.ca is not None], float)
    if len(pep) < 2:
        raise ValueError("peptide chain with ≥2 Cα required")
    cdr_cen = np.array([r.ca for c in structure.chains if c.chain_type in _TCR_TYPES
                        for reg in c.regions if reg.region_type.startswith("CDR")
                        for r in (c.by_seq_index(i) for i in range(reg.start_seq_index, reg.end_seq_index + 1))
                        if r is not None and r.ca is not None], float).mean(0)

    mhc_stub = _symmetry_stub(mhc_ca, point_towards=pep.mean(0))
    tcr_stub = _symmetry_stub(tcr_ca, point_towards=cdr_cen)

    torsion = _dihedral(mhc_stub["origin"] + mhc_stub["axes"][1], mhc_stub["origin"],
                        tcr_stub["origin"], tcr_stub["origin"] + tcr_stub["axes"][2])
    torsion = (torsion + 2 * np.pi) % (2 * np.pi)
    tcr_unit = _global2local(mhc_stub, tcr_stub["origin"]); tcr_unit /= np.linalg.norm(tcr_unit)
    mhc_unit = _global2local(tcr_stub, mhc_stub["origin"]); mhc_unit /= np.linalg.norm(mhc_unit)
    d = float(np.linalg.norm(mhc_stub["origin"] - tcr_stub["origin"]))
    return DockingGeometry(d, torsion, float(tcr_unit[1]), float(tcr_unit[2]),
                           float(mhc_unit[1]), float(mhc_unit[2]))
