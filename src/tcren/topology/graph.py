"""The interface as a graph, and as a matrix.

Two families, both read off the same complex, both free of the two free parameters the older
footprint measures carry.

**The graph.** :mod:`tcren.footprint` measures coverage by tallying residue-pair contacts into a
fixed partition -- twelve cells (six CDR loops x {peptide, MHC}) or twenty-four -- and then taking
the diversity of the tally. But the contact map at 5 A already *is* a bipartite graph, with the CDR
loop residues on one side and the pMHC residues they touch on the other, and the cell partition
throws that incidence structure away: it records how many contacts each bin holds, never which
residue touched which. Everything here is a functional of the biadjacency matrix ``B`` alone, so
there is no binning to choose. Degree evenness replaces cell entropy, the component count of the
**contact** graph replaces Betti-0 of a Calpha flag complex at an arbitrary 7 or 8 A, and the
normalised cyclomatic number replaces the raw one the footprint docstring rejects for tracking
interface size.

**The matrix.** A Calpha map between one region and another is an ``L x M`` matrix whose shape is
the two regions' lengths, so no entry of it can be compared across a 9-mer and a 15-mer, or across
CDR3 loops of different length. Its **singular values** can. Turn the distance map into a soft
adjacency ``K = exp(-D^2 / 2 sigma^2)`` and the normalised spectrum of ``K`` is a shape descriptor
that does not know how long either side was: the effective rank fraction says how many independent
approach modes the interface has, and ``s2/s1`` says how far it is from a rank-one (separable)
approach, which is what a receptor that leans on a surface rather than reading it produces. This is
**not** the graphon registration of :func:`tcren.contactmap.registered_map`, which resamples the map
onto a fixed grid and whose signal is epitope-identity provenance; a singular value is an invariant
of the map, not a resampling of it.

Alongside them, the one thing Calpha cannot see. ``d_Calpha - d_Cbeta`` over the contacting pairs is
positive when two side chains point at each other and negative when the backbones are close and the
side chains point away -- the shape a pose forced to satisfy a contact-count objective takes.

    >>> from tcren.interface_graph import graph_features, matrix_features
    >>> graph_features(structure)["g_loop_overlap"]      # doctest: +SKIP
    >>> matrix_features(structure)["m_erank_tp"]         # doctest: +SKIP
"""

from __future__ import annotations

import numpy as np
import polars as pl

from ..contactmap import ContactMap
from ..structure.model import MHC_TYPES, PEPTIDE_TYPE, Structure

__all__ = [
    "GRAPH_FEATURES",
    "MATRIX_FEATURES",
    "MHC_HELIX_REGIONS",
    "PROMOTED_POSE_FEATURES",
    "graph_features",
    "matrix_features",
]

#: The six CDR loops, as :data:`tcren.footprint.CELL_LOOPS` names them. Imported lazily in the
#: functions rather than at module scope so ``footprint`` may import this module without a cycle.
_CELL_LOOPS: tuple[str, ...] = ("TRA:CDR1", "TRA:CDR2", "TRA:CDR3",
                                "TRB:CDR1", "TRB:CDR2", "TRB:CDR3")

GRAPH_FEATURES: tuple[str, ...] = (
    "g_even_tcr", "g_even_pmhc", "g_comp_frac", "g_alg_conn",
    "g_cyclo_frac", "g_loop_even", "g_loop_overlap", "g_assort",
)

#: The Calpha-against-Cbeta block. ``m_face_tp``, ``ca_cb_agreement_tp`` and ``ca_cb_agreement_tm``
#: are **promotions**: :mod:`tcren.pose` has computed all three since the pose layer was written and
#: no catalogued family reached them. ``m_face_tp`` was called ``sidechain_toward`` there and is
#: renamed rather than duplicated -- catalguing a second name for a number the package already
#: computes is the defect the 2026-07-28 descriptor audit removed. Only ``m_face_tm`` is new.
MATRIX_FEATURES: tuple[str, ...] = (
    "m_erank_tp", "m_gap_tp", "m_erank_tm", "m_gap_tm",
    "m_face_tp", "m_face_tm", "ca_cb_agreement_tp", "ca_cb_agreement_tm",
)

#: The MHC groove regions the receptor meets. The floor is excluded: a receptor approaching from
#: above meets the two helix crests, and the floor residues it reaches at all it reaches through
#: the peptide.
MHC_HELIX_REGIONS: tuple[str, ...] = ("HELIX_A1", "HELIX_A2", "HELIX_B1")

#: Length scale of the Gaussian proximity kernel, Angstrom. The contact criterion's own 5 A, so the
#: binary contact map is this kernel's ``sigma -> 0`` limit rather than an independent choice.
_SIGMA = 5.0

#: Two descriptors that :mod:`tcren.pose` has computed since the pose layer was written and that no
#: catalogued family reached, because ``POSE_FEATURES`` is not ``DESCRIPTORS``. They are the
#: order-2 (participation ratio) reading of the same receptor-side degrees ``g_even_tcr`` reads at
#: order 1, over TCR:peptide alone rather than the whole pMHC. Computed here from the contact map
#: this module already builds, through :func:`tcren.pose._degree_descriptors` so there is exactly
#: one formula. ``max_degree_tp``, which that function also returns, is deliberately not promoted.
PROMOTED_POSE_FEATURES: tuple[str, ...] = ("degree_evenness_tp", "frac_well_coordinated_tp")


# --- the bipartite contact graph -----------------------------------------------------------------

def _biadjacency(structure: Structure, cutoff: float):
    """``(B, loops, n_rows, n_cols)`` -- the CDR-loop x pMHC residue biadjacency at ``cutoff``.

    ``B[i, j]`` is 1 when TCR residue ``i`` (of a CDR loop) has a heavy atom within ``cutoff`` of
    pMHC residue ``j``. Rows and columns are the residues that engage at all, in sorted key order so
    the matrix is byte-reproducible between runs; ``loops[i]`` is row ``i``'s loop name.

    Returns ``(None, None, 0, 0)`` when the structure makes no CDR-loop contact with the pMHC.
    """
    cm = ContactMap.from_structure(structure, cutoff=cutoff)
    frames = []
    for iface in ("tcr_peptide", "tcr_mhc"):
        d = cm.interface(iface)
        if d.is_empty():
            continue
        frames.append(d.select(
            pl.concat_str([pl.col("chain.type.from"), pl.col("region.type.from")],
                          separator=":").alias("loop"),
            pl.col("chain.id.from"), pl.col("residue.index.from"),
            pl.col("chain.id.to"), pl.col("residue.index.to")))
    if not frames:
        return None, None, 0, 0
    t = (pl.concat(frames, how="vertical")
         .filter(pl.col("loop").is_in(list(_CELL_LOOPS)))
         .unique(subset=["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]))
    if t.is_empty():
        return None, None, 0, 0

    # sorted, not merely grouped: an unordered key order permutes B's rows, which leaves every
    # symmetric functional here unchanged but makes the intermediate irreproducible to inspect.
    src = sorted({(c, int(r), lp) for c, r, lp in
                  zip(t["chain.id.from"], t["residue.index.from"], t["loop"])})
    dst = sorted({(c, int(r)) for c, r in zip(t["chain.id.to"], t["residue.index.to"])})
    ri = {(c, r): i for i, (c, r, _) in enumerate(src)}
    ci = {k: j for j, k in enumerate(dst)}
    B = np.zeros((len(src), len(dst)), dtype=np.int8)
    for cf, rf, ct, rt in zip(t["chain.id.from"], t["residue.index.from"],
                              t["chain.id.to"], t["residue.index.to"]):
        B[ri[(cf, int(rf))], ci[(ct, int(rt))]] = 1
    return B, [lp for _, _, lp in src], len(src), len(dst)


def _pielou(counts: np.ndarray, base: float) -> float:
    """Shannon entropy of ``counts`` as a distribution, divided by ``ln base``. NaN if undefined."""
    tot = float(counts.sum())
    if tot <= 0 or base <= 1:
        return float("nan")
    p = counts[counts > 0] / tot
    return float(-(p * np.log(p)).sum() / np.log(base))


def _algebraic_connectivity(B: np.ndarray) -> float:
    """Second-smallest eigenvalue of the normalised Laplacian on the largest connected component.

    Normalised, not combinatorial: ``L_sym = I - D^-1/2 A D^-1/2`` has its spectrum in ``[0, 2]``
    whatever the interface's size, so a large footprint and a small one are on one scale. Zero means
    the component is about to fall in two; the bipartite structure caps a connected component's
    largest eigenvalue at exactly 2, which is why the pair is worth reading together with
    ``g_comp_frac``.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    m, n = B.shape
    A = np.zeros((m + n, m + n))
    A[:m, m:] = B
    A[m:, :m] = B.T
    ncomp, lab = connected_components(coo_matrix(A), directed=False)
    keep = lab == np.bincount(lab).argmax()
    if keep.sum() < 3:
        return float("nan")
    sub = A[np.ix_(keep, keep)]
    deg = sub.sum(axis=1)
    if (deg <= 0).any():
        return float("nan")
    dinv = 1.0 / np.sqrt(deg)
    lsym = np.eye(len(sub)) - (sub * dinv[:, None]) * dinv[None, :]
    w = np.linalg.eigvalsh((lsym + lsym.T) / 2.0)     # symmetrise: eigvalsh reads one triangle
    return float(w[1])


def graph_features(structure: Structure, *, cutoff: float = 5.0) -> dict[str, float]:
    """Every graph functional of the 5 A contact map, as a flat row.

    Args:
        structure: a chain-typed, CDR-region-annotated TCR-pMHC structure. Run
            :func:`tcren.mhc.annotate_mhc` first or the MHC half of the graph is unreachable and
            every measure here is computed on peptide contacts alone.
        cutoff: heavy-atom contact threshold in Angstrom. The only parameter in the family.

    Returns:
        ``{feature: value}`` over :data:`GRAPH_FEATURES` and :data:`PROMOTED_POSE_FEATURES`. Values
        are ``nan`` where the structure gives them no support, never 0 -- an interface that makes no
        contact has no evenness, and reporting 0 would rank it below a bad one.
    """
    row: dict[str, float] = dict.fromkeys(GRAPH_FEATURES + PROMOTED_POSE_FEATURES, float("nan"))
    B, loops, m, n = _biadjacency(structure, cutoff)
    if B is None:
        return row

    a = B.sum(axis=1).astype(float)                    # TCR-side degrees
    b = B.sum(axis=0).astype(float)                    # pMHC-side degrees
    E = float(B.sum())
    V = m + n
    if E <= 0:
        return row

    # Evenness of the degree sequence on each side. Base is the engaged count, so this is evenness
    # proper and stays orthogonal to extent, which `n_pep_contacted` and `pep_cov_frac` already
    # carry. Undefined on one engaged residue -- one residue is not unevenly covered, it is the
    # whole cover -- so NaN rather than the 0 that `ln 1` would divide into.
    row["g_even_tcr"] = _pielou(a, m) if m > 1 else float("nan")
    row["g_even_pmhc"] = _pielou(b, n) if n > 1 else float("nan")

    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    A = np.zeros((V, V))
    A[:m, m:] = B
    A[m:, :m] = B.T
    ncomp = int(connected_components(coo_matrix(A), directed=False)[0])
    row["g_comp_frac"] = ncomp / V
    # For a graph, b1 = E - V + C exactly, so this is the hole count as a share of the contacts that
    # could carry it -- the size-free form the raw cyclomatic number is not.
    row["g_cyclo_frac"] = (E - V + ncomp) / E
    row["g_alg_conn"] = _algebraic_connectivity(B)

    # Per-loop partner sets: which pMHC residues each loop reaches. Counting distinct partners
    # rather than contacts is what makes this different from `H_loop`, whose tally rises with a
    # residue's atom count.
    partners: dict[str, set[int]] = {}
    for i, lp in enumerate(loops):
        partners.setdefault(lp, set()).update(np.flatnonzero(B[i]).tolist())
    engaged = [partners[k] for k in _CELL_LOOPS if partners.get(k)]
    if engaged:
        row["g_loop_even"] = _pielou(np.array([len(s) for s in engaged], float), len(_CELL_LOOPS))
    if len(engaged) > 1:
        jac = [len(x & y) / len(x | y) for i, x in enumerate(engaged) for y in engaged[i + 1:]]
        row["g_loop_overlap"] = float(np.mean(jac))

    # Degree assortativity over the E edges. Constant degree on either side leaves it undefined
    # rather than 0: there is no correlation to measure, which is not the same as no correlation.
    ei, ej = np.nonzero(B)
    if len(ei) > 2:
        x, y = a[ei], b[ej]
        if x.std() > 1e-12 and y.std() > 1e-12:
            row["g_assort"] = float(np.corrcoef(x, y)[0, 1])

    row.update(_promoted_degree(structure, cutoff))
    return row


def _promoted_degree(structure: Structure, cutoff: float) -> dict[str, float]:
    """:data:`PROMOTED_POSE_FEATURES`, through ``pose._degree_descriptors`` so there is one formula.

    That function counts on the receptor side of TCR:peptide, keyed by ``key.tcr.*``. The contact
    map puts the TCR on the ``from`` side of that interface by construction
    (:meth:`tcren.contactmap.ContactMap.interface`), so the rename is the whole adaptation and no
    orientation logic is repeated here.
    """
    from .pose import _degree_descriptors

    d = ContactMap.from_structure(structure, cutoff=cutoff).interface("tcr_peptide")
    got = _degree_descriptors(
        d.select(pl.col("chain.id.from").alias("key.tcr.chain"),
                 pl.col("residue.index.from").alias("key.tcr.res")), "_tp")
    return {k: got[k] for k in PROMOTED_POSE_FEATURES}


# --- the Calpha / Cbeta maps ---------------------------------------------------------------------
#
# Two readings, and the difference between them is measured rather than asserted. Both were screened
# against peptide length on 148 class I Native2026 crystals (peptide 8-13), reporting the rank
# correlation with length and the share of variance surviving after length is regressed out:
#
#     m_gap_tm             rho +0.029   100.0% beyond length   <- the cleanest column in the module
#     m_erank_tm           rho -0.188    93.1%
#     ca_cb_agreement_tp   rho +0.163    98.1%
#     m_face_tp / _tm      rho -0.146 / -0.065   97.7% / 99.3%
#     m_gap_tp             rho -0.306    94.4%
#     ca_cb_agreement_tm   rho -0.354    83.5%
#     m_erank_tp           rho -0.518    73.1%
#
# The peptide-facing spectral columns carry length, and the mechanism is the author's (2026-09-01):
# in class I the groove is closed at both ends, so a longer peptide must BULGE, and a bulged peptide
# sits closer to the receptor and presents more independent approach modes. That is a real property
# an end user asks about -- is this a bulged epitope -- not an artefact to normalise away. It is
# also why none of them is dropped: every one keeps most of its variance after length is removed
# (73% in the worst case), so each carries something length does not, and a model given both a
# length-coupled and a length-free column can form the contrast that cancels the shared part. The
# `STATUS` entries say which is which so the choice is the caller's.

def _map_features(structure: Structure, cutoff: float) -> dict[str, float]:
    """The Calpha-against-Cbeta comparison on both receptor arms.

    Reads :func:`tcren.pose._interface_layers`, which stacks the three representative-atom layers
    -- ``d1`` closest heavy atom, ``d2`` Cbeta, ``d3`` Calpha -- onto one row per residue pair. Going
    through it rather than re-deriving the layers is what keeps ``m_face_tp`` one number with one
    implementation: :func:`tcren.pose.pose_consistency` reads the same column and the two agree to
    4.4e-16 on 196 crystals, which a test pins.

    ``m_face_*`` is the mean over **contacting** pairs; ``ca_cb_agreement_*`` is the rank correlation
    over the whole **approach shell**, which is about ten times larger and includes the pairs that
    come close and form nothing.
    """
    from ..contacts.definitions import _interface_layers
    from .pose import _CA_SHELL_RADIUS, _spearman
    from ..structure.model import MHC_TYPES as _MHC

    out: dict[str, float] = dict.fromkeys(MATRIX_FEATURES, float("nan"))
    for suffix, partner in (("_tp", (PEPTIDE_TYPE,)), ("_tm", tuple(_MHC))):
        try:
            wide = _interface_layers(structure, cutoff, partner=partner)
        except Exception:  # noqa: BLE001 - an arm with no partner chain is not an error
            continue
        if wide.is_empty():
            continue
        # face: the contacting pairs only. A pair whose backbones are close while its side chains
        # point away reads negative, which is the shape a contact-count objective produces.
        both = wide.filter(pl.col("d1").is_not_null()
                           & pl.col("d2").is_not_null() & pl.col("d3").is_not_null())
        if both.height:
            out[f"m_face{suffix}"] = float(
                (both["d3"].to_numpy() - both["d2"].to_numpy()).mean())
        # agreement: the whole approach shell. Do the side chains track the backbone, as they do
        # in a real interface, or has the packing been solved for one and not the other?
        near = wide.filter(pl.col("d3").is_not_null() & (pl.col("d3") <= _CA_SHELL_RADIUS)
                           & pl.col("d2").is_not_null())
        if near.height >= 3:
            out[f"ca_cb_agreement{suffix}"] = _spearman(near["d3"].to_numpy(),
                                                        near["d2"].to_numpy())
    return out


def _region_residues(structure: Structure, chain_types, region_types=None) -> list:
    """Residues of the named chain types, restricted to the named region markup when given."""
    out = []
    for c in structure.chains:
        if c.chain_type not in chain_types:
            continue
        if region_types is None:
            out.extend(c.residues)
            continue
        for reg in getattr(c, "regions", []) or []:
            if reg.region_type in region_types:
                out.extend(reg.residues)
    return out


def _spectral(rows: list, cols: list, sigma: float) -> tuple[float, float]:
    """``(effective rank fraction, s2/s1)`` of the Gaussian proximity kernel between two residue sets.

    ``K = exp(-D^2 / 2 sigma^2)`` over the Calpha distance map. The effective rank is the
    participation ratio of the singular values divided by ``min(L, M)``; the gap is ``s2/s1``, which
    is 0 for a separable (rank-one) approach -- a receptor leaning on a surface rather than pairing
    specific residues with it.
    """
    ra = np.array([r.ca for r in rows if r.ca is not None], float)
    ca = np.array([r.ca for r in cols if r.ca is not None], float)
    if len(ra) < 2 or len(ca) < 2:
        return float("nan"), float("nan")
    d = np.linalg.norm(ra[:, None, :] - ca[None, :, :], axis=2)
    s = np.linalg.svd(np.exp(-(d ** 2) / (2.0 * sigma ** 2)), compute_uv=False)
    s = s[s > 0]
    p = min(len(ra), len(ca))
    if len(s) < 2 or p < 2:
        return float("nan"), float("nan")
    return float(s.sum() ** 2 / (s ** 2).sum() / p), float(s[1] / s[0])


def matrix_features(structure: Structure, *, cutoff: float = 5.0,
                    sigma: float = _SIGMA) -> dict[str, float]:
    """Length-agnostic comparisons of the interface Calpha and Cbeta maps.

    Args:
        structure: a chain-typed TCR-pMHC structure. The ``_tm`` columns need
            :func:`tcren.mhc.annotate_mhc` to have run, since ``MHC_TYPES`` matches the supertype
            it assigns; without it they are ``nan``.
        cutoff: heavy-atom contact threshold in Angstrom, for the ``m_face_*`` pair.

    Returns:
        ``{feature: value}`` over :data:`MATRIX_FEATURES`, ``nan`` where unsupported.

    Note:
        Measured against peptide length on 196 Native2026 crystals (peptide 6-20 residues):
        ``m_face_tm`` +0.023, ``ca_cb_agreement_tp`` -0.043, ``m_face_tp`` +0.161,
        ``ca_cb_agreement_tm`` -0.272. The first two are the cleanest columns in this module.
    """
    out = _map_features(structure, cutoff)
    loops = _region_residues(structure, ("TRA", "TRB"), ("CDR1", "CDR2", "CDR3"))
    out["m_erank_tp"], out["m_gap_tp"] = _spectral(
        loops, _region_residues(structure, (PEPTIDE_TYPE,)), sigma)
    out["m_erank_tm"], out["m_gap_tm"] = _spectral(
        loops, _region_residues(structure, MHC_TYPES, MHC_HELIX_REGIONS), sigma)
    return out


def _selfcheck() -> None:  # pragma: no cover - exercised by tests/unit/test_interface_graph.py
    """Assert the invariants that make these numbers meaningful."""
    # Pielou: uniform is exactly 1, concentrated is below it, and the base is what is passed.
    assert abs(_pielou(np.array([5.0, 5, 5, 5]), 4) - 1.0) < 1e-12
    assert _pielou(np.array([17.0, 1, 1, 1]), 4) < 0.7
    assert np.isnan(_pielou(np.zeros(4), 4))

    # A perfect matching (each loop residue reads exactly one distinct partner) is one component
    # per edge and carries no cycle; a complete bipartite block is one component and E - V + 1.
    eye = np.eye(4, dtype=np.int8)
    assert _biadj_cyclo(eye) == 0.0, _biadj_cyclo(eye)
    full = np.ones((3, 4), dtype=np.int8)
    assert abs(_biadj_cyclo(full) - (12 - 7 + 1) / 12) < 1e-12, _biadj_cyclo(full)

    # Algebraic connectivity: a complete bipartite graph is maximally connected (lambda2 = 1 for
    # K_{m,n} under the normalised Laplacian); a barbell of two blocks joined by one edge is near 0.
    assert abs(_algebraic_connectivity(full) - 1.0) < 1e-9, _algebraic_connectivity(full)
    bar = np.zeros((4, 4), dtype=np.int8)
    bar[:2, :2] = 1
    bar[2:, 2:] = 1
    bar[1, 2] = 1
    assert 0.0 < _algebraic_connectivity(bar) < 0.5, _algebraic_connectivity(bar)



def _biadj_cyclo(B: np.ndarray) -> float:
    """``(E - V + C) / E`` for a biadjacency matrix; the self-check's route to the same algebra."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    m, n = B.shape
    A = np.zeros((m + n, m + n))
    A[:m, m:] = B
    A[m:, :m] = B.T
    c = int(connected_components(coo_matrix(A), directed=False)[0])
    e = float(B.sum())
    return (e - (m + n) + c) / e


if __name__ == "__main__":  # pragma: no cover
    _selfcheck()
    print("interface_graph: self-check passed")
