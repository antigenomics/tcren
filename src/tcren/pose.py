"""Per-structure pose consistency: do the tight contacts carry the favourable chemistry?

:func:`tcren.cohort.coupling` measures the forced-pose signature **across a cohort** as
``C* = corr(Q, dPhi)``: in a genuine complex a better interface holds more favourable contacts, so
the two channels rise together, while a generator that manufactures a pose optimises contacts
without the interface and breaks that tie. It is the right diagnostic and the wrong estimator for a
user with two or three models --- at ``n = 2`` the sample correlation is +-1 by construction, and its
sign is wrong in roughly a third to a half of draws (``bench/scripts/coupling_smalln.py``).

The tie it measures also holds **within one structure**, over that structure's own contacts. In a
crystal the residue pairs that sit tightest are the ones whose identities are complementary, because
that is what selected the pose; a pose built to satisfy a contact-density prior has no such
alignment. Correlating contact tightness against contact favourability *inside* a single complex
therefore reads the same physics from ``n = 1`` structure, over its ~20--120 interface pairs.

Three superimposable maps over one residue-pair index carry it, and
:func:`tcren.contacts.multi_contacts` already returns all three stacked (``layer`` column):

* ``d1`` --- closest heavy-atom distance (the 5 A contact definition used everywhere else);
* ``d2`` --- Cbeta distance (Calpha for glycine), i.e. where the side chains point;
* ``d3`` --- Calpha distance, i.e. where the backbones sit.

The chemistry axis is ``J``, not the raw potential entry. A contact energy splits as
``e(a,b) = mean + H_tcr(a) + H_pep(b) + J(a,b)``, where the two ``H`` terms depend on one residue
each and ``J`` is the double-centred remainder. Correlating distance against raw ``e`` would partly
measure *which* residues happen to sit at the interface rather than whether they suit each other;
``J`` is the pair-specific part, and complementarity lives there.

:meth:`tcren.Potential.decompose` performs that split but only for a *symmetric* matrix, and TCRen2
is deliberately directional (TCR residue by peptide residue; symmetrising it costs measurable
accuracy). :func:`_double_centred` therefore applies the same two-way centring to the matrix as
given, which is well defined whether or not it is symmetric and reduces to ``decompose`` when it is.

Every descriptor is oriented **higher = more crystal-like**, so they compose with
:func:`tcren.cohort.q_score` under the same all-descriptors-higher-is-better convention.

Evaluation (ROC/PR/CI) belongs downstream in the benchmark repo, not here.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from functools import lru_cache

from .contacts.definitions import ContactDefinition, multi_contacts
from .structure.model import PEPTIDE_TYPE, RECEPTOR_TYPES, Structure

__all__ = ["pose_consistency", "c_score", "pose_native_reference", "pose_af_reference",
           "POSE_FEATURES"]

#: The cross-map descriptors :func:`pose_consistency` returns, each oriented positive-is-crystal-like.
#: These are the ``k`` terms a pose score standardizes against a native-crystal reference.
POSE_FEATURES = (
    "c_local",
    "e_tight_minus_loose",
    "frac_close_favourable",
    "frac_cb_close_engaged",
    "sidechain_toward",
    "margin_energy_slope",
)

# The Cbeta layer threshold for "side chains are near enough that an interaction is expected".
# Distinct from the layer build cutoff below, which is deliberately generous so that every d1
# contact also has a Cbeta and a Calpha distance available to pair against.
_CB_CLOSE = 8.0
# Long side chains (Arg, Lys, Trp) put two heavy atoms within 5 A while their Calpha atoms sit far
# apart, so the representative-atom layers are built well past their nominal 8/12 A defaults; the
# per-descriptor thresholds are applied afterwards by filtering.
_REP_BUILD_CUTOFF = 18.0

_KEY = ["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]


def _double_centred(potential):
    """``(J, index)`` --- the pair-specific part of a potential, one-body terms removed.

    ``J(a,b) = e(a,b) - rowmean(a) - colmean(b) + grandmean`` over the residue axes, so every row and
    column of ``J`` sums to zero and what remains depends on the *pair*, not on either residue alone.
    Unlike :meth:`tcren.Potential.decompose` this does not require a symmetric matrix, which matters
    because the shipped TCRen2 is directional. Never-observed cells are ``nan`` and are ignored by
    the means (and stay ``nan`` in ``J``, so they drop out of every descriptor).
    """
    import warnings

    m, index = potential.as_matrix()
    with warnings.catch_warnings():
        # An alphabet entry with no observed contact at all (e.g. the gap symbol) is an all-nan
        # row: its mean is genuinely undefined, J stays nan there, and those pairs drop out.
        warnings.simplefilter("ignore", RuntimeWarning)
        row = np.nanmean(m, axis=1, keepdims=True)
        col = np.nanmean(m, axis=0, keepdims=True)
        return m - row - col + np.nanmean(m), index


def _pair_j(aa_from, aa_to, jmat, index) -> np.ndarray:
    """Vectorised gather of ``J(a, b)``; NaN outside the alphabet or for an unobserved cell."""
    i = np.array([index.get(a, -1) for a in aa_from], dtype=np.int64)
    j = np.array([index.get(b, -1) for b in aa_to], dtype=np.int64)
    out = np.full(len(i), np.nan)
    ok = (i >= 0) & (j >= 0)
    out[ok] = jmat[i[ok], j[ok]]
    return out


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rho over the finite pairs; NaN when fewer than 3 remain or either side is constant."""
    from scipy.stats import spearmanr

    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.std(x[ok]) < 1e-12 or np.std(y[ok]) < 1e-12:
        return float("nan")  # guard before the call: scipy warns rather than returning quietly
    return float(spearmanr(x[ok], y[ok]).statistic)


def _interface_layers(structure: Structure, cutoff: float) -> pl.DataFrame:
    """The d1/d2/d3 layers pivoted onto one row per TCR:peptide residue pair.

    Returns a frame keyed by the residue pair with columns ``d1``/``d2``/``d3`` (Angstrom, null where
    that layer does not reach), ``aa.tcr``/``aa.pep`` and the d1 atom names. The TCR side is
    normalised to ``aa.tcr`` regardless of which side the canonical chain ordering put it on.
    """
    stacked = multi_contacts(
        structure,
        ContactDefinition(d1=cutoff, d2=_REP_BUILD_CUTOFF, d3=_REP_BUILD_CUTOFF),
    )
    ctype = {c.chain_id: c.chain_type for c in structure.chains}
    stacked = stacked.with_columns(
        pl.col("chain.id.from").replace_strict(ctype, default=None).alias("type.from"),
        pl.col("chain.id.to").replace_strict(ctype, default=None).alias("type.to"),
    )
    tcr, pep = list(RECEPTOR_TYPES), [PEPTIDE_TYPE]
    fwd = pl.col("type.from").is_in(tcr) & pl.col("type.to").is_in(pep)
    rev = pl.col("type.from").is_in(pep) & pl.col("type.to").is_in(tcr)
    stacked = stacked.filter(fwd | rev).with_columns(
        pl.when(fwd).then(pl.col("residue.aa.from")).otherwise(pl.col("residue.aa.to")).alias("aa.tcr"),
        pl.when(fwd).then(pl.col("residue.aa.to")).otherwise(pl.col("residue.aa.from")).alias("aa.pep"),
    )
    if stacked.is_empty():
        return stacked.select(*_KEY, "aa.tcr", "aa.pep").with_columns(
            pl.lit(None, dtype=pl.Float64).alias(c) for c in ("d1", "d2", "d3")
        )
    # Residue identity is carried by the d1 layer alone; d2/d3 contribute only their distance, so
    # the three frames share no column but the key and nothing collides on the join. A pair present
    # in d2 but not d1 still appears (outer join) --- that is exactly the unengaged pair
    # `frac_cb_close_engaged` counts, and it needs no residue identity.
    wide = (stacked.filter(pl.col("layer") == "d1")
            .select(*_KEY, "aa.tcr", "aa.pep", pl.col("dist").alias("d1")))
    for layer in ("d2", "d3"):
        part = (stacked.filter(pl.col("layer") == layer)
                .select(*_KEY, pl.col("dist").alias(layer)))
        # A residue pair appears at most once per layer (each keeps its closest atom pair).
        wide = wide.join(part, on=_KEY, how="full", coalesce=True)
    return wide


def pose_consistency(
    structure: Structure, potential=None, cutoff: float = 5.0
) -> dict[str, float]:
    """Cross-map consistency descriptors of one TCR:peptide interface.

    Reads whether the structure's *tight* contacts are its *complementary* ones --- the
    within-structure analogue of :func:`tcren.cohort.coupling`, and unlike it defined for a single
    complex. Every value is oriented so that **higher is more crystal-like**.

    Args:
        structure: a chain-typed complex (:func:`tcren.annotation.classify_chains` run) with a
            peptide chain and at least one receptor chain.
        potential: the residue-pair potential whose double-centred ``J`` supplies the chemistry
            axis; defaults to the bundled TCRen2 matrix.
        cutoff: the heavy-atom contact cutoff (A) defining the d1 layer and the contact margin.

    Returns:
        A dict with :data:`POSE_FEATURES` plus ``n_contacts`` (the pair count every value rests on)
        and ``n_cb_close``. Descriptors that cannot be estimated --- fewer than three contacts, a
        constant axis, no Cbeta-close pairs --- come back as ``nan`` rather than a made-up number.

    Note:
        This is a *pose* readout, not a binder score: it says whether the geometry and the chemistry
        of one model agree, not whether the receptor binds.
    """
    if potential is None:
        from .potential import tcren2

        potential = tcren2()
    jmat, jindex = _double_centred(potential)

    wide = _interface_layers(structure, cutoff)
    contacts = wide.filter(pl.col("d1").is_not_null())
    n = contacts.height
    out: dict[str, float] = {k: float("nan") for k in POSE_FEATURES}
    out["n_contacts"] = float(n)
    out["n_cb_close"] = float("nan")

    # --- the Cbeta-engagement descriptor lives on the d2 layer, not on the contacts ---------------
    cb_close = wide.filter(pl.col("d2").is_not_null() & (pl.col("d2") <= _CB_CLOSE))
    out["n_cb_close"] = float(cb_close.height)
    if cb_close.height:
        out["frac_cb_close_engaged"] = float(cb_close["d1"].is_not_null().mean())

    # --- do the side chains lean in? Cbeta closer than Calpha means they point at each other ------
    # A mean over pairs, so unlike the correlations below it is defined for a single contact.
    both = contacts.filter(pl.col("d2").is_not_null() & pl.col("d3").is_not_null())
    if both.height:
        out["sidechain_toward"] = float(
            (both["d3"].to_numpy() - both["d2"].to_numpy()).mean()
        )

    # Everything below is a correlation, a slope or a tercile split, and needs at least three pairs.
    if n < 3:
        return out

    d1 = contacts["d1"].to_numpy()
    j = _pair_j(contacts["aa.tcr"].to_list(), contacts["aa.pep"].to_list(), jmat, jindex)
    margin = cutoff - d1          # positional slack: higher = the pair sits deeper than the cutoff
    fav = -j                      # favourability: J is an energy, so lower J is better

    out["c_local"] = _spearman(margin, fav)

    ok = np.isfinite(fav)
    if ok.sum() >= 3 and np.std(margin[ok]) > 1e-12:
        # per-Angstrom slope of favourability on slack; the signed, physically-scaled companion
        out["margin_energy_slope"] = float(np.polyfit(margin[ok], fav[ok], 1)[0])

    if ok.sum() >= 3:
        lo, hi = np.quantile(d1[ok], [1 / 3, 2 / 3])
        tight, loose = fav[ok][d1[ok] <= lo], fav[ok][d1[ok] >= hi]
        if len(tight) and len(loose):
            out["e_tight_minus_loose"] = float(tight.mean() - loose.mean())
        med = np.median(d1[ok])
        close = j[ok][d1[ok] <= med]
        if len(close):
            out["frac_close_favourable"] = float((close < 0).mean())
    return out


def _selfcheck() -> None:
    """Assert the descriptors read the sign they claim, on a hand-built two-pair interface."""
    index = {"A": 0, "B": 1, "C": 2}
    jm = np.array([[0.0, -1.0, 1.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    # three pairs spanning J = -1 (complementary), 0, +1 (repulsive)
    j = _pair_j(["A", "B", "A"], ["B", "C", "C"], jm, index)
    assert np.allclose(j, [-1.0, 0.0, 1.0]), j
    # crystal-like: the complementary pair sits tightest, the repulsive one loosest
    rho = _spearman(5.0 - np.array([3.0, 4.0, 4.9]), -j)
    assert rho == 1.0, f"tight+favourable must give c_local = +1, got {rho}"
    # forced: the same chemistry with the distance order inverted
    rho_forced = _spearman(5.0 - np.array([4.9, 4.0, 3.0]), -j)
    assert rho_forced == -1.0, f"inverted pose must give c_local = -1, got {rho_forced}"
    # an unknown residue must not silently score as zero
    assert np.isnan(_pair_j(["A"], ["X"], jm, index)[0])
    assert np.isnan(_spearman(np.array([1.0, 2.0]), np.array([1.0, 2.0])))  # n < 3
    assert np.isnan(_spearman(np.ones(5), np.arange(5.0)))                  # constant axis

    # double-centring: rows and columns of J sum to zero, for an ASYMMETRIC matrix too
    class _P:
        _m = np.array([[1.0, 2.0, 9.0], [3.0, 0.0, 1.0], [5.0, 4.0, 2.0]])
        def as_matrix(self):
            return self._m, {"A": 0, "B": 1, "C": 2}

    J, _ = _double_centred(_P())
    assert np.allclose(J.sum(axis=0), 0) and np.allclose(J.sum(axis=1), 0), J
    # a nan cell must not poison the whole row
    class _Pn(_P):
        _m = np.array([[1.0, 2.0, np.nan], [3.0, 0.0, 1.0], [5.0, 4.0, 2.0]])

    Jn, _ = _double_centred(_Pn())
    assert np.isnan(Jn[0, 2]) and np.isfinite(Jn[0, 0])
    print("pose selfcheck ok")


if __name__ == "__main__":
    _selfcheck()


def _load_reference(name: str) -> dict:
    """Load a bundled reference CSV as a dict of column arrays."""
    import csv
    from importlib import resources

    path = resources.files("tcren.data") / name
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    cols = [c for c in rows[0] if c != "pdb.id"]
    return {c: np.array([float(r[c]) for r in rows]) for c in cols}


@lru_cache(maxsize=1)
def pose_af_reference() -> dict:
    """The pose descriptors over 1,018 **AlphaFold** TCR:pMHC models --- the reference to use when
    the structure being scored is itself generated.

    Every generated pose carries a constant offset from the crystal manifold (median ``c_local``
    0.031 here against 0.085 on crystals, ``frac_close_favourable`` 0.651 against 0.764). Scored
    against :func:`pose_native_reference`, that shared offset dominates and ``C`` reads *provenance*
    --- how model-like the structure is --- which is the right question for spotting a fabricated
    complex and the wrong one for ranking models against each other. Standardizing against the
    generated manifold instead removes the offset, so ``C`` reads which model is an outlier *among
    its own kind*.

    Fitted **label-blind** on the whole ``vdjdb_binder_benchmark`` deposit (523 real + 566 mock
    pairings): it defines what an AlphaFold TCR:pMHC model looks like, not what a binder looks like,
    so no binder label enters it. Disjoint from the TCRvdb cohort (0 of 618 ids shared).

    Provenance: ``scripts/fit_pose_reference.py --struct-dir <vdjdb_binder_benchmark>``.
    """
    return _load_reference("pose_af_reference.csv")


@lru_cache(maxsize=1)
def pose_native_reference() -> dict:
    """The pose descriptors over the Native2026 alpha-beta crystals, bundled so a **single** user
    structure can be standardized against the natural interface manifold rather than against itself.

    Returns a dict of column arrays usable as the ``reference`` argument of :func:`c_score` (and of
    :func:`tcren.cohort.zscore`). Provenance: ``scripts/fit_pose_reference.py`` over
    ``data/Native2026``, restricted by the derivation's hard rule --- both CDR3 loops resolved and a
    20-letter peptide --- and to rows finite on every descriptor.

    This is a **second** reference file. :func:`tcren.cohort.native_reference` is untouched, so no
    published ``Q`` moves.
    """
    return _load_reference("pose_native_reference.csv")


def c_score(table, reference=None, features=POSE_FEATURES, method="z", decorrelate=True):
    r"""Pose-consistency score ``C`` --- fit-free, single-structure-capable; higher = more crystal-like.

    The same construction as :func:`tcren.cohort.q_score`, on the cross-map descriptors instead of
    the interface-quality ones:

    .. math::  C(x) = z(x)^{\top} \hat{C}^{-1} \mathbf{1},
       \qquad z(x)_k = \frac{d_k(x) - \mu_k}{\sigma_k}

    with :math:`\mu, \sigma` and the descriptor covariance estimated on the **native crystal
    reference** (:func:`pose_native_reference`), and :math:`\mathbf 1` the fixed
    every-descriptor-higher-is-more-crystal-like direction. No fitted coefficient, no negative set,
    defined for a single row.

    ``C`` answers "do this model's geometry and chemistry agree?", which is *not* "does this receptor
    bind". It is the pose channel: compose it with a binder score, do not substitute it for one.

    Args:
        table: a mapping / polars / pandas frame carrying :data:`POSE_FEATURES` --- e.g. rows built
            from :func:`pose_consistency`.
        reference: cohort defining :math:`\mu, \sigma` and the covariance; defaults to
            :func:`pose_native_reference`.
        features: the descriptors to combine; defaults to all of :data:`POSE_FEATURES`.
        method: per-descriptor standardization, ``"z"`` (default) or ``"rank"``.
        decorrelate: whiten by the native covariance (default); ``False`` gives the equal-weight mean.
    """
    from .cohort import q_score

    if reference is None or isinstance(reference, str):
        key = reference or "native"
        try:
            reference = {"native": pose_native_reference, "af": pose_af_reference}[key]()
        except KeyError:
            raise ValueError(f"reference must be 'native', 'af' or a mapping; got {key!r}") from None
    return q_score(table, reference, features=features, method=method, decorrelate=decorrelate)
