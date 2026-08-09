"""Cohort-relative recognition scores — the **recommended, fit-free** screening layer.

Prefer these over the fitted :func:`tcren.binder.binder_score` (``p_bind``) and
:func:`tcren.recognition.forced_pose_score` (``p_forced``). Those carry trained coefficients; the
functions here carry none — no logistic, no fit, no training set — so they cannot leak, cannot go
stale, and there is nothing to re-derive. The benchmark repo settled the trade-off empirically
(ledger C24/C25/C26):

* :func:`q_score` matches or beats the fitted ``p_bind`` and, unlike it, **generalises across
  cohorts** — a logistic trained on one cohort learns that cohort's epitope composition and does not
  transfer, whereas ``Q`` has nothing to transfer. With ipTM it reproduces the headline synergy
  fit-free: ``z(ipTM) + z(Q)`` reaches macro ROC 0.83 on TCRvdb against ipTM's 0.79.
* :func:`strain_z` grades pose forcedness (crystal < AF-real < AF-decoy) reproducibly, unlike
  ``FORCED_POSE_MODEL`` whose training rows are lost.

They are **cohort-relative** by default: each standardizes a feature over *the set being ranked*.
For a candidate set, score the whole batch together (``tcren recognize`` over a directory). For a
**single structure**, or a small/heterogeneous user set where the batch is not a fair reference, pass
``reference=native_reference()`` (with ``features=Q_FEATURES_GEOM``): the descriptors are then
standardized against the shipped Native2026 crystal manifold, so ``Q`` is defined for one structure and
transfers across inputs. The descriptors are counts and bounded ratios (mildly non-normal), so
``method="rank"`` gives a robust, assumption-free percentile standardization; on the benchmarks it
agrees with the default ``z`` to ρ≈0.98. The division of labour is scores in ``tcren``, evaluation
(ROC/PR/CI) downstream.

All functions take the table ``tcren recognize --full`` emits (a mapping of column name to
sequence, a ``polars``/``pandas`` frame, or a dict of arrays) and return one value per row.

Sign convention: every term is oriented so that **higher = more binder-like** for
:func:`q_score`, and **higher = more forced/strained** for :func:`strain_z`.

.. note::
   :func:`phi_bind` is **deprecated** — every term it adds to ``Q`` lowers ranking accuracy
   (benchmark ledger C19b), and its ``z(-pitch)`` term is both below chance on its own and derived
   from an AlphaFold-contaminated angle. Use :func:`q_score`.
"""

from __future__ import annotations

import warnings
from functools import lru_cache

import numpy as np

__all__ = ["zscore", "q_score", "q_iptm", "f_score", "q_f", "q_f_iptm", "f_invert_by_iptm", "phi_bind",
           "q_coupled", "coupling", "strain_z", "native_reference", "Q_FEATURES", "Q_FEATURES_CORE",
           "Q_FEATURES_GEOM", "F_TERMS", "STRAIN_TERMS"]

#: The five interface-quality descriptors, equal-weighted in :func:`q_score`. Each is oriented
#: positive-is-better as given. ``pp_combo`` is the CDR1/2-vs-CDR3alpha TCRen contrast — the one
#: energy term robust to the forced-pose inversion (benchmark ledger C27), since it is a *contrast*
#: rather than an absolute contact energy. Per-term macro AUROC on TCRvdb: burial 0.73, n_hbond 0.69,
#: pp_combo 0.66, n_pep_contacted 0.62, chain_balance 0.61; the terms are near-independent
#: (mean absolute Spearman 0.20).
Q_FEATURES = ("burial", "n_pep_contacted", "chain_balance", "n_hbond", "pp_combo")

#: The four load-bearing descriptors. ``n_pep_contacted`` is dropped: it is the weakest term and
#: removing it *raises* macro AUROC 0.795 -> 0.801 on TCRvdb (benchmark ledger, energy memo). Pass
#: ``features=Q_FEATURES_CORE`` to :func:`q_score` for the simpler, marginally better score.
Q_FEATURES_CORE = ("burial", "chain_balance", "n_hbond", "pp_combo")

#: The four **geometry-only** descriptors — :data:`Q_FEATURES` without the ``pp_combo`` energy contrast.
#: This is ``Q_geom``, the AF-orthogonal channel that survives the forced-pose regime where the contact
#: energy inverts (benchmark ledger C27/C42): ``z(ipTM) + z(q_score(..., features=Q_FEATURES_GEOM))``
#: beats raw-AF ipTM on well-modelled ("template-covered") epitopes on both ROC and PR, while the energy
#: term is used only conditioned on pose quality. Pass to :func:`q_score` / :func:`q_iptm`.
Q_FEATURES_GEOM = ("burial", "n_pep_contacted", "chain_balance", "n_hbond")

#: The TCRen contact-energy terms summed into the binder-oriented :func:`f_score`. ``F_tcr_pep`` is the
#: TCR:peptide TCRen energy, ``F_tcr_mhc`` the TCR:MHC energy; both are emitted by ``tcren recognize``.
#: They are raw energies (lower = tighter), so :func:`f_score` negates the sum to make higher = more
#: binder-like. **This term is pose-conditional** — it reads real binding chemistry on well-modelled
#: (crystal-templated) poses and *inverts* on forced ones (benchmark ledger C27/C42): on the forced
#: GLCTLVAML TCRvdb pose ``-F_tcr_pep`` ranks binders at AUROC 0.36 (backwards), on the clean YLQPRTFLL
#: pose at 0.59. Use it only conditioned on pose quality — gate with :func:`strain_z`, or read
#: ``z(Q)-z(F)`` on forced poses and ``z(Q)+z(F)`` on clean ones.
F_TERMS = ("F_tcr_pep", "F_tcr_mhc")

#: Crystal-calibrated interface-strain terms with their physical signs. A forced pose reaches
#: further from the peptide with a thinner, less balanced interface.
STRAIN_TERMS = (("cdr3b_topep", +1.0), ("cdr3b_reach", +1.0),
                ("extent_per_ct", +1.0), ("chain_balance", -1.0))


def _col(table, name):
    """Fetch a column from a dict / pandas / polars frame as a float array."""
    if hasattr(table, "columns") and not isinstance(table, dict):
        if name not in table.columns:
            raise KeyError(f"column {name!r} not in table; build it with "
                           f"recognition_table(items, full=True) or `tcren recognize --full`")
        col = table[name]
        return np.asarray(col.to_numpy() if hasattr(col, "to_numpy") else col, float)
    if name not in table:
        raise KeyError(f"column {name!r} not in table; build it with "
                       f"recognition_table(items, full=True) or `tcren recognize --full`")
    return np.asarray(table[name], float)


def _derive(table, name):
    """Columns that ``recognize`` does not emit directly but are one division away."""
    if name == "extent_per_ct":  # interface thinness
        return _col(table, "extent") / np.maximum(_col(table, "n_contacts_tp"), 1.0)
    if name == "pp_combo":       # z(sum J CDR1/2) - z(sum J CDR3alpha)
        return zscore(_col(table, "F_cdr12")) - zscore(_col(table, "F_cdr3a"))
    return _col(table, name)


def zscore(x, reference=None, method="z") -> np.ndarray:
    """NaN-aware standardization. ``reference`` calibrates against another cohort.

    Passing ``reference`` is what makes :func:`strain_z` *crystal-calibrated*: the mean and sd come
    from the crystallographic ensemble, so the score reads ~0 on crystals by construction and grows
    as a pose departs from the natural manifold. Without it, a cohort of uniformly forced poses
    would standardize to zero mean and the shift would be invisible.

    Args:
        x: values to standardize.
        reference: cohort defining the location/scale; defaults to ``x`` itself (cohort-relative).
        method: ``"z"`` (mean/sd, the default) or ``"rank"`` — the percentile of each ``x`` against
            the reference, mapped to ``[-1, 1]`` (``2·percentile − 1``). ``"rank"`` is scale-free and
            makes **no normality assumption**, so it is the robust choice for the bounded/count
            descriptors of ``Q`` (chain balance, H-bond and contact counts are not normal). On the
            benchmarks ``z`` and ``rank`` agree to Spearman ρ≈0.98 and differ by <0.005 AUROC, so
            ``z`` is kept as the default; use ``rank`` when a heavy-tailed user descriptor could
            distort the mean/sd.
    """
    x = np.asarray(x, float)
    ref = x if reference is None else np.asarray(reference, float)
    if method == "rank":
        r = np.sort(ref[np.isfinite(ref)])
        if r.size == 0:
            return np.full_like(x, np.nan)
        pct = np.searchsorted(r, x, side="right") / r.size   # fraction of reference <= x
        return np.where(np.isfinite(x), 2.0 * pct - 1.0, np.nan)
    mu = np.nanmean(ref)
    sd = np.nanstd(ref)
    # A constant column does not give sd == 0 exactly: np.nanstd(np.full(20, 3.7)) is 4.4e-16.
    # Testing `sd > 0` therefore divides by float residue and amplifies noise by ~1e16. Scale the
    # tolerance to the data so a genuinely constant (or degenerate) descriptor contributes nothing.
    if not np.isfinite(sd) or sd <= 1e-12 * max(1.0, abs(float(mu))):
        return np.zeros_like(x)
    return (x - mu) / sd


@lru_cache(maxsize=1)
def native_reference() -> dict:
    """The interface-geometry descriptors over the 374 Native2026 crystal complexes, shipped so a
    **single** user structure (or any small cohort) can be standardized against the natural interface
    manifold instead of against itself — the deployment path for generic input::

        from tcren import cohort
        q = cohort.q_score(user_table, reference=cohort.native_reference(),
                           features=cohort.Q_FEATURES_GEOM)

    Use :data:`Q_FEATURES_GEOM` (the four geometry terms) for one structure: the fifth term
    ``pp_combo`` is a within-cohort z-contrast and is undefined for a single row. Returns a dict of
    column arrays (``burial, n_pep_contacted, chain_balance, n_hbond, F_cdr12, F_cdr3a``) usable as
    the ``reference`` argument. Provenance: ``tcren recognize --full`` over the Native2026 set.
    """
    import csv
    from importlib import resources
    path = resources.files("tcren.data") / "q_native_reference.csv"
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    return {c: np.array([float(r[c]) for r in rows]) for c in rows[0]}


def q_score(table, reference=None, features=Q_FEATURES_GEOM, method="z", decorrelate=True) -> np.ndarray:
    r"""Interface-quality score ``Q`` — fit-free, single-structure-capable; the default binder score.

    The default is the **directional, decorrelated** one-class score over ``k = len(features)`` terms

    .. math::  Q(x) \;=\; z(x)^{\top}\, C^{-1}\, \mathbf{1},
       \qquad z(x)_k = \frac{d_k(x)-\mu_k}{\sigma_k},

    where each descriptor is standardized against the **native crystal reference**
    (``reference``, default :func:`native_reference` — its :math:`\mu_k,\sigma_k`), :math:`C` is the
    native descriptor correlation matrix, and :math:`\mathbf 1` is the biophysical
    *every-descriptor-higher-is-better* direction. Whitening by :math:`C^{-1}` stops correlated
    descriptors from double-counting. The score carries **no fitted coefficient** (only the native
    covariance is estimated; :math:`\mathbf 1` is fixed), needs **no negative set**, is calibrated on
    natives so it **transfers** across inputs, and is defined for a **single structure**. It reduces to
    the equal-weight mean when the descriptors are uncorrelated (:math:`C=I`). See the manuscript
    Methods (\S Scores).

    ``decorrelate=False`` recovers the legacy equal-weight mean :math:`Q=\frac1k\sum_k z(d_k)`.

    Args:
        table: the ``tcren recognize --full`` table (dict / pandas / polars).
        reference: cohort defining :math:`\mu,\sigma,C`. ``None`` uses :func:`native_reference` when
            ``decorrelate`` (so the covariance is defined for any input, incl. one structure); with
            ``decorrelate=False`` it means cohort-relative (the ``table`` itself).
        features: the ``k`` descriptors. Default :data:`Q_FEATURES_GEOM` (``k=4``, geometry only) — the
            validated default. Adding the ``pp_combo`` energy term (``k=5``, :data:`Q_FEATURES`)
            *degrades* ranking on generated poses, because that term inverts on forced ones (ledger C42).
        method: per-descriptor standardization, ``"z"`` (default) or ``"rank"`` — see :func:`zscore`.
        decorrelate: whiten by the native covariance and project onto :math:`\mathbf 1` (default); else
            the equal-weight mean.
    """
    ref = native_reference() if (decorrelate and reference is None) else reference
    Z = np.vstack([zscore(_derive(table, f), None if ref is None else _derive(ref, f), method=method)
                   for f in features])                                   # k x n, standardized to ref
    if not decorrelate:
        return np.nanmean(Z, axis=0)
    Zref = np.vstack([zscore(_derive(ref, f), None, method=method) for f in features])   # k x n_ref
    C = np.atleast_2d(np.cov(Zref))                                      # native descriptor correlation
    w = np.linalg.pinv(C) @ np.ones(len(features))                      # C^{-1} 1: the decorrelated weights
    return np.nansum(w[:, None] * Z, axis=0)


def q_iptm(table, iptm, reference=None, features=Q_FEATURES_GEOM, decorrelate=True) -> np.ndarray:
    """Fit-free synergy score ``z(ipTM) + z(Q)`` — the interface-quality score composed with the
    generator's own confidence.

    ``Q`` (interface geometry) and the AlphaFold/TCRmodel2 ipTM are near-orthogonal (they fail in
    different pose regimes, benchmark ledger C26/C35), so their standardized sum out-ranks either alone:
    macro ROC 0.83 / PR 0.83 on TCRvdb vs ipTM 0.79, and on well-modelled epitopes it beats raw-AF ipTM
    on both metrics (ledger C42). Both terms are standardized over the same candidate set, so pass an
    ``iptm`` vector aligned row-for-row with ``table``. Use ``features=Q_FEATURES_GEOM`` for the
    geometry-only ``Q_geom`` variant that is robust to the forced-pose energy inversion.

    Args:
        table: the ``tcren recognize --full`` table (dict / pandas / polars).
        iptm: per-structure ipTM, aligned to ``table`` rows. Structures whose ipTM is missing (``NaN``)
            fall back to ``z(Q)`` alone, so the score always ranks; an all-missing ``iptm`` returns
            plain ``z(Q)`` — i.e. rank by the model geometry when no generator confidence is available.
        reference: optional cohort to standardize against (see :func:`zscore`).
        features: descriptors for ``Q``; defaults to the five :data:`Q_FEATURES`.
    """
    zq = zscore(q_score(table, reference, features, decorrelate=decorrelate))
    out = zscore(np.asarray(iptm, float)) + zq
    missing = ~np.isfinite(out)                       # ipTM absent for a structure -> rank by Q alone
    out[missing] = zq[missing]
    return out


def f_score(table, reference=None, terms=F_TERMS) -> np.ndarray:
    """Binder-oriented TCRen contact energy ``F = z(-(F_tcr_pep + F_tcr_mhc))`` — the chemistry channel.

    The standardized, sign-flipped sum of the :data:`F_TERMS` contact energies, so **higher = more
    binder-like** and it is on the same z-scale as :func:`q_score`. Unlike ``Q`` (interface geometry),
    ``F`` reads the actual contact chemistry — and unlike ``Q`` it is **pose-conditional**: it works on
    well-modelled poses and *inverts* on forced ones (benchmark ledger C27/C42). Do not use it
    unconditioned on pose quality; see :data:`F_TERMS` and :func:`q_f`.

    Cohort-relative (standardized over the ranked set); pass ``reference`` to standardize against another
    cohort (see :func:`zscore`).
    """
    e = sum(_col(table, t) for t in terms)
    ref = None if reference is None else sum(_col(reference, t) for t in terms)
    return zscore(-e, None if ref is None else -ref)


def q_f(table, reference=None, sign=1.0, features=Q_FEATURES_GEOM, terms=F_TERMS,
        decorrelate=True) -> np.ndarray:
    """Pure-tcren combiner ``z(Q_geom) + sign * z(F)`` — geometry plus contact energy, no deep learning.

    With ``sign=+1`` this is ``z(Q)+z(F)``; with ``sign=-1`` it is ``z(Q)-z(F)``. On **clean
    (template-covered) poses** ``z(Q)+z(F)`` beats raw-AF ipTM on both ROC and PR with no DL term
    (benchmark ledger C42: macro 0.759 ROC / 0.725 PR vs ipTM 0.692 / 0.693). On **forced poses** the
    energy inverts, so ``z(Q)-z(F)`` is the one that ranks (C27: on the forced GLCTLVAML pose
    ``z(Q)-z(F)``=0.71 vs ``z(Q)+z(F)``=0.52). Pick the sign from pose quality — grade it with
    :func:`strain_z` — or prefer :func:`q_iptm` (``z(ipTM)+z(Q)``), the geometry channel that is robust
    to the inversion without needing the energy at all.

    Args:
        table: the ``tcren recognize --full --scores`` table (dict / pandas / polars).
        reference: optional cohort to standardize against (see :func:`zscore`).
        sign: ``+1`` for ``z(Q)+z(F)`` (clean poses), ``-1`` for ``z(Q)-z(F)`` (forced poses).
        features: ``Q`` descriptors; defaults to the geometry-only :data:`Q_FEATURES_GEOM`.
        terms: energy terms for ``F``; defaults to :data:`F_TERMS`.
        decorrelate: passed to :func:`q_score`. ``False`` recovers the legacy equal-weight ``Q``.
    """
    return (zscore(q_score(table, reference, features, decorrelate=decorrelate))
            + sign * f_score(table, reference, terms))


def q_f_iptm(table, iptm, threshold=0.5, reference=None, features=Q_FEATURES_GEOM,
             terms=F_TERMS) -> np.ndarray:
    """Pose-adaptive ``z(Q) + s·z(F)`` where the F sign ``s`` is chosen per structure from ipTM.

    Automates the forced-pose inversion: a **confident** pose (``ipTM >= threshold``) keeps ``+z(F)``
    because the contact energy is trustworthy there; a **forced** pose (``ipTM < threshold``) flips to
    ``-z(F)`` because the energy inverts on forced poses (benchmark ledger C27/C42). A structure with no
    ipTM (``NaN``) keeps ``+z(F)`` — nothing marks it as forced. See :func:`f_invert_by_iptm` for the
    boolean flag alone.

    ipTM is a *pose-confidence* proxy, not a calibrated forced-pose detector — grading forced-ness with
    :func:`strain_z` is the principled alternative (C27), and :func:`q_iptm` (``z(ipTM)+z(Q)``) sidesteps
    the energy entirely. Provided because it is the single-call pose-adaptive combiner.

    Args:
        table: the ``tcren recognize --full --scores`` table (dict / pandas / polars).
        iptm: per-structure ipTM aligned to ``table`` rows.
        threshold: ipTM below which a pose is treated as forced and ``F`` is inverted (default 0.5).
        reference / features / terms: as in :func:`q_f`.
    """
    sign = np.where(f_invert_by_iptm(iptm, threshold), -1.0, 1.0)
    return q_f(table, reference, sign=sign, features=features, terms=terms)


def f_invert_by_iptm(iptm, threshold=0.5) -> np.ndarray:
    """Boolean per-structure flag: invert ``F`` where ``ipTM < threshold`` (a forced pose). ``NaN`` ipTM
    is not inverted. This is the ``F_invert`` column :func:`q_f_iptm` acts on."""
    iptm = np.asarray(iptm, float)
    return np.isfinite(iptm) & (iptm < threshold)


def phi_bind(table, reference=None) -> np.ndarray:
    """Deprecated screening score ``Phi_bind = Q + 0.5 * [z(-pitch) + z(-F_tcr_mhc)]``.

    .. deprecated::
       Use :func:`q_score`. Both terms this adds to ``Q`` *lower* ranking accuracy — on TCRvdb
       macro ROC falls from Q's 0.795 to 0.653, and ``z(-pitch)`` alone is below chance (0.43)
       (benchmark ledger C19b). The ``pitch`` axis also carries AlphaFold-confidence leakage
       (ledger C19). It is retained only to reproduce older figures; do not use it for new work.
    """
    warnings.warn("phi_bind is deprecated and degrades ranking vs q_score (benchmark ledger C19b); "
                  "use q_score", DeprecationWarning, stacklevel=2)
    ref_pitch = None if reference is None else -_col(reference, "pitch")
    ref_tm = None if reference is None else -_col(reference, "F_tcr_mhc")
    return (q_score(table, reference)
            + 0.5 * (zscore(-_col(table, "pitch"), ref_pitch)
                     + zscore(-_col(table, "F_tcr_mhc"), ref_tm)))


def coupling(q, energy) -> float:
    r"""Interface–energy coupling :math:`r=\mathrm{corr}(Q,\,\Delta\Phi)` over a cohort — the
    label-free forced-pose diagnostic.

    In a genuine complex the two channels are physically tied: a larger, better-packed interface
    holds more contacts, so favourable contact energy and good interface geometry rise together and
    :math:`r>0`. A structure generator that manufactures a pose optimises contacts *without* the
    interface, breaking the tie — the two channels decouple or run opposite, and :math:`r<0`.

    So the sign and size of :math:`r` say how far the energy of this cohort can be trusted, using
    **no labels and no reference set**. It is the weight :func:`q_coupled` admits the energy with.

    Args:
        q: interface-quality scores (e.g. :func:`q_score` output) for the cohort.
        energy: binder-oriented referenced contact energy for the same rows.

    Returns:
        Pearson :math:`r` over the rows where both are finite; ``0.0`` if fewer than three remain
        (an uninformative cohort contributes no energy evidence rather than a spurious weight).
    """
    q, energy = np.asarray(q, float), np.asarray(energy, float)
    ok = np.isfinite(q) & np.isfinite(energy)
    if ok.sum() < 3 or np.std(q[ok]) < 1e-12 or np.std(energy[ok]) < 1e-12:
        return 0.0
    return float(np.corrcoef(q[ok], energy[ok])[0, 1])


def q_coupled(q, energy) -> np.ndarray:
    r"""Parameter-free binder score: interface geometry **and** coupling-weighted contact energy.

    .. math::
       S(x) \;=\; \tfrac14\Big[1+\operatorname{erf}\tfrac{z(Q(x))}{\sqrt2}\Big]
                  \Big[1+\operatorname{erf}\tfrac{r\,z(\Delta\Phi(x))}{\sqrt2}\Big],
       \qquad r=\operatorname{corr}(Q,\Delta\Phi)

    Each bracket is :math:`2\times` a Gaussian tail probability, written with ``erf`` rather than the
    normal-CDF symbol so nothing collides with the potential :math:`\Phi`. Three biophysical
    statements, no free parameter — :math:`z` is standardization and :math:`r` is measured, not
    chosen.

    1. **Binding needs both.** A complex forms only if there *is* an interface and the residues in
       it are favourable. Each factor is the one-class probability that the candidate is native-like
       on that channel, and the product is the conjunction of two pieces of evidence — the smooth
       AND, with no threshold and no softness constant.
    2. **The energy is admitted in proportion to its coupling.** Under joint normality
       :math:`\mathbb E[z(Q)\mid z(\Delta\Phi)] = r\,z(\Delta\Phi)`, so :math:`r\,z(\Delta\Phi)` is
       exactly the part of the energy that is evidence about interface nativeness. Nothing is
       discarded and nothing is over-trusted.
    3. **A forced pose disarms itself.** :math:`r<0` on a fabricated cohort (:func:`coupling`), which
       flips the energy's sign automatically; :math:`r\approx0` shrinks the factor to
       :math:`\Phi(0)=\tfrac12`, a constant, leaving the geometry alone. The failure mode that makes
       raw :math:`\Phi` inverting and dangerous (ledger C27) is handled by the same :math:`r` that
       measures it.

    On TCRvdb this reaches macro ROC 0.799 / PR 0.817 / precision-at-10 %-recall 0.949, ahead of
    every TCRmodel2 confidence (best 0.795 / 0.800 / 0.916) with no generative term, and it is
    balanced across the two epitopes rather than trading one for the other.

    Args:
        q: interface-quality scores for the cohort, e.g. :func:`q_score` output.
        energy: binder-oriented referenced contact energy for the same rows — for receptor ranking
            use the **TCR**-referenced :math:`\Delta_{\mathrm{TCR}}\Phi` (the peptide is fixed there,
            so the peptide reference carries no signal); for peptide ranking use
            :math:`\Delta_{\mathrm{pep}}\Phi` (:func:`tcren.reference_delta`).

    Returns:
        Scores in :math:`(0,1)`; higher is more binder-like. Cohort-relative — rank within the set
        you scored, do not compare across cohorts.
    """
    from scipy.special import erf

    zq, ze = zscore(q), zscore(energy)
    g = lambda v: 0.5 * (1.0 + erf(v / np.sqrt(2.0)))          # noqa: E731 - Gaussian tail prob
    return g(zq) * g(coupling(zq, ze) * ze)


def strain_z(table, reference=None) -> np.ndarray:
    """Crystal-calibrated interface strain; higher = more forced. The recommended forced-pose score.

    Directional mean-z of :data:`STRAIN_TERMS` with fixed physical signs. Pass the crystal cohort
    as ``reference`` to reproduce the provenance gradient (crystal +0.02 < generated-real +0.40 <
    generated-decoy +0.81); without it the score is only relative within the input set.

    Prefer this over :func:`tcren.recognition.forced_pose_score` (``p_forced``): it is unfitted —
    no logistic, no coefficients, just signed standardization — so it carries no training set and
    is fully reproducible, whereas ``FORCED_POSE_MODEL``'s coefficients were frozen from a training
    set that no longer exists (benchmark ledger C23). It also grades forced-ness continuously, which
    is what pairs with :func:`q_score` to catch the forced poses where the contact energy inverts
    (ledger C27).
    """
    z = [sign * zscore(_derive(table, f),
                       None if reference is None else _derive(reference, f))
         for f, sign in STRAIN_TERMS]
    return np.nanmean(np.vstack(z), axis=0)
