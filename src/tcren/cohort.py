"""Cohort-relative recognition scores — the **recommended, fit-free** screening layer.

These carry no trained coefficients — no logistic, no fit, no training set — so they cannot leak,
cannot go stale, and there is nothing to re-derive. The fitted composites that stood beside them
were removed in 2.26.0; their coefficients were frozen against training sets that no longer exist,
which made them the one part of the package a reader could not reproduce.

* :func:`q_score` **generalises across cohorts**: a logistic trained on one cohort learns that
  cohort's epitope composition and does not transfer, whereas ``Q`` has nothing to transfer. With
  ipTM it reproduces the headline synergy fit-free.
* :func:`strain_z` grades pose forcedness (crystal < AF-real < AF-decoy) reproducibly.

They are **cohort-relative** by default: each standardizes a feature over *the set being ranked*.
For a candidate set, score the whole batch together (``tcren recognize`` over a directory). For a
**single structure**, or a small/heterogeneous user set where the batch is not a fair reference, pass
``reference=native_reference()`` (with ``features=Q_FEATURES_GEOM``): the descriptors are then
standardized against the shipped Native2026 crystal manifold, so ``Q`` is defined for one structure and
transfers across inputs. The descriptors are counts and bounded ratios (mildly non-normal), so
``method="rank"`` gives a robust, assumption-free percentile standardization; on the benchmarks it
agrees with the default ``z`` to ρ≈0.98. The division of labour is scores in ``tcren``, evaluation
(ROC/PR/CI) downstream.

All functions take the table ``tcren features`` emits (a mapping of column name to sequence, a
``polars``/``pandas`` frame, or a dict of arrays) and return one value per row.

**Where the line is drawn.** Every *score* the TCRen2 manuscript reports is computed here or in
:mod:`tcren.footprint`, :mod:`tcren.mechanics`, :mod:`tcren.ddg` and :mod:`tcren.potential` --- a
benchmark script that recomputes one of them by hand is a bug, not a shortcut. What stays outside
the library is *evaluation*: ROC/PR/AUC, bootstrap intervals, macro-averaging over cohorts, and any
protocol that consumes a binder label (leave-one-epitope-out anchoring, an in-sample GLM against a
generator's confidence). Those need the labels this library is built to do without, so they live in
the benchmark repo next to the data that carries them.

Sign convention: every term is oriented so that **higher = more binder-like** for
:func:`q_score`, and **higher = more forced/strained** for :func:`strain_z`.

.. note::
   The hand-written combination rules this module once exposed were removed in 2.12.0, and the
   fitted cohort posterior that replaced them was itself discarded in 2.26.0. Use
   :func:`q_score` for the single-structure interface-quality score, and
   :func:`tcren.reliability.s_free` for the composition.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

__all__ = ["zscore", "q_score", "phi_score",
           "q_coupled", "coupling", "strain_z", "native_reference", "Q_FEATURES", "Q_FEATURES_CORE",
           "Q_FEATURES_GEOM", "PHI_TERMS", "STRAIN_TERMS"]

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
#: term is used only conditioned on pose quality. Pass to :func:`q_score`.
Q_FEATURES_GEOM = ("burial", "n_pep_contacted", "chain_balance", "n_hbond")

#: The TCRen contact-energy terms summed into the binder-oriented :func:`phi_score`. ``Phi_tcr_pep`` is the
#: TCR:peptide TCRen energy, ``Phi_tcr_mhc`` the TCR:MHC energy; both are emitted by ``tcren recognize``.
#: They are raw energies (lower = tighter), so :func:`phi_score` negates the sum to make higher = more
#: binder-like. **This term is pose-conditional** — it reads real binding chemistry on well-modelled
#: (crystal-templated) poses and *inverts* on forced ones (benchmark ledger C27/C42): on the forced
#: GLCTLVAML TCRvdb pose ``-Phi_tcr_pep`` ranks binders at AUROC 0.36 (backwards), on the clean YLQPRTFLL
#: pose at 0.59. Use it only conditioned on pose quality — gate with :func:`strain_z`, or read
#: ``z(Q)-z(F)`` on forced poses and ``z(Q)+z(F)`` on clean ones.
PHI_TERMS = ("Phi_tcr_pep", "Phi_tcr_mhc")

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
        return zscore(_col(table, "Phi_cdr12")) - zscore(_col(table, "Phi_cdr3a"))
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
    column arrays (``burial, n_pep_contacted, chain_balance, n_hbond, Phi_cdr12, Phi_cdr3a``) usable as
    the ``reference`` argument. Provenance: ``tcren recognize --full`` over the Native2026 set.
    """
    import csv
    from importlib import resources
    path = resources.files("tcren.data") / "q_native_reference.csv"
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    return {c: np.array([float(r[c]) for r in rows]) for c in rows[0]}


def q_score(table, reference=None, features=Q_FEATURES_GEOM, method="z", decorrelate=True,
            signs=None) -> np.ndarray:
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
        signs: per-descriptor orientation replacing :math:`\mathbf 1`, for a block whose terms are not
            all "higher = more native" — the topology block's footprint fraction runs the other way.
            Length must match ``features``. ``None`` keeps :math:`\mathbf 1`.
    """
    ref = native_reference() if (decorrelate and reference is None) else reference
    Z = np.vstack([zscore(_derive(table, f), None if ref is None else _derive(ref, f), method=method)
                   for f in features])                                   # k x n, standardized to ref
    if not decorrelate:
        return np.nanmean(Z, axis=0)
    Zref = np.vstack([zscore(_derive(ref, f), None, method=method) for f in features])   # k x n_ref
    C = np.atleast_2d(np.cov(Zref))                                      # native descriptor correlation
    sgn = np.ones(len(features)) if signs is None else np.asarray(signs, float)
    if len(sgn) != len(features):
        raise ValueError(f"signs has {len(sgn)} entries for {len(features)} features")
    w = np.linalg.pinv(C) @ sgn                                          # C^{-1} s: the decorrelated weights
    return np.nansum(w[:, None] * Z, axis=0)



def phi_score(table, reference=None, terms=PHI_TERMS) -> np.ndarray:
    """Binder-oriented TCRen contact energy ``F = z(-(Phi_tcr_pep + Phi_tcr_mhc))`` — the chemistry channel.

    The standardized, sign-flipped sum of the :data:`PHI_TERMS` contact energies, so **higher = more
    binder-like** and it is on the same z-scale as :func:`q_score`. Unlike ``Q`` (interface geometry),
    ``F`` reads the actual contact chemistry — and unlike ``Q`` it is **pose-conditional**: it works on
    well-modelled poses and *inverts* on forced ones (benchmark ledger C27/C42). Do not use it
    unconditioned on pose quality; see :data:`PHI_TERMS`.

    Cohort-relative (standardized over the ranked set); pass ``reference`` to standardize against another
    cohort (see :func:`zscore`).
    """
    e = sum(_col(table, t) for t in terms)
    ref = None if reference is None else sum(_col(reference, t) for t in terms)
    return zscore(-e, None if ref is None else -ref)



def coupling(q, energy) -> float:
    r"""Interface–energy coupling :math:`r=\mathrm{corr}(Q,\,\Delta\Phi)` over a cohort — the
    label-free forced-pose diagnostic.

    .. deprecated:: 2.12
       No longer a component of any recommended score; keep it as a **diagnostic** you report, not
       a weight you apply. What it measures is real and worth knowing — on the heavily crystallised
       GLCTLVAML cohort it reads −0.2617 and the referenced energy ranks binders at AUROC 0.338
       [0.250, 0.433], entirely below chance, while on the sparsely templated YLQPRTFLL it reads
       +0.4784 and the same energy reads 0.776 [0.728, 0.820].

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



def q_coupled(q, energy, r=None) -> np.ndarray:
    r"""Parameter-free binder score: interface geometry **and** coupling-weighted contact energy.

    .. deprecated:: 2.12
       No longer a component of any recommended score. Nothing here changes: this function
       returns exactly what it always has, and the numbers it produces stand (TCRvdb macro ROC
       0.802 / PR 0.817). What changed is that the footprint-shape channel makes the gate
       unnecessary — the energy is one input among four rather than a term that has to be disarmed.

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
        r: coupling weight. ``None`` (default) measures it from the cohort with :func:`coupling`,
            which is what every published number uses. Pass a scalar or a per-row array to supply
            it from outside — a **predicted** :math:`\hat C` from a single structure is defined at
            ``n = 1``, where the cohort estimator is not (at ``n = 2`` it returns
            :math:`|\hat r| = 1` by construction, with the wrong sign in 43.4 % of GLCTLVAML draws).

    Returns:
        Scores in :math:`(0,1)`; higher is more binder-like. Cohort-relative — rank within the set
        you scored, do not compare across cohorts.
    """
    from scipy.special import erf

    zq, ze = zscore(q), zscore(energy)
    g = lambda v: 0.5 * (1.0 + erf(v / np.sqrt(2.0)))          # noqa: E731 - Gaussian tail prob
    w = coupling(zq, ze) if r is None else np.asarray(r, float)
    return g(zq) * g(w * ze)



def strain_z(table, reference=None) -> np.ndarray:
    """Crystal-calibrated interface strain; higher = more forced. The recommended forced-pose score.

    Directional mean-z of :data:`STRAIN_TERMS` with fixed physical signs. Pass the crystal cohort
    as ``reference`` to reproduce the provenance gradient (crystal +0.02 < generated-real +0.40 <
    generated-decoy +0.81); without it the score is only relative within the input set.

    Unfitted — no logistic, no coefficients, just signed standardization — so it carries no training
    set and is fully reproducible. It grades forced-ness continuously, which is what pairs with
    :func:`q_score` to catch the forced poses where the contact energy inverts.
    """
    z = [sign * zscore(_derive(table, f),
                       None if reference is None else _derive(reference, f))
         for f, sign in STRAIN_TERMS]
    return np.nanmean(np.vstack(z), axis=0)
