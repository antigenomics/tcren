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

They are **cohort-relative**: each standardizes a feature over *the set being ranked*, so they are
defined for a candidate set, not for one structure. Score a whole batch together
(``tcren recognize`` over a directory), never one structure at a time. The division of labour is
scores in ``tcren``, evaluation (ROC/PR/CI) downstream.

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

import numpy as np

__all__ = ["zscore", "q_score", "phi_bind", "strain_z", "Q_FEATURES", "Q_FEATURES_CORE",
           "STRAIN_TERMS"]

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

#: Crystal-calibrated interface-strain terms with their physical signs. A forced pose reaches
#: further from the peptide with a thinner, less balanced interface.
STRAIN_TERMS = (("cdr3b_topep", +1.0), ("cdr3b_reach", +1.0),
                ("extent_per_ct", +1.0), ("chain_balance", -1.0))


def _col(table, name):
    """Fetch a column from a dict / pandas / polars frame as a float array."""
    if hasattr(table, "columns") and not isinstance(table, dict):
        if name not in table.columns:
            raise KeyError(f"column {name!r} not in table; run `tcren recognize --full`")
        col = table[name]
        return np.asarray(col.to_numpy() if hasattr(col, "to_numpy") else col, float)
    if name not in table:
        raise KeyError(f"column {name!r} not in table; run `tcren recognize --full`")
    return np.asarray(table[name], float)


def _derive(table, name):
    """Columns that ``recognize`` does not emit directly but are one division away."""
    if name == "extent_per_ct":  # interface thinness
        return _col(table, "extent") / np.maximum(_col(table, "n_contacts_tp"), 1.0)
    if name == "pp_combo":       # z(sum J CDR1/2) - z(sum J CDR3alpha)
        return zscore(_col(table, "e_cdr12")) - zscore(_col(table, "e_cdr3a"))
    return _col(table, name)


def zscore(x, reference=None) -> np.ndarray:
    """NaN-aware standardization. ``reference`` calibrates against another cohort.

    Passing ``reference`` is what makes :func:`strain_z` *crystal-calibrated*: the mean and sd come
    from the crystallographic ensemble, so the score reads ~0 on crystals by construction and grows
    as a pose departs from the natural manifold. Without it, a cohort of uniformly forced poses
    would standardize to zero mean and the shift would be invisible.
    """
    x = np.asarray(x, float)
    ref = x if reference is None else np.asarray(reference, float)
    mu = np.nanmean(ref)
    sd = np.nanstd(ref)
    # A constant column does not give sd == 0 exactly: np.nanstd(np.full(20, 3.7)) is 4.4e-16.
    # Testing `sd > 0` therefore divides by float residue and amplifies noise by ~1e16. Scale the
    # tolerance to the data so a genuinely constant (or degenerate) descriptor contributes nothing.
    if not np.isfinite(sd) or sd <= 1e-12 * max(1.0, abs(float(mu))):
        return np.zeros_like(x)
    return (x - mu) / sd


def q_score(table, reference=None, features=Q_FEATURES) -> np.ndarray:
    """Equal-weight interface-quality score ``Q = mean_k z(d_k)`` — the recommended binder score.

    No fitting and no labels: the descriptors enter with equal weight and a label-free
    standardization over the candidate set. Reproduces the shipped in-sample ``p_bind`` logistic at
    r ~ 0.92 while carrying no training set, and — unlike ``p_bind`` — generalises across cohorts
    (benchmark ledger C25). With ipTM, ``z(ipTM) + z(q_score(...))`` is the fit-free synergy score
    (macro ROC 0.83 on TCRvdb vs ipTM 0.79).

    Args:
        table: the ``tcren recognize --full`` table (dict / pandas / polars).
        reference: optional cohort to standardize against (see :func:`zscore`).
        features: which descriptors to average. Defaults to the five :data:`Q_FEATURES`; pass
            :data:`Q_FEATURES_CORE` for the simpler four-term score that is marginally better.
    """
    z = [zscore(_derive(table, f), None if reference is None else _derive(reference, f))
         for f in features]
    return np.nanmean(np.vstack(z), axis=0)


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
