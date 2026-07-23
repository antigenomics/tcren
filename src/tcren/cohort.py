"""Cohort-relative recognition scores: ``Q``, ``Phi_bind``, ``S_strain``, ``kit_score``.

These are the manuscript's headline screening scores. They are **cohort-relative**: each
standardizes a feature over *the set being ranked*, so they are defined for a candidate set, not
for one structure. That is precisely why they used to live in analysis scratch — and why the
paper's headline numbers were not regenerable from ``tcren`` alone. They are here now; the
division of labour is scores in ``tcren``, evaluation (ROC/PR/CI) downstream.

All functions take the table ``tcren recognize --full`` emits (a mapping of column name to
sequence, a ``polars``/``pandas`` frame, or a dict of arrays) and return one value per row.

Sign convention: every term is oriented so that **higher = more binder-like** for
:func:`q_score`/:func:`phi_bind`, and **higher = more forced/strained** for :func:`strain_z`.
"""

from __future__ import annotations

import numpy as np

__all__ = ["zscore", "q_score", "phi_bind", "strain_z", "Q_FEATURES", "STRAIN_TERMS"]

#: The five interface-quality descriptors, equal-weighted in :func:`q_score`. Each is oriented
#: positive-is-better as given. ``pp_combo`` is the CDR1/2-vs-CDR3alpha potential contrast.
Q_FEATURES = ("burial", "n_pep_contacted", "chain_balance", "n_hbond", "pp_combo")

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


def q_score(table, reference=None) -> np.ndarray:
    """Equal-weight interface-quality score ``Q = (1/5) * sum z(d_k)``.

    No fitting and no labels: the five descriptors enter with equal weight and a label-free
    standardization over the candidate set. Reproduces the shipped in-sample binder logistic at
    r ~ 0.92, so nothing is trained on the benchmark it is evaluated on.
    """
    z = [zscore(_derive(table, f), None if reference is None else _derive(reference, f))
         for f in Q_FEATURES]
    return np.nanmean(np.vstack(z), axis=0)


def phi_bind(table, reference=None) -> np.ndarray:
    """Screening score ``Phi_bind = Q + 0.5 * [z(-pitch) + z(-F_tcr_mhc)]``.

    Adds an orthogonal docking-correctness axis to :func:`q_score`: binders dock canonically (low
    pitch) and make favourable TCR:MHC contact (low, i.e. more negative, energy).

    .. warning::
       ``pitch`` here is the **recomputed** geometric incident angle from
       :mod:`tcren.orient.docking`. Do **not** substitute a generator-cached ``pitch_angle``
       column: on validation it matched no clean geometric angle (best r ~ 0.42) while
       out-discriminating every clean docking feature, i.e. it carries generator-confidence
       leakage. Using it would silently void the predictor-independence of this score.
    """
    ref_pitch = None if reference is None else -_col(reference, "pitch")
    ref_tm = None if reference is None else -_col(reference, "F_tcr_mhc")
    return (q_score(table, reference)
            + 0.5 * (zscore(-_col(table, "pitch"), ref_pitch)
                     + zscore(-_col(table, "F_tcr_mhc"), ref_tm)))


def strain_z(table, reference=None) -> np.ndarray:
    """Crystal-calibrated interface strain; higher = more forced.

    Directional mean-z of :data:`STRAIN_TERMS` with fixed physical signs. Pass the crystal cohort
    as ``reference`` to reproduce the provenance gradient (crystal < generated-real <
    generated-decoy); without it the score is only relative within the input set.

    Unlike :func:`tcren.recognition.forced_pose_score` this is unfitted — no logistic, no
    coefficients, just signed standardization — so it carries no training set at all.
    """
    z = [sign * zscore(_derive(table, f),
                       None if reference is None else _derive(reference, f))
         for f, sign in STRAIN_TERMS]
    return np.nanmean(np.vstack(z), axis=0)
