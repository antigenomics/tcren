r"""Single-structure reliability: ``S`` and the generator diagnostic.

The cohort-fitted posterior this module replaced refitted a latent-class model per call and raised
when a cohort had fewer rows than features, so it was undefined for one structure and its numbers
depended on which rows the fit was anchored on. It was discarded in 2.26.0. ``S`` has neither
property:

.. math::  S \;=\; \frac{Q}{\sigma_Q} \;+\; \frac{T}{\sigma_T}
           \;+\; \frac{\Pi - \mu_\Pi}{\sigma_\Pi}

Three blocks, each a fit-free directional score :math:`z(x)^\top C^{-1} s` over the **Native2026
crystals**, divided by that block's native spread. No cohort, no EM, no anchors, no label anywhere,
and every term is defined for a single row.

**The outer transform is one divide, not a z.** A block score's native mean is 0 by construction
(:math:`z` is centred on the reference), but its variance is :math:`s^\top C^{-1} s`, which is not 1
— native spreads run 1.43 (``Q``), 1.61 (``T``) and 14.13 (:math:`\Pi`), so without the divide the
energy would carry ten times the weight of the geometry. Equal weight *in native-sd units* is the
claim, and the division is what makes it true.

:math:`\Pi` is the interface energy read against the partition function rather than against a
poly-alanine reference — ``neg_energy`` from :mod:`tcren.potts`, which is
:math:`-E(\sigma^{\mathrm{obs}}) = \log Z + \mathcal L`. It is the least redundant with ``Q`` of the
five ways of spending that decomposition (native Pearson +0.33, against +0.75 for the contact
count), which is why it is the frozen choice.

Provenance of every frozen constant: ``data/reliability_moments.json``, written by the benchmark's
``bench/eda/blockcov.py`` and ``bench/eda/calibrate.py``.
"""
from __future__ import annotations

import json
from functools import lru_cache

import numpy as np

from .cohort import Q_FEATURES_GEOM, q_score

#: The topology block: the SHAPE of the contact set, free of its size.
T_FEATURES_TOPO = ("D2_pep24", "fp_b0_frac_r7", "H_cell", "L_canon", "ab_imb")
#: Its orientation. Every term rises towards a native interface except the footprint's
#: connected-component fraction at 7 A, which falls.
T_SIGNS = (1.0, -1.0, 1.0, 1.0, 1.0)
#: The partition-function-referenced energy S spends. See the module docstring.
PI_FROZEN = "neg_energy"

#: Columns only the ``potts`` family emits. A feature table carrying ``n_contacts`` alongside none
#: of these was written by tcren <= 2.19.0 without ``-i potts``, so its ``n_contacts`` is the
#: footprint's CDR-loop tally (now ``n_loop_contacts``) rather than the engaged-pair count these
#: moments were estimated on. See :func:`_check_potts_contacts`.
_POTTS_MARKERS = (PI_FROZEN, "log_z", "log_lik", "psi", "n_sites", "mu_star")


def _column_names(table) -> set[str]:
    """The column names of a dict / pandas / polars table, the way :func:`tcren.cohort._col` reads them."""
    if hasattr(table, "columns") and not isinstance(table, dict):
        return set(table.columns)
    return set(table)


def _check_potts_contacts(table) -> None:
    """Raise unless a table's ``n_contacts`` came from :mod:`tcren.potts`.

A ``n_contacts`` read from the footprint rather than from the Potts model differs by a factor of
    two to four on real structures (1ao7: 66 footprint contacts against 29 Potts ones), so any term
    standardized against the Potts population would shift by several native sd with no error and no
    NaN. A table with no
    ``n_contacts`` column at all is fine — the caller joined the count from somewhere this cannot
    see, which is the documented path.
    """
    cols = _column_names(table)
    if "n_contacts" in cols and not cols & set(_POTTS_MARKERS):
        raise ValueError(
            "this table's 'n_contacts' did not come from the potts family: it carries none of "
            f"{list(_POTTS_MARKERS)}, so the column is tcren <= 2.19.0's footprint CDR-loop tally "
            "(renamed 'n_loop_contacts'), which the frozen correction is not standardized on. "
            "Rebuild with `tcren features -i placement,interface,topology,potts`, or drop the "
            "column and pass contacts= from `tcren potts score`.")


@lru_cache(maxsize=1)
def reliability_reference() -> dict:
    """The Q and T descriptors plus the Potts energies over the Native2026 crystals.

    Same role as :func:`tcren.cohort.native_reference`, extended to the topology block and the
    partition-function terms so a single user structure can be standardized against the crystal
    manifold for all three blocks at once.
    """
    import csv
    from importlib import resources

    with (resources.files("tcren.data") / "reliability_reference.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for k in rows[0]:
        if k == "complex.id":
            continue
        out[k] = np.array([float(r[k]) if r[k] not in ("", "null") else np.nan for r in rows])
    return out


@lru_cache(maxsize=1)
def moments() -> dict:
    """Every frozen constant this module reads, from ``data/reliability_moments.json``.

Four keys, and **none of them is a fit against a binding label.** ``blocks`` carries the native
    mean and spread of each block, keyed by ``Q``, ``T`` and the :math:`\\Pi` column named in
    ``pi_frozen``; those spreads are the divisors in :func:`s_score` and are measured on the
    Native2026 crystals. ``af_bands`` carries the confidence-band tables :func:`af_band` reads,
    which are quantile bins of the generator's confidence with the observed non-binder fraction and
    its Wilson interval in each. ``phi`` carries the per-interface energy spreads.

    The out-of-fold-fitted sections that used to sit here — the Platt links and the confidence
    correction — were removed in 2.28.0; see the manuscript repository's ``LEGACY.md``.

    Returns:
        The parsed JSON. Cached, so callers may treat it as read-only.
    """
    from importlib import resources

    return json.loads((resources.files("tcren.data") / "reliability_moments.json").read_text())


def t_score(table, reference=None) -> np.ndarray:
    """Topology block ``T`` — the same construction as :func:`tcren.cohort.q_score`, other terms.

    ``T`` reads the SHAPE of the footprint rather than its size, which is why it survives where the
    geometry block does not: on the balanced VDJdb panel it loses 0.06 ROC-AUC when the epitope has
    no solved complex to template on, against ``Q``'s 0.24.
    """
    return q_score(table, reference=reliability_reference() if reference is None else reference,
                   features=T_FEATURES_TOPO, signs=T_SIGNS)


def s_score(table, reference=None, energy=None) -> np.ndarray:
    """``S``, the recommended single-structure binder score. Higher = more native-like.

    Args:
        table: a ``tcren features`` table (dict / pandas / polars) carrying the four geometry and
            five topology descriptors.
        reference: overrides :func:`reliability_reference`.
        energy: :math:`\\Pi` per row, from :func:`tcren.potts.bound_unbound`'s ``neg_energy``. When
            ``None`` the energy term is dropped and the score is ``Q/sd + T/sd`` — still defined,
            and reported as such rather than silently imputed.

    Returns:
        One value per row. NaN only where a whole block is unavailable.
    """
    ref = reliability_reference() if reference is None else reference
    m = moments()["blocks"]
    q = q_score(table, reference=ref, features=Q_FEATURES_GEOM) / m["Q"]["sd"]
    t = t_score(table, reference=ref) / m["T"]["sd"]
    if energy is None:
        return q + t
    e = np.asarray(energy, float)
    return q + t + (e - m[PI_FROZEN]["mean"]) / m[PI_FROZEN]["sd"]
def af_band(iptm, reference: str = "binder_bm|ipTM") -> list[dict]:
    """Look each confidence up in the frozen band table: how often is a model this confident wrong?

    The bands are deciles of the benchmark's own confidence distribution, never scanned for an
    effect. Each entry carries ``p_nonbinder`` with a Wilson interval and ``s_roc_in_band`` —
    what ``S`` still separates INSIDE that band, which is the actionable half: on the balanced
    VDJdb panel the top ipTM decile is 26% non-binders and is also where ``S`` reads highest.

    Values outside the reference range clamp to the end bands rather than extrapolating.
    """
    tab = moments()["af_bands"]
    if reference not in tab:
        raise KeyError(f"no band table {reference!r}; have {sorted(tab)}")
    b = tab[reference]
    out = []
    for v in np.atleast_1d(np.asarray(iptm, float)):
        if not np.isfinite(v):
            out.append({})
            continue
        k = next((i for i, e in enumerate(b) if v < e["hi"]), len(b) - 1)
        out.append(dict(b[k], band=k))
    return out


def available_bands() -> list[str]:
    """The frozen confidence-band tables, ``<benchmark>|<confidence>``."""
    return sorted(moments()["af_bands"])


def inversion_flag(table, reference=None, energy=None) -> np.ndarray:
    """Is this model's recognition energy running *backwards*? The forced-pose detector.

    A co-folding model that has been pushed into a confident but wrong pose does not produce a
    random interface. To seat the chains it selects residue pairs it believes are favourable, so
    the recognition energy comes out **good** — better, often, than a genuine complex of the same
    epitope. The energy therefore inverts under forcing rather than degrading, and the inversion is
    a signal with its sign flipped, not an absence of signal.

    Measured on a 24-cohort forced-pose panel (1,707 structures), the TCR:peptide energy reads
    macro ROC-AUC 0.4952 against the forced poses and is below 0.5 in **15 of the 24 cohorts**,
    where the generator's own ipTM reads 0.6093. That is what this flag reads off one structure:
    an energy far better than the crystal manifold's while the geometry and topology blocks are
    not is the forced-pose signature, because a generator can fake favourable contacts far more
    easily than it can fake a well-formed footprint.

    Args:
        table: a ``tcren features`` table, as for :func:`s_score`.
        reference: overrides :func:`reliability_reference`.
        energy: :math:`\\Pi` per row, from :func:`tcren.potts.score_sites`' ``neg_energy``.
            Required — with no energy term there is nothing to invert, and the flag is NaN.

    Returns:
        One value per row: the energy block's native-sd score minus the mean of the geometry and
        topology blocks'. Large positive means the energy is vouching for a structure the shape
        does not, which is the pattern to distrust. It is a diagnostic to rank and inspect by, and
        the package ships no probability to threshold on.
    """
    ref = reliability_reference() if reference is None else reference
    m = moments()["blocks"]
    q = q_score(table, reference=ref, features=Q_FEATURES_GEOM) / m["Q"]["sd"]
    t = t_score(table, reference=ref) / m["T"]["sd"]
    if energy is None:
        return np.full(np.shape(q), np.nan)
    e = (np.asarray(energy, float) - m[PI_FROZEN]["mean"]) / m[PI_FROZEN]["sd"]
    return e - 0.5 * (q + t)


def screening_yield(score, budget: float = 0.1, prevalence: float | None = None) -> dict:
    """What a caller gets for testing the top ``budget`` fraction of a scored set.

    ROC-AUC answers "does the score order the set". It does not answer the question an
    experimenter asks, which is "I can test ten of these hundred models -- how many binders do I
    get?". This is that number.

    Args:
        score: the score for every candidate in the set, higher = more likely to bind.
        budget: fraction of the set that can be tested, in (0, 1].
        prevalence: the hit rate expected if the set were tested at random. When given, the
            returned ``expected_hits`` is what testing that slice blindly would yield -- the number
            the score has to beat to have been worth computing.

    Returns:
        ``n_tested``, the ``threshold`` score at the cut, and ``rank_cut``, the percentile it
        corresponds to; plus ``expected_hits`` when ``prevalence`` is given.

    Note:
        Enrichment -- hits over the random baseline -- is deliberately not returned. It needs the
        labels this function does not have, and returning it as NaN would read like a measurement.
        The benchmark computes it where the labels are, in ``bench/metrics.screening_yield``.
    """
    if not 0 < budget <= 1:
        raise ValueError(f"budget must be in (0, 1], got {budget}")
    s = np.asarray(score, float)
    ok = np.isfinite(s)
    n = int(ok.sum())
    if not n:
        return {"n_tested": 0, "threshold": float("nan"), "rank_cut": float("nan")}
    k = max(1, int(np.ceil(budget * n)))
    cut = float(np.sort(s[ok])[::-1][k - 1])
    out = {"n_tested": k, "threshold": cut, "rank_cut": 1.0 - k / n}
    if prevalence is not None:
        out["expected_hits"] = float(k * prevalence)
    return out




# ---------------------------------------------------------------------------------------------
# The artefact test: which directions belong to the generator rather than to the interface.
# ---------------------------------------------------------------------------------------------

#: A binder direction tighter than this (residual s.d., transformed reference units) is where the
#: crystal test finds the generator's own regularity rather than interface physics. Measured
#: 2026-09-03: real crystals break the sub-0.05 directions 4.55x their binder spread against
#: decoys' 3.18x, and that band scores 0.504 over 16 template-free cohorts -- a coin -- where the
#: directions above 0.15 read 0.606. Ledoit-Wolf shrinkage floors the shipped model's smallest
#: direction at 0.0797, so `tcren.score.pose_score` cannot read this band at all.
ARTEFACT_SD = 0.05


def artefact_directions(binder_coords, crystal_coords, *, tolerance: float = ARTEFACT_SD) -> dict:
    """Which directions of a binder manifold are the generator's regularity, not the interface's.

    **No binding label enters this.** The argument is that a direction real complexes are free to
    vary along, while modelled binders are pinned to it, is a property of the model that produced
    them. A physical constraint holds on a crystal at least as tightly as on a prediction of one;
    a generator's habit does not.

    Args:
        binder_coords: ``(n, p)`` transformed coordinates for modelled binders -- the population
            whose covariance defines the directions.
        crystal_coords: ``(m, p)`` the same coordinates for experimentally solved complexes.
        tolerance: the residual s.d. below which a direction counts as tight.

    Returns:
        ``{"index", "sd_binder", "sd_crystal", "ratio", "is_artefact"}``, one entry per direction,
        ordered stiffest first. ``is_artefact`` is a tight direction that crystals break harder
        than binders do.

    Reading it: a large ``ratio`` on a tight direction means the constraint is the generator's.
    Drop those coordinates, or use a covariance estimator that suppresses them -- shrinkage does it
    for free, which is why the shipped model does not need this call to be safe.
    """
    import numpy as np

    B = np.asarray(binder_coords, float)
    C = np.asarray(crystal_coords, float)
    if B.ndim != 2 or C.ndim != 2 or B.shape[1] != C.shape[1]:
        raise ValueError(f"coordinate shapes disagree: {B.shape} against {C.shape}")
    mu = B.mean(0)
    lam, U = np.linalg.eigh(np.cov(B - mu, rowvar=False))
    lam, U = np.maximum(lam, 1e-12), U
    sd_b = np.sqrt(lam)
    sd_c = ((C - mu) @ U).std(0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = sd_c / sd_b
    return {"index": np.arange(len(lam)), "sd_binder": sd_b, "sd_crystal": sd_c,
            "ratio": ratio, "is_artefact": (sd_b < tolerance) & (ratio > 1.0)}
