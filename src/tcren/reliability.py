r"""Single-structure reliability: ``S_free``, its calibration, and the generator diagnostic.

``P_native`` refits a latent-class model per call and **raises when a cohort has fewer rows than
features**, so it is undefined for one structure and its published numbers depend on which rows the
fit was anchored on. ``S_free`` has neither property:

.. math::  S_{\mathrm{free}} \;=\; \frac{Q}{\sigma_Q} \;+\; \frac{T}{\sigma_T}
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
#: The partition-function-referenced energy S_free spends. See the module docstring.
PI_FROZEN = "neg_energy"


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

    Four keys. ``blocks`` carries the native mean and spread of each block, keyed by ``Q``, ``T``
    and the :math:`\\Pi` column named in ``pi_frozen``; those spreads are the divisors in
    :func:`s_free`. ``calibration`` carries the frozen Platt links :func:`p_binder` selects with
    ``link=``, and ``af_bands`` the confidence-band tables :func:`af_band` reads. Every one was
    fitted out of fold on the benchmark and is shipped rather than refitted, so a score computed
    today means what it meant when the paper was written.

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


def s_free(table, reference=None, energy=None) -> np.ndarray:
    """``S_free``, the recommended single-structure binder score. Higher = more native-like.

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


def p_binder(score, link: str = "binder_bm|S_nat") -> np.ndarray:
    """Map a score onto a probability through a frozen Platt link.

    The links were fitted OUT OF FOLD on the benchmarks — leave-one-epitope-out on the 22-cohort
    VDJdb panel, within-epitope 5-fold on TCRvdb — and the coefficients are the fold means. Names
    are ``<benchmark>|<score>``; :func:`available_links` lists them.

    A probability is a stronger claim than a rank, so read the expected calibration error beside it:
    the composed score reaches ECE 0.020 on the panel where ipTM alone reads 0.065.
    """
    cal = moments()["calibration"]
    if link not in cal:
        raise KeyError(f"no frozen link {link!r}; have {sorted(cal)}")
    c = cal[link]
    z = c["slope"] * np.asarray(score, float) + c["intercept"]
    return 1.0 / (1.0 + np.exp(-z))


def available_links() -> list[str]:
    """The frozen calibration links, ``<benchmark>|<score>``."""
    return sorted(moments()["calibration"])


def af_band(iptm, reference: str = "binder_bm|ipTM") -> list[dict]:
    """Look each confidence up in the frozen band table: how often is a model this confident wrong?

    The bands are deciles of the benchmark's own confidence distribution, never scanned for an
    effect. Each entry carries ``p_nonbinder`` with a Wilson interval and ``s_free_roc_in_band`` —
    what ``S_free`` still separates INSIDE that band, which is the actionable half: on the balanced
    VDJdb panel the top ipTM decile is 26% non-binders and is also where ``S_free`` reads highest.

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
        table: a ``tcren features`` table, as for :func:`s_free`.
        reference: overrides :func:`reliability_reference`.
        energy: :math:`\\Pi` per row, from :func:`tcren.potts.score_sites`' ``neg_energy``.
            Required — with no energy term there is nothing to invert, and the flag is NaN.

    Returns:
        One value per row: the energy block's native-sd score minus the mean of the geometry and
        topology blocks'. Large positive means the energy is vouching for a structure the shape
        does not, which is the pattern to distrust. It is a diagnostic to rank and inspect by, not
        a calibrated probability; ``p_binder`` is the calibrated read-out.
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
