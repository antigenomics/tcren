"""Sampling contact maps, and the partition function.

With couplings on, ``Z`` no longer factorises. Two pieces do the work:

**Block Gibbs.** Colour the coupling graph (:func:`tcren.potts.colour`); sites of one colour are
conditionally independent, so a whole colour updates in one vectorised step and a full sweep costs
one pass over the colours rather than one pass over the sites.

**Annealed importance sampling** (Neal, *Stat Comput* 11:125–139, 2001), annealing *only* the
coupling term from ``beta = 0`` to ``1``. At ``beta = 0`` the model is the factorised one, whose
partition function is exact and closed form, ``log Z_0 = sum_a log(1 + exp(eta_a))`` — so the
reference is not an approximation but the uncoupled model itself, and ``log Z`` is estimated
unbiasedly in ``Z``.
"""

from __future__ import annotations

import numpy as np


def _logsumexp(x: np.ndarray) -> float:
    m = float(np.max(x))
    return m + float(np.log(np.sum(np.exp(x - m))))


def factorised_log_z(eta: np.ndarray) -> float:
    """``log Z`` of the uncoupled model — exact, in closed form."""
    return float(np.logaddexp(0.0, eta).sum())


def gibbs(eta: np.ndarray, A: np.ndarray, colours, rng, *, chains: int = 64,
          burn: int = 100, draws: int = 100, thin: int = 3, observer=None):
    """Block Gibbs at ``beta = 1``.

    Args:
        eta: One-body log-odds per site.
        A: Symmetric coupling matrix with zero diagonal.
        colours: Colour classes from :func:`tcren.potts.colour`.
        rng: A ``numpy.random.Generator``.
        chains: Independent chains run in parallel (vectorised, not looped).
        burn: Sweeps discarded before the first draw.
        draws: Draws kept per chain.
        thin: Sweeps between kept draws.
        observer: Optional callable invoked on each kept draw with the ``(chains, n)`` matrix of
            configurations. Use it to accumulate a statistic of whole configurations without
            materialising every draw — the matrix is reused between calls, so copy what you keep.

    Returns:
        ``(occupancy, totals)`` — the mean of ``sigma`` per site over all kept draws, and the
        contact total of every kept draw.
    """
    n = len(eta)
    sig = (rng.random((chains, n)) < 1.0 / (1.0 + np.exp(-eta))).astype(np.float64)
    acc, tot, keep = np.zeros(n), 0, []
    coupled = A.any()
    for it in range(burn + draws * thin):
        for cidx in colours:
            f = eta[cidx] + (sig @ A[:, cidx] if coupled else 0.0)
            sig[:, cidx] = (rng.random((chains, len(cidx)))
                            < 1.0 / (1.0 + np.exp(-f))).astype(np.float64)
        if it >= burn and (it - burn) % thin == 0:
            acc += sig.sum(0)
            tot += chains
            keep.append(sig.sum(1).copy())
            if observer is not None:
                observer(sig)
    return acc / max(tot, 1), (np.concatenate(keep) if keep else np.zeros(0))


def ais_log_z(eta: np.ndarray, A: np.ndarray, colours, rng, *,
              particles: int = 64, steps: int = 256) -> tuple[float, float]:
    """``log Z`` by annealed importance sampling. Returns ``(log_z, effective_sample_size)``.

    The ``beta = 0`` reference is the factorised model, sampled exactly, so no burn-in is needed
    and the estimator is unbiased in ``Z``. An effective sample size close to ``particles`` means
    the annealing schedule is long enough; a small one means it is not.
    """
    n = len(eta)
    log_z0 = factorised_log_z(eta)
    if not colours or not A.any():
        return log_z0, float(particles)
    sig = (rng.random((particles, n)) < 1.0 / (1.0 + np.exp(-eta))).astype(np.float64)
    betas = np.linspace(0.0, 1.0, steps + 1)
    w = np.zeros(particles)
    for m in range(1, steps + 1):
        w += (betas[m] - betas[m - 1]) * 0.5 * ((sig @ A) * sig).sum(1)
        for cidx in colours:
            f = eta[cidx] + betas[m] * (sig @ A[:, cidx])
            sig[:, cidx] = (rng.random((particles, len(cidx)))
                            < 1.0 / (1.0 + np.exp(-f))).astype(np.float64)
    lse = _logsumexp(w)
    return log_z0 + lse - float(np.log(particles)), float(np.exp(2 * lse - _logsumexp(2 * w)))


def exact_log_z(eta: np.ndarray, A: np.ndarray) -> float:
    """``log Z`` by brute force over all ``2^n`` configurations. For testing; ``n <= 22``."""
    n = len(eta)
    if n > 22:
        raise ValueError(f"exact_log_z enumerates 2^n states; n = {n} is too many")
    s = ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1).astype(np.float64)
    return _logsumexp(s @ eta + 0.5 * ((s @ A) * s).sum(1))


def energy(sigma: np.ndarray, eta: np.ndarray, A: np.ndarray) -> float:
    """``E(sigma) = -eta . sigma - 1/2 sigma' A sigma``. Lower is more favourable."""
    return float(-(sigma @ eta) - 0.5 * (sigma @ A @ sigma))


def tilt_mean(totals: np.ndarray, mu: float) -> float:
    """``<N>`` under the tilted model ``E_mu = E - mu N``, by reweighting draws taken at ``mu = 0``.

    Exact in the tilted family, because every tilt considered here depends on ``sigma`` only
    through ``N(sigma)``: the importance weight is ``exp(mu N)`` and needs no new sampling.
    """
    w = mu * totals
    p = np.exp(w - w.max())
    return float((p * totals).sum() / p.sum())


def mu_star(totals: np.ndarray, n_obs: float, *, span: float = 40.0,
            tol: float = 1e-8) -> float:
    """The chemical potential at which the model's mean contact count matches ``n_obs``.

    ``<N>_mu`` is non-decreasing in ``mu`` (its derivative is the tilted variance), so a bisection
    is exact. Returns ``nan`` when ``n_obs`` lies outside the sampled support, where reweighting
    would be extrapolating rather than reweighting.
    """
    if not len(totals) or not (totals.min() < n_obs < totals.max()):
        return float("nan")
    lo, hi = -span, span
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if tilt_mean(totals, mid) < n_obs:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def count_free_energy(totals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(N, F(N))`` with ``F(N) = -log p(N)`` from the sampled contact-count histogram.

    The Legendre partner of ``log Z(mu)``, and the same object the protein-folding literature
    draws along the native-contact coordinate. Only the sampled range is returned: the histogram
    cannot see the unbound basin of a docked pose, which is what :func:`delta_f_empty` is for.
    """
    if not len(totals):
        return np.zeros(0), np.zeros(0)
    n = np.arange(int(totals.min()), int(totals.max()) + 1)
    c = np.array([(totals == k).sum() for k in n], float)
    with np.errstate(divide="ignore"):
        return n, -np.log(c / c.sum())


def delta_f_empty(log_z: float) -> float:
    """``log[P(N >= 1) / P(N = 0)] = log(Z - 1)`` — exact, since ``E(empty) = 0``.

    The whole-interface reading of the two-state contrast that ``eta_a`` is for a single site: the
    free energy of making any contact at all against making none.
    """
    return float(log_z + np.log1p(-np.exp(-log_z)))


def delta_f_threshold(totals: np.ndarray, x: int) -> float:
    """``log[P(N >= x) / P(N < x)]`` under the model. ``Z`` cancels, so this needs no AIS."""
    if not len(totals):
        return float("nan")
    hi = float((totals >= x).mean())
    if hi <= 0.0 or hi >= 1.0:
        return float("inf") if hi >= 1.0 else float("-inf")
    return float(np.log(hi / (1.0 - hi)))
