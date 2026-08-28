"""Scoring a structure under a fitted model: energy, ``log Z``, likelihood, contact probabilities.

Every quantity here is per structure and needs the model's parameters plus that structure's own
available-pair set. Nothing is fitted.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .kernel import bucket_edges, colour, coupling_matrix, edges
from .model import PottsModel, kernel_names
from .sample import (ais_log_z, count_free_energy, delta_f_empty, delta_f_threshold, energy,
                     factorised_log_z, gibbs, mu_star)
from .sites import available_pairs, eta, site_codes


def _prepare(sites: pl.DataFrame, model: PottsModel):
    """Per-structure arrays: ``(pdb ids, starts, eta, sigma, edge table, offsets, coefficients)``."""
    sites = sites.filter(pl.col("d_ca") <= model.radius).sort("pdb.id", maintain_order=True)
    codes, _, q = site_codes(sites, model)
    e = eta(codes, model)
    sigma = q["sigma"].to_numpy()
    upid, gid = np.unique(q["pdb.id"].to_numpy(), return_inverse=True)
    starts = np.searchsorted(gid, np.arange(len(upid) + 1))
    kv = np.asarray(model.kernel, dtype=float)
    edge_lists = edges(q, joint=model.joint) if kv.size else []
    E, eoff = bucket_edges(edge_lists, starts[:-1], kv) if edge_lists else (
        np.zeros((0, 3), np.int64), np.zeros(len(upid) + 1, np.int64))
    return upid, starts, e, sigma, E, eoff, kv, q


def _one(lo: int, hi: int, e, sigma, E, eoff, kv, s: int):
    """``(eta, sigma, A, colours)`` for one structure, with the colouring asserted valid."""
    n = hi - lo
    Es = E[eoff[s]:eoff[s + 1]]
    ea, eb, ec = Es[:, 0] - lo, Es[:, 1] - lo, Es[:, 2]
    A = coupling_matrix(n, ea, eb, ec, kv) if len(ea) else np.zeros((n, n))
    cols = colour(n, ea, eb) if len(ea) else ([np.arange(n)] if n else [])
    for cidx in cols:
        assert not A[np.ix_(cidx, cidx)].any(), "same-colour sites are coupled"
    return e[lo:hi], sigma[lo:hi], A, cols


def score_sites(sites: pl.DataFrame, model: PottsModel, *, particles: int = 64,
                steps: int = 256, seed: int = 0) -> pl.DataFrame:
    """Energy, ``log Z`` and likelihoods for every structure in a table of available pairs.

    Columns, one row per structure:

    ``n_sites`` / ``n_contacts``
        available pairs and the contacts among them.
    ``energy`` / ``neg_energy``
        ``E(sigma_obs)``, the Hamiltonian of the observed contact map, and its negation. Lower
        ``energy`` is more favourable, so ``neg_energy`` is the one that runs the same way as
        every other score here: higher is more native-like. It is the :math:`\\Pi` block of
        :func:`tcren.reliability.s_free`.
    ``log_z`` / ``log_z0``
        the coupled partition function by AIS, and the factorised one in closed form. ``log_z0``
        is the free energy of the available set *before* any contact is observed.
    ``log_lik`` / ``psi``
        ``log P(sigma_obs)`` and the same divided by ``n_sites``, so it compares across
        interfaces of different size.
    ``pseudo_log_lik`` / ``psi_pseudo``
        the exact, MCMC-free cross-check.
    ``ais_ess``
        effective sample size of the AIS weights, out of ``particles``. Close to ``particles``
        means the annealing schedule was long enough.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`, one or many structures.
        model: A fitted :class:`PottsModel`.
        particles: AIS particles per structure.
        steps: AIS annealing steps.
        seed: Seed for the sampler.
    """
    rng = np.random.default_rng(seed)
    upid, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    rows = []
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        et, sg, A, cols = _one(lo, hi, e, sigma, E, eoff, kv, s)
        log_z, ess = ais_log_z(et, A, cols, rng, particles=particles, steps=steps)
        en = energy(sg, et, A)
        f = et + A @ sg
        pll = float(np.sum(sg * f - np.logaddexp(0.0, f)))
        rows.append({"pdb.id": upid[s], "n_sites": int(hi - lo), "n_contacts": int(sg.sum()),
                     "energy": en, "neg_energy": -en,
                     "log_z": log_z, "log_z0": factorised_log_z(et),
                     "log_lik": -en - log_z, "psi": (-en - log_z) / (hi - lo),
                     "pseudo_log_lik": pll, "psi_pseudo": pll / (hi - lo),
                     "ais_ess": ess, "n_colours": len(cols)})
    return pl.DataFrame(rows)


def contact_probabilities(sites: pl.DataFrame, model: PottsModel, *, chains: int = 64,
                          burn: int = 100, draws: int = 100, thin: int = 3,
                          seed: int = 0) -> pl.DataFrame:
    """Per-site contact probability under the model, beside the observed indicator.

    Three probabilities per site, because the difference between them *is* the couplings:

    ``p_independent``
        ``1 / (1 + exp(-eta))`` — the one-body model alone, ignoring every other site.
    ``p_model``
        the marginal ``<sigma_a>`` of the full coupled model, by block Gibbs. This is the number
        to use.
    ``p_conditional``
        ``P(sigma_a = 1 | the observed rest)`` — the pseudolikelihood conditional, which reads
        "given everything else this structure actually did, how likely was this contact?".

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`.
        model: A fitted :class:`PottsModel`.
        chains, burn, draws, thin: Gibbs settings; see :func:`tcren.potts.gibbs`.
        seed: Seed for the sampler.

    Returns:
        The input frame with the three probabilities appended, one row per site.
    """
    rng = np.random.default_rng(seed)
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    p_mod = np.zeros(len(e))
    p_cond = np.zeros(len(e))
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        et, sg, A, cols = _one(lo, hi, e, sigma, E, eoff, kv, s)
        occ, _ = gibbs(et, A, cols, rng, chains=chains, burn=burn, draws=draws, thin=thin)
        p_mod[lo:hi] = occ
        p_cond[lo:hi] = 1.0 / (1.0 + np.exp(-(et + A @ sg)))
    return q.drop("sid", "loop", "pchain").with_columns(
        pl.Series("eta", e),
        pl.Series("p_independent", 1.0 / (1.0 + np.exp(-e))),
        pl.Series("p_model", p_mod),
        pl.Series("p_conditional", p_cond))


def bound_unbound(sites: pl.DataFrame, model: PottsModel, *, threshold: int | None = None,
                  chains: int = 64, burn: int = 100, draws: int = 100, thin: int = 3,
                  particles: int = 64, steps: int = 256, seed: int = 0) -> pl.DataFrame:
    """The whole-interface two-state free energy, in its three readings.

    A single site has two states, ``sigma_a = 0`` and ``1``, and ``eta_a`` is the free-energy
    difference between them. The same contrast for the whole interface needs a macrostate, and the
    contact count ``N(sigma)`` is the coordinate that defines one. Because every tilt in ``N`` is
    an exponential family, one Gibbs pass serves all three readings.

    Columns, one row per structure:

    ``neg_energy``
        :math:`-E(\\sigma_{\\mathrm{obs}})`, higher is more native-like. The :math:`\\Pi` block of
        :func:`tcren.reliability.s_free`, emitted here so one call supplies both the free-energy
        readings and the energy term the reliability score needs.
    ``df_empty``
        ``log[P(N >= 1) / P(N = 0)] = log(Z - 1)``, exact, from the AIS ``log Z``. ``E(empty) = 0``,
        so the empty configuration needs no separate estimate. This is the *capacity* of the
        interface: what it would gain by making any contact at all.
    ``df_threshold``
        ``log[P(N >= x) / P(N < x)]`` at ``threshold``, from the sampled histogram. ``Z`` cancels,
        so this needs no AIS — but it is only finite while ``x`` sits inside the sampled range.
    ``mu_star``
        the tilt at which ``<N>_mu`` equals the observed count: how much extra chemical potential
        the model needs to explain the map it was shown. Positive means the structure made more
        contacts than its fields and couplings warrant. ``nan`` outside the sampled support.
    ``n_lo``, ``n_hi``, ``n_mean``, ``n_var``
        the sampled contact-count range, mean and variance, so every ``nan`` above is auditable.

    Note:
        ``df_empty`` and ``df_threshold`` are not two estimates of one number. The unbound basin of
        a *docked* pose is astronomically improbable — the model is conditioned on an available set
        that already holds the receptor against the peptide — so no sampler reaches ``N = 0`` and
        only the ``log Z`` route gives it. The threshold reading is a local contrast inside the
        bound basin.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`.
        model: A fitted :class:`PottsModel`.
        threshold: ``x`` for ``df_threshold``. ``None`` disables that column.
        chains, burn, draws, thin: Gibbs settings; see :func:`tcren.potts.gibbs`.
        particles, steps: AIS settings; see :func:`tcren.potts.ais_log_z`.
        seed: Seed for both samplers.

    Returns:
        One row per structure.
    """
    rng = np.random.default_rng(seed)
    upid, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    out = []
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        et, sg, A, cols = _one(lo, hi, e, sigma, E, eoff, kv, s)
        _, tot = gibbs(et, A, cols, rng, chains=chains, burn=burn, draws=draws, thin=thin)
        lz, _ = ais_log_z(et, A, cols, rng, particles=particles, steps=steps)
        n_obs = float(sg.sum())
        out.append({
            "pdb.id": upid[s], "n_sites": hi - lo, "n_contacts": n_obs,
            "neg_energy": -energy(sg, et, A),
            "log_z": lz, "df_empty": delta_f_empty(lz),
            "df_threshold": (float("nan") if threshold is None
                             else delta_f_threshold(tot, threshold)),
            "mu_star": mu_star(tot, n_obs),
            "n_lo": float(tot.min()) if len(tot) else float("nan"),
            "n_hi": float(tot.max()) if len(tot) else float("nan"),
            "n_mean": float(tot.mean()) if len(tot) else float("nan"),
            "n_var": float(tot.var()) if len(tot) else float("nan")})
    return pl.DataFrame(out)


def count_profile(sites: pl.DataFrame, model: PottsModel, *, chains: int = 64, burn: int = 100,
                  draws: int = 100, thin: int = 3, seed: int = 0) -> pl.DataFrame:
    """Pooled ``F(N) = -log p(N)`` over every structure in ``sites``, plus the observed counts.

    The free-energy profile along the contact-count coordinate. Use it to see whether the model's
    contact-count landscape has a barrier -- if it does not, a threshold reading of the two-state
    contrast has nothing to key on and ``mu_star`` is the meaningful statistic.
    """
    rng = np.random.default_rng(seed)
    upid, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    keep, obs = [], []
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        et, sg, A, cols = _one(lo, hi, e, sigma, E, eoff, kv, s)
        _, tot = gibbs(et, A, cols, rng, chains=chains, burn=burn, draws=draws, thin=thin)
        keep.append(tot)
        obs.append(sg.sum())
    n, f = count_free_energy(np.concatenate(keep))
    o = np.asarray(obs, float)
    c = np.array([(o == k).sum() for k in n], float)
    with np.errstate(divide="ignore"):
        f_obs = -np.log(c / max(c.sum(), 1.0))
    return pl.DataFrame({"n_contacts": n, "f_model": f, "f_observed": f_obs,
                         "count_observed": c})


def sample_maps(sites: pl.DataFrame, model: PottsModel, *, chains: int = 64, burn: int = 100,
                draws: int = 100, thin: int = 3, seed: int = 0) -> pl.DataFrame:
    """Contact totals of maps drawn from the model, one row per structure per draw.

    The generative check: a model that reproduces a real interface must reproduce the *spread* of
    its contact count, not only the mean.
    """
    rng = np.random.default_rng(seed)
    upid, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    out = []
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        et, sg, A, cols = _one(lo, hi, e, sigma, E, eoff, kv, s)
        _, totals = gibbs(et, A, cols, rng, chains=chains, burn=burn, draws=draws, thin=thin)
        out.append(pl.DataFrame({"pdb.id": [upid[s]] * len(totals),
                                 "n_contacts_observed": [float(sg.sum())] * len(totals),
                                 "n_contacts_sampled": totals}))
    return pl.concat(out) if out else pl.DataFrame()


def connected_correlations(sites: pl.DataFrame, model: PottsModel, *, chains: int = 64,
                           burn: int = 100, draws: int = 100, thin: int = 3,
                           seed: int = 0) -> pl.DataFrame:
    r"""The generative test: two-point correlations in the data against in maps sampled from it.

    For each coupling class, pooled over every edge in it and every structure,

    .. math::

       c_\Delta = \langle \sigma_a \sigma_b \rangle
                   - \langle \sigma_a \rangle \langle \sigma_b \rangle

    computed once from the observed contact maps and once from maps drawn from the model. This is
    the criterion Boltzmann-machine DCA trains to, and a pseudolikelihood fit never sees it — the
    conditionals it maximises are one-site quantities — so agreement here is a real test rather
    than a restatement of the objective.

    Returns:
        One row per coupling class: ``class``, ``n_edges``, ``c_data``, ``c_model``.
    """
    rng = np.random.default_rng(seed)
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    n_k = len(kv)
    if not n_k:
        return pl.DataFrame(schema={"class": pl.Utf8, "n_edges": pl.Int64,
                                    "c_data": pl.Float64, "c_model": pl.Float64})
    acc = {k: np.zeros(n_k) for k in ("m2", "m1a", "m1b", "d2", "d1a", "d1b", "dn")}
    m_n = np.zeros(n_k)
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        et, sg, A, cols = _one(lo, hi, e, sigma, E, eoff, kv, s)
        Es = E[eoff[s]:eoff[s + 1]]
        ea, eb, ec = Es[:, 0] - lo, Es[:, 1] - lo, Es[:, 2]
        if not len(ea):
            continue
        acc["d2"] += np.bincount(ec, weights=sg[ea] * sg[eb], minlength=n_k)
        acc["d1a"] += np.bincount(ec, weights=sg[ea], minlength=n_k)
        acc["d1b"] += np.bincount(ec, weights=sg[eb], minlength=n_k)
        acc["dn"] += np.bincount(ec, minlength=n_k)
        two, one_a, one_b, tot = _two_point(et, A, cols, rng, ea, eb, ec, n_k,
                                            chains, burn, draws, thin)
        acc["m2"] += two
        acc["m1a"] += one_a
        acc["m1b"] += one_b
        m_n += np.bincount(ec, minlength=n_k) * tot
    dn = np.where(acc["dn"] > 0, acc["dn"], np.nan)
    mn = np.where(m_n > 0, m_n, np.nan)
    return pl.DataFrame({
        "class": kernel_names(model.joint)[:n_k],
        "n_edges": acc["dn"].astype(np.int64),
        "c_data": acc["d2"] / dn - (acc["d1a"] / dn) * (acc["d1b"] / dn),
        "c_model": acc["m2"] / mn - (acc["m1a"] / mn) * (acc["m1b"] / mn),
    })


def _two_point(eta_s, A, cols, rng, ea, eb, ec, n_k, chains, burn, draws, thin):
    """Accumulate ``sum sigma_a sigma_b`` per class over sampled maps, without storing them."""
    n = len(eta_s)
    sig = (rng.random((chains, n)) < 1.0 / (1.0 + np.exp(-eta_s))).astype(np.float64)
    two, one_a, one_b, tot = np.zeros(n_k), np.zeros(n_k), np.zeros(n_k), 0
    coupled = A.any()
    for it in range(burn + draws * thin):
        for cidx in cols:
            f = eta_s[cidx] + (sig @ A[:, cidx] if coupled else 0.0)
            sig[:, cidx] = (rng.random((chains, len(cidx)))
                            < 1.0 / (1.0 + np.exp(-f))).astype(np.float64)
        if it >= burn and (it - burn) % thin == 0:
            np.add.at(two, ec, (sig[:, ea] * sig[:, eb]).sum(0))
            np.add.at(one_a, ec, sig[:, ea].sum(0))
            np.add.at(one_b, ec, sig[:, eb].sum(0))
            tot += chains
    return two, one_a, one_b, tot


def score_structure(structure, model: PottsModel | None = None, *, partner: str = "peptide",
                    **kwargs) -> dict:
    """Convenience: enumerate one structure's available pairs and score them.

    Args:
        structure: A chain-typed structure (MHC-annotated as well, for ``partner="mhc"``).
        model: A fitted model; defaults to the bundled TCR:peptide one.
        partner: ``"peptide"`` or ``"mhc"``.
        **kwargs: Passed to :func:`score_sites`.
    """
    model = model or PottsModel.bundled()
    pairs = available_pairs(structure, partner, radius=model.radius, cutoff=model.cutoff)
    if pairs.is_empty():
        return {"pdb.id": structure.pdb_id, "n_sites": 0, "n_contacts": 0}
    return score_sites(pairs, model, **kwargs).to_dicts()[0]
