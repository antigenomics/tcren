"""Scoring a structure under a fitted model: energy, ``log Z``, likelihood, contact probabilities.

Every quantity here is per structure and needs the model's parameters plus that structure's own
available-pair set. Nothing is fitted.
"""

from __future__ import annotations

import hashlib
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import polars as pl

from .kernel import bucket_edges, colour, coupling_matrix, edges
from .model import PottsModel, kernel_names
from .sample import (ais_log_z, count_free_energy, delta_f_empty, delta_f_threshold, energy,
                     factorised_log_z, gibbs, mu_star)
from .sites import available_pairs, eta, site_codes


def _rng_for(seed: int, pdb_id) -> np.random.Generator:
    """A generator determined by ``(seed, pdb_id)`` alone — never by table position.

    Every sampler below used to share one generator across its loop over structures, so a
    structure's ``log Z`` depended on **how many structures preceded it in the table**: the same
    PDB scored on its own, in a reordered frame, or in a subset came back with a different value.
    Deriving the stream from the structure's own identifier makes every sampled quantity a function
    of the structure, which is what makes a per-structure score reproducible. ``blake2b`` rather
    than :func:`hash`, whose seed is randomised per interpreter.
    """
    h = int.from_bytes(hashlib.blake2b(str(pdb_id).encode(), digest_size=8).digest(), "big")
    return np.random.default_rng(np.random.SeedSequence([int(seed), h]))


def _prepare(sites: pl.DataFrame, model: PottsModel):
    """Per-structure arrays: ``(pdb ids, starts, eta, sigma, edge table, offsets, coefficients)``."""
    # Sorted on the site's own identity, not on arrival order. The colouring, the edge indices and
    # therefore every sampled quantity are functions of this order, so leaving it to the caller made
    # `log Z` depend on how the pair table happened to be concatenated. With `_rng_for` below, this
    # is the second half of "a per-structure score is a function of the structure".
    sites = (sites.filter(pl.col("d_ca") <= model.radius)
             .sort("pdb.id", "chain.rec", "region.rec", "pos.rec", "pos.par"))
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


#: Filled once per worker process by :func:`_init_workers`, so the prepared arrays cross the
#: process boundary once rather than once per structure.
_SHARED: dict = {}


def _init_workers(payload: dict) -> None:
    _SHARED.clear()
    _SHARED.update(payload)


def _run_chunk(args):
    body, idx, kw = args
    return [_BODIES[body](s, **kw) for s in idx]


def _map_structures(body: str, n: int, payload: dict, kw: dict, workers: int | None):
    """Run one per-structure body over every structure, in contiguous chunks.

    Safe to parallelise only because each result is a function of ``(seed, pdb.id)`` and of that
    structure's own sites — never of its position in the frame — so the output is identical however
    the work is split. Chunks are contiguous and as few as there are workers: a pool of thousands of
    one-structure tasks spends its time in dispatch and pickling.

    Measured on an M3 (16 cores), 64 structures of 20--80 sites: serial 48.9 ms per structure,
    3.6x at 8 processes. Threads make it *slower* (0.33x) — the arrays are too small for numpy to
    release the GIL usefully.
    """
    if workers is None:
        workers = os.cpu_count() or 1
    if workers <= 1 or n < 2 * workers:
        _init_workers(payload)
        return _run_chunk((body, list(range(n)), kw))
    bounds = np.linspace(0, n, workers + 1).astype(int)
    chunks = [(body, list(range(bounds[i], bounds[i + 1])), kw)
              for i in range(workers) if bounds[i + 1] > bounds[i]]
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_workers,
                             initargs=(payload,)) as ex:
        return [row for part in ex.map(_run_chunk, chunks) for row in part]


def _payload(upid, starts, e, sigma, E, eoff, kv) -> dict:
    return {"upid": upid, "starts": starts, "e": e, "sigma": sigma, "E": E, "eoff": eoff, "kv": kv}


def _unpack(s: int):
    """``(pdb id, lo, hi, rng, eta, sigma, A, colours)`` for structure ``s`` from the shared arrays."""
    d = _SHARED
    lo, hi = d["starts"][s], d["starts"][s + 1]
    et, sg, A, cols = _one(lo, hi, d["e"], d["sigma"], d["E"], d["eoff"], d["kv"], s)
    return d["upid"][s], lo, hi, et, sg, A, cols


def _body_score_sites(s: int, *, particles: int, steps: int, seed: int) -> dict:
    pid, lo, hi, et, sg, A, cols = _unpack(s)
    log_z, ess = ais_log_z(et, A, cols, _rng_for(seed, pid), particles=particles, steps=steps)
    en = energy(sg, et, A)
    f = et + A @ sg
    pll = float(np.sum(sg * f - np.logaddexp(0.0, f)))
    return {"pdb.id": pid, "n_sites": int(hi - lo), "n_contacts": int(sg.sum()),
            "energy": en, "neg_energy": -en,
            "log_z": log_z, "log_z0": factorised_log_z(et),
            "log_lik": -en - log_z, "psi": (-en - log_z) / (hi - lo),
            "pseudo_log_lik": pll, "psi_pseudo": pll / (hi - lo),
            "ais_ess": ess, "n_colours": len(cols)}


def _body_bound_unbound(s: int, *, threshold, chains: int, burn: int, draws: int, thin: int,
                        particles: int, steps: int, seed: int) -> dict:
    pid, lo, hi, et, sg, A, cols = _unpack(s)
    rng = _rng_for(seed, pid)
    _, tot = gibbs(et, A, cols, rng, chains=chains, burn=burn, draws=draws, thin=thin)
    lz, _ = ais_log_z(et, A, cols, rng, particles=particles, steps=steps)
    n_obs = float(sg.sum())
    return {"pdb.id": pid, "n_sites": hi - lo, "n_contacts": n_obs,
            "neg_energy": -energy(sg, et, A),
            "log_z": lz, "df_empty": delta_f_empty(lz),
            "df_threshold": (float("nan") if threshold is None
                             else delta_f_threshold(tot, threshold)),
            "mu_star": mu_star(tot, n_obs),
            "n_lo": float(tot.min()) if len(tot) else float("nan"),
            "n_hi": float(tot.max()) if len(tot) else float("nan"),
            "n_mean": float(tot.mean()) if len(tot) else float("nan"),
            "n_var": float(tot.var()) if len(tot) else float("nan")}


def _body_contact_probabilities(s: int, *, chains: int, burn: int, draws: int, thin: int,
                                seed: int) -> dict:
    pid, lo, hi, et, sg, A, cols = _unpack(s)
    occ, _ = gibbs(et, A, cols, _rng_for(seed, pid), chains=chains, burn=burn, draws=draws,
                   thin=thin)
    # `lo`/`hi` travel back with the arrays so the caller scatters them into the right slice of
    # the full-table buffer; a worker cannot see the frame it came from.
    return {"lo": int(lo), "hi": int(hi), "p_model": occ,
            "p_conditional": 1.0 / (1.0 + np.exp(-(et + A @ sg)))}


_BODIES = {"score_sites": _body_score_sites, "bound_unbound": _body_bound_unbound,
           "contact_probabilities": _body_contact_probabilities}


def score_sites(sites: pl.DataFrame, model: PottsModel, *, particles: int = 64,
                steps: int = 256, seed: int = 0, workers: int | None = None) -> pl.DataFrame:
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
        workers: Processes to split the structures across. ``None`` takes every core, ``1``
            runs serially. The result is identical either way — each structure's numbers
            depend on ``(seed, pdb.id)`` alone, never on how the work was split.
    """
    upid, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    rows = _map_structures("score_sites", len(upid), _payload(upid, starts, e, sigma, E, eoff, kv),
                           {"particles": particles, "steps": steps, "seed": seed}, workers)
    return pl.DataFrame(rows)


def contact_probabilities(sites: pl.DataFrame, model: PottsModel, *, chains: int = 64,
                          burn: int = 100, draws: int = 100, thin: int = 3,
                          seed: int = 0, workers: int | None = None) -> pl.DataFrame:
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
        workers: Processes to split the structures across. ``None`` (default) uses every core,
            ``1`` runs serially. Bit-identical either way; see :func:`tcren.potts.score_sites`.

    Returns:
        The input frame with the three probabilities appended, one row per site.
    """
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    rows = _map_structures("contact_probabilities", len(upid),
                           _payload(upid, starts, e, sigma, E, eoff, kv),
                           dict(chains=chains, burn=burn, draws=draws, thin=thin, seed=seed),
                           workers)
    p_mod, p_cond = np.zeros(len(e)), np.zeros(len(e))
    for r in rows:
        p_mod[r["lo"]:r["hi"]] = r["p_model"]
        p_cond[r["lo"]:r["hi"]] = r["p_conditional"]
    return q.drop("sid", "loop", "pchain").with_columns(
        pl.Series("eta", e),
        pl.Series("p_independent", 1.0 / (1.0 + np.exp(-e))),
        pl.Series("p_model", p_mod),
        pl.Series("p_conditional", p_cond))



#: What each ``by=`` groups on, beyond ``pdb.id``. ``"pair"`` is the ungrouped passthrough.
_MAP_KEYS = {"pair": (),
             "loop": ("region.rec", "pos.par", "aa.par"),
             "position": ("pos.par", "aa.par")}


def contact_map(sites: pl.DataFrame, model: PottsModel, *, by: str = "loop", chains: int = 64,
                burn: int = 100, draws: int = 100, thin: int = 3, seed: int = 0,
                workers: int | None = None) -> pl.DataFrame:
    r"""Predicted contact frequency, closed from per-pair probabilities onto a coarser grid.

    :func:`contact_probabilities` gives ``p_model``, the coupled model's marginal for a single
    receptor-residue : peptide-residue pair. Two coarser readings are what an experiment actually
    measures:

    ``by="loop"``
        one row per (structure, CDR loop, peptide position) — the **contact-frequency map**, the
        grid a molecular-dynamics trajectory reports as the fraction of frames in which any residue
        of that loop touches that peptide position.
    ``by="position"``
        one row per (structure, peptide position) — **peptide residue importance**: how engaged the
        model expects that position to be, before any residue identity is scored.
    ``by="pair"``
        the ungrouped table, exactly :func:`contact_probabilities`.

    The residues of a loop are distinct pairs with different probabilities, so the number of
    simultaneous contacts is Poisson-binomially distributed and has no closed form. The event "at
    least one" does, and it is what a frequency map measures:

    .. math:: P(N \ge 1) \;=\; 1 - \prod_j (1 - p_j)

    where :math:`p_j` is the model marginal of pair :math:`j` in the group and :math:`N` the number
    of contacts the group makes. It is accumulated in :math:`\log(1 - p)` so a twelve-residue loop
    does not underflow, and a pair at :math:`p_j = 1` forces the group to 1 exactly rather than to
    ``nan``.

    These are contact frequencies: dimensionless, in :math:`[0, 1]`, higher meaning more often in
    contact. They are **not** energies and carry no :math:`k_\mathrm{B}T`, so nothing here belongs
    in an energy block — this is the diagnostic and importance side of the model.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`.
        model: A fitted :class:`PottsModel`.
        by: ``"loop"`` (default), ``"position"`` or ``"pair"``.
        chains, burn, draws, thin: Gibbs settings; see :func:`tcren.potts.gibbs`.
        seed: Seed for the sampler.
        workers: Processes to split the structures across; ``None`` uses every core, ``1`` runs
            serially. Bit-identical either way.

    Returns:
        For ``"pair"``, :func:`contact_probabilities`' frame unchanged. Otherwise the grouping
        columns plus ``p_any`` (the predicted frequency above), ``p_expected`` (the expected number
        of contacts in the group, :math:`\sum_j p_j`), ``n_pairs`` (available pairs in the group),
        ``n_observed`` (how many of them this structure made) and ``observed`` (1 if any did, else
        0 — the indicator the prediction is scored against).

    Example:
        >>> from tcren.potts import available_pairs, contact_map, PottsModel  # doctest: +SKIP
        >>> m = contact_map(available_pairs(structure), PottsModel.bundled())  # doctest: +SKIP
    """
    if by not in _MAP_KEYS:
        raise ValueError(f"by must be one of {'|'.join(_MAP_KEYS)}, got {by!r}")
    p = contact_probabilities(sites, model, chains=chains, burn=burn, draws=draws, thin=thin,
                              seed=seed, workers=workers)
    if by == "pair" or p.is_empty():
        return p
    keys = ["pdb.id", *_MAP_KEYS[by]]
    return (p.group_by(keys)
            .agg([(1.0 - (-pl.col("p_model")).log1p().sum().exp()).alias("p_any"),
                  pl.col("p_model").sum().alias("p_expected"),
                  pl.len().cast(pl.Int64).alias("n_pairs"),
                  pl.col("sigma").sum().cast(pl.Int64).alias("n_observed"),
                  (pl.col("sigma").max() > 0).cast(pl.Int64).alias("observed")])
            .sort(keys))


def peptide_free_energy(sites: pl.DataFrame, model: PottsModel, *, coupled: bool = False,
                        marginals: pl.DataFrame | None = None, chains: int = 64, burn: int = 100,
                        draws: int = 100, thin: int = 3, seed: int = 0,
                        workers: int | None = None) -> pl.DataFrame:
    r"""Free energy of the interface with each residue threaded through each partner position.

    :func:`contact_map` reads how engaged a position is expected to be *before any residue identity
    is scored*. This reads what happens when the identity changes. The partner residue enters the
    one-body field twice — through the partner propensity :math:`h^{\mathrm{par}}` and through the
    pair term :math:`J` — so substituting position :math:`i` shifts :math:`\eta` at every available
    pair carrying that position, and the interface free energy moves with it:

    .. math::

       \Phi^{\mathrm{Potts}}(x) \;=\; \log Z_0\big(\eta(x)\big)
       \;=\; \sum_a \log\!\big(1 + e^{\eta_a(x)}\big),
       \qquad
       \Delta F_i(a) \;=\; \Phi^{\mathrm{Potts}}(x_{i \to a})
                          \;-\; \tfrac{1}{20}\sum_b \Phi^{\mathrm{Potts}}(x_{i \to b})

    Higher is more favourable: :math:`\log Z_0` is the interface's capacity to make contacts at all,
    so a residue that raises it engages more. The reference is the **equimolar** one — the mean over
    the twenty residues at that position, not the residue the structure happens to carry — which is
    the null a positional-scanning library actually holds the other positions at.

    Unlike :func:`contact_map`'s frequencies this **is** an energy: :math:`\log Z_0` carries
    :math:`k_\mathrm{B}T` and belongs in an energy block.

    Two readings, from the same fields:

    ``coupled=False`` (default)
        :math:`\log Z_0` for the coupling-free model, which is exact and closed form — no sampling.
    ``coupled=True``
        linear response about the observed sequence. Since
        :math:`\partial \log Z / \partial \eta_a = \langle\sigma_a\rangle`, the coupled free energy
        moves as :math:`\Delta \log Z \approx \sum_a p_a \Delta\eta_a` with :math:`p_a` the marginal
        of :func:`contact_probabilities` — one Gibbs pass, then a dot product per cell.

    Only ``aa.par`` changes: the backbone, the Cα distances, the receptor residues and the partner
    roles are the structure's own and are held fixed, which is the same fixed-backbone approximation
    every threading score in the package makes.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`.
        model: A fitted :class:`PottsModel`.
        coupled: Take the linear-response path against the coupled marginals.
        marginals: A frame from :func:`contact_probabilities` to reuse when ``coupled=True``.
            ``None`` computes it.
        chains, burn, draws, thin: Gibbs settings, used only when ``coupled=True``.
        seed: Seed for the sampler.
        workers: Processes to split the structures across when ``coupled=True``.

    Returns:
        One row per (``pdb.id``, ``pos.par``, ``aa.par``) with ``log_z0`` (the whole-interface
        :math:`\log Z_0` under that substitution), ``dF`` (its equimolar-referenced effect),
        ``n_pairs`` (available pairs carrying that position) and ``is_observed`` (1 for the residue
        the structure carries). ``dF`` sums to zero over the twenty residues at every position.

    Example:
        >>> from tcren.potts import available_pairs, peptide_free_energy, PottsModel  # doctest: +SKIP
        >>> peptide_free_energy(available_pairs(structure), PottsModel.bundled())  # doctest: +SKIP
    """
    if sites.is_empty():
        return pl.DataFrame(schema={"pdb.id": pl.String, "pos.par": pl.Int64, "aa.par": pl.String,
                                    "log_z0": pl.Float64, "dF": pl.Float64,
                                    "n_pairs": pl.Int64, "is_observed": pl.Int64})
    aa = tuple(model.alphabet)
    codes, _sizes, q = site_codes(sites, model)
    # A partner position is one residue of one chain, so every site carrying it must agree on the
    # identity. Averaging over a disagreement would silently score a sequence that does not exist.
    clash = (q.group_by("pdb.id", "pos.par").agg(pl.col("aa.par").n_unique().alias("n"))
             .filter(pl.col("n") > 1))
    if not clash.is_empty():
        raise ValueError(f"{clash.height} partner position(s) carry more than one residue, so "
                         f"there is no sequence to substitute into: {clash.head(3).to_dicts()}")
    # eta splits into a part the partner residue does not touch and a part it does. Recovering the
    # first as eta(observed) minus the second keeps ONE definition of eta in the package: any change
    # to `sites.eta` propagates here rather than being silently re-implemented.
    J, h_par = model.coupling_array(), np.asarray(model.h_par)
    rest = eta(codes, model) - h_par[codes[1]] - J[codes[0], codes[1]]
    # (n_sites, 20): the field at every site under every candidate partner residue
    eta_all = rest[:, None] + h_par[None, :] + J[codes[0], :]

    if coupled:
        p = (marginals if marginals is not None else
             contact_probabilities(sites, model, chains=chains, burn=burn, draws=draws, thin=thin,
                                   seed=seed, workers=workers))["p_model"].to_numpy()
        if p.shape != (q.height,):
            raise ValueError(f"marginals must have one row per site ({q.height}), got {p.shape}")
        # log Z is linear in eta to first order with slope <sigma>; the constant cancels in dF
        per_site = p[:, None] * eta_all
    else:
        per_site = np.logaddexp(0.0, eta_all)

    # what each site contributes at the residue the structure actually carries there
    at_observed = per_site[np.arange(q.height), codes[1]]

    out = []
    for (pid,), g in q.with_row_index("_r").group_by(["pdb.id"], maintain_order=True):
        whole = float(at_observed[g["_r"].to_numpy()].sum())
        for pos, gp in g.group_by(["pos.par"], maintain_order=True):
            rp = gp["_r"].to_numpy()
            # sites away from this position keep their own residues and contribute the same term
            # whatever this position carries, so they enter as one constant and drop out of dF
            here = per_site[rp].sum(axis=0)
            z = whole - float(at_observed[rp].sum()) + here
            obs = gp["aa.par"][0]
            out += [{"pdb.id": pid, "pos.par": int(pos[0]), "aa.par": a,
                     "log_z0": float(z[k]), "dF": float(here[k] - here.mean()),
                     "n_pairs": len(rp), "is_observed": int(a == obs)}
                    for k, a in enumerate(aa)]
    return pl.DataFrame(out).sort("pdb.id", "pos.par", "aa.par")


def bound_unbound(sites: pl.DataFrame, model: PottsModel, *, threshold: int | None = None,
                  chains: int = 64, burn: int = 100, draws: int = 100, thin: int = 3,
                  particles: int = 64, steps: int = 256, seed: int = 0,
                  workers: int | None = None) -> pl.DataFrame:
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
        workers: Processes to split the structures across. ``None`` takes every core, ``1``
            runs serially. The result is identical either way — each structure's numbers
            depend on ``(seed, pdb.id)`` alone, never on how the work was split.

    Returns:
        One row per structure.
    """
    upid, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    out = _map_structures("bound_unbound", len(upid), _payload(upid, starts, e, sigma, E, eoff, kv),
                          {"threshold": threshold, "chains": chains, "burn": burn, "draws": draws,
                           "thin": thin, "particles": particles, "steps": steps, "seed": seed},
                          workers)
    return pl.DataFrame(out)


def count_profile(sites: pl.DataFrame, model: PottsModel, *, chains: int = 64, burn: int = 100,
                  draws: int = 100, thin: int = 3, seed: int = 0) -> pl.DataFrame:
    """Pooled ``F(N) = -log p(N)`` over every structure in ``sites``, plus the observed counts.

    The free-energy profile along the contact-count coordinate. Use it to see whether the model's
    contact-count landscape has a barrier -- if it does not, a threshold reading of the two-state
    contrast has nothing to key on and ``mu_star`` is the meaningful statistic.
    """
    upid, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    keep, obs = [], []
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        rng = _rng_for(seed, upid[s])
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
    upid, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    out = []
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        rng = _rng_for(seed, upid[s])
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
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    n_k = len(kv)
    if not n_k:
        return pl.DataFrame(schema={"class": pl.Utf8, "n_edges": pl.Int64,
                                    "c_data": pl.Float64, "c_model": pl.Float64})
    acc = {k: np.zeros(n_k) for k in ("m2", "m1a", "m1b", "d2", "d1a", "d1b", "dn")}
    m_n = np.zeros(n_k)
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        rng = _rng_for(seed, upid[s])
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
