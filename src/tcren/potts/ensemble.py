r"""The contact map as an ensemble: non-linear energy read-outs and the configurations themselves.

:func:`tcren.potts.contact_probabilities` returns the per-site marginal
:math:`p_a = \langle\sigma_a\rangle`. That is enough for the *mean* interaction energy and for
nothing else. Because :math:`\Phi(\sigma) = \sum_a \sigma_a \varphi_a` is linear in
:math:`\sigma`, linearity of expectation gives

.. math::  \langle\Phi\rangle \;=\; \sum_a \langle\sigma_a\rangle\, \varphi_a

for **any** joint distribution, coupled or not. Dependence does not break the identity; it fixes
what :math:`\langle\sigma_a\rangle` has to be — the marginal of the *coupled* model
(``p_model``), not the factorised one (``p_independent``) — and it fixes the sequence at which the
marginal is evaluated. Everything else the ensemble knows is invisible to the marginals:

``phi_var``
    :math:`\operatorname{Var}\Phi`, which needs the two-point function
    :math:`\langle\sigma_a\sigma_b\rangle` and is therefore zero under any independent model with
    the same marginals only when the couplings are.
``phi_free``
    :math:`-\beta^{-1}\log\langle e^{-\beta\Phi}\rangle`, a genuine free energy of the threading.
    It is a soft minimum over configurations: :math:`\beta \to 0` recovers
    :math:`\langle\Phi\rangle` and :math:`\beta \to \infty` the best map the geometry admits.
``map_entropy``
    :math:`-\langle \log P(\sigma)\rangle = \langle E\rangle + \log Z`, how many contact maps the
    interface is really choosing between.

Typical use::

    from tcren.potts import PottsModel, available_pairs
    from tcren.potts.ensemble import energy_ensemble, map_ensemble, threaded_ensemble

    pairs = available_pairs(structure)
    model = PottsModel.bundled()
    energy_ensemble(pairs, model)                    # one row per structure
    map_ensemble(pairs, model, top_k=20)             # the configurations, with probabilities
    threaded_ensemble(pairs, model, candidates)      # per candidate, marginals recomputed

The CLI mirrors it: ``tcren potts ensemble``.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from ..potential import Potential, tcren2
from .model import PottsModel
from .sample import ais_log_z, gibbs
from .sites import site_codes
from .score import (_BODIES, _map_structures, _one as _one_slice, _payload, _prepare,
                    _rng_for, _unpack)

__all__ = ["energy_ensemble", "map_ensemble", "threaded_ensemble", "site_energies"]

#: Columns `_prepare` adds to the frame it returns; dropped before re-deriving codes so that
#: `site_codes`' own `with_row_index("sid")` does not collide with the one already there.
_ADDED = ("sid", "loop", "pchain")


def site_energies(frame: pl.DataFrame, potential: Potential) -> tuple[np.ndarray, int]:
    r"""Per-site interaction energy :math:`\varphi_a` = ``potential(aa.rec, aa.par)``.

    A pair the potential does not define — TCRen2 carries no ``from = C`` row, so it is 380 cells
    and not 400 — contributes **zero**, matching :func:`tcren.energetics.scoring._interface_energy`'s
    ``nansum``. The count of such sites is returned rather than swallowed, because a silent zero on
    a large fraction of the interface is the failure mode this whole path exists to avoid.

    Returns:
        ``(phi, undefined)`` — the energies, and the per-site boolean mask of the sites the
        potential does not define, so a caller can count them per structure rather than in bulk.
    """
    dense, index = potential.as_matrix()
    rec = np.array([index.get(a, -1) for a in frame["aa.rec"]], dtype=np.int64)
    par = np.array([index.get(b, -1) for b in frame["aa.par"]], dtype=np.int64)
    ok = (rec >= 0) & (par >= 0)
    phi = np.zeros(len(rec))
    phi[ok] = dense[rec[ok], par[ok]]
    undefined = ~ok | ~np.isfinite(phi)
    phi[undefined] = 0.0
    return phi, undefined


def _accumulate(et, A, cols, rng, phi, sigma_obs, *, chains, burn, draws, thin, init,
                groups=None, n_groups=0):
    """Run one Gibbs chain set, accumulating whole-configuration statistics as they are drawn.

    The draws are never materialised: `gibbs`' observer hands back the ``(chains, n)`` matrix of
    each kept sweep and every statistic here is a reduction over it. `groups` additionally packs a
    **coarse-grained** pattern per draw — whether each group of sites engaged at all — which is the
    only granularity at which a configuration recurs often enough to carry a probability (see
    :func:`map_ensemble`).
    """
    phi_vals: list[np.ndarray] = []
    e_vals: list[np.ndarray] = []
    ham: list[np.ndarray] = []
    packed: list[np.ndarray] = []
    coupled = A.any()

    def observe(sig):
        phi_vals.append(sig @ phi)
        e_vals.append(-(sig @ et) - (0.5 * ((sig @ A) * sig).sum(1) if coupled else 0.0))
        ham.append(np.abs(sig - sigma_obs).sum(1))
        if groups is not None:
            g = np.zeros((sig.shape[0], n_groups))
            np.add.at(g.T, groups, sig.T)
            packed.append(np.packbits((g > 0).astype(np.uint8), axis=1))

    occ, totals = gibbs(et, A, cols, rng, chains=chains, burn=burn, draws=draws, thin=thin,
                        observer=observe, init=init)
    return (occ, totals,
            np.concatenate(phi_vals) if phi_vals else np.zeros(0),
            np.concatenate(e_vals) if e_vals else np.zeros(0),
            np.concatenate(ham) if ham else np.zeros(0),
            np.concatenate(packed) if packed else np.zeros((0, 0), np.uint8))


def _free_energy(phi_vals: np.ndarray, beta: float) -> float:
    r""":math:`-\beta^{-1}\log\langle e^{-\beta\Phi}\rangle`, by logsumexp so it survives the tails."""
    if not len(phi_vals) or beta <= 0:
        return float(phi_vals.mean()) if len(phi_vals) else float("nan")
    x = -beta * phi_vals
    m = float(x.max())
    return -(m + float(np.log(np.mean(np.exp(x - m))))) / beta


def _summary(pid, n_sites, phi_vals, e_vals, ham, sigma_obs, log_z, beta,
             n_undefined) -> dict:
    return {
        "pdb.id": pid,
        "n_sites": int(n_sites),
        "n_contacts_observed": float(sigma_obs.sum()),
        "n_undefined": n_undefined,
        "phi_mean": float(phi_vals.mean()) if len(phi_vals) else float("nan"),
        "phi_var": float(phi_vals.var()) if len(phi_vals) else float("nan"),
        "phi_free": _free_energy(phi_vals, beta),
        "map_entropy": (float(e_vals.mean()) + log_z) if len(e_vals) else float("nan"),
        "n_draws": int(len(phi_vals)),
        "hamming_obs": float(ham.mean()) if len(ham) else float("nan"),
    }


def _body_energy_ensemble(s: int, *, beta: float, anchor: str, chains: int, burn: int,
                          draws: int, thin: int, particles: int, steps: int, seed: int) -> dict:
    pid, lo, hi, et, sg, A, cols = _unpack(s)
    phi = _SHARED_PHI()[lo:hi]
    rng = _rng_for(seed, pid)
    init = sg if anchor == "observed" else None
    occ, _, phi_vals, e_vals, ham, _ = _accumulate(
        et, A, cols, rng, phi, sg, chains=chains, burn=burn, draws=draws, thin=thin, init=init)
    log_z, _ = ais_log_z(et, A, cols, rng, particles=particles, steps=steps)
    row = _summary(pid, hi - lo, phi_vals, e_vals, ham, sg, log_z, beta, _SHARED_UNDEF(lo, hi))
    # The identity the module exists to enforce: <Phi> must equal sum_a p_a phi_a with p_a the
    # coupled marginal. Carried out so a caller can check it without a second Gibbs run.
    row["phi_mean_from_marginals"] = float(occ @ phi)
    row["log_z"] = log_z
    return row


def _SHARED_PHI():
    from .score import _SHARED
    return _SHARED["phi"]


def _SHARED_UNDEF(lo: int, hi: int) -> int:
    from .score import _SHARED
    return int(_SHARED["undefined"][lo:hi].sum())


_BODIES["energy_ensemble"] = _body_energy_ensemble


def energy_ensemble(sites: pl.DataFrame, model: PottsModel, potential: Potential | None = None, *,
                    beta: float = 1.0, anchor: str = "observed", chains: int = 64,
                    burn: int = 100, draws: int = 100, thin: int = 3, particles: int = 64,
                    steps: int = 256, seed: int = 0,
                    workers: int | None = None) -> pl.DataFrame:
    r"""Ensemble read-outs of the interaction energy, one row per structure.

    The mean is the number :func:`tcren.potts.contact_probabilities` already determines; the
    variance, the free energy and the entropy are not, and are what a joint model buys over a set
    of marginals.

    ================================ ====================== =====================================
    column                           units                  what it is
    ================================ ====================== =====================================
    ``phi_mean``                     potential units        :math:`\langle\Phi\rangle`
    ``phi_mean_from_marginals``      potential units        :math:`\sum_a p_a \varphi_a`; equal to
                                                            ``phi_mean`` up to Monte-Carlo error
    ``phi_var``                      (potential units)²     :math:`\operatorname{Var}\Phi`
    ``phi_free``                     potential units        :math:`-\beta^{-1}\log\langle e^{-\beta\Phi}\rangle`
    ``map_entropy``                  nats                   :math:`\langle E\rangle + \log Z`
    ``hamming_obs``                  count                  mean distance from the observed map
    ``n_undefined``                  count                  sites the potential does not define
    ================================ ====================== =====================================

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`.
        model: A fitted :class:`~tcren.potts.PottsModel`.
        potential: The interaction matrix :math:`\varphi` is read from. Default TCRen2.
        beta: Inverse temperature in ``phi_free``. ``beta -> 0`` recovers the mean.
        anchor: ``"observed"`` starts every chain at the structure's own contact map — the
            ensemble of maps that "do not differ much from :math:`\delta_{ij}`" — and
            ``hamming_obs`` reports how far it travels. ``"random"`` starts from the factorised
            model. The stationary distribution is the same; only the burn-in differs.
        chains, burn, draws, thin: Gibbs settings; see :func:`tcren.potts.gibbs`.
        particles, steps: Annealed-importance settings for ``log Z``; see
            :func:`tcren.potts.ais_log_z`.
        seed: Seed for the sampler.
        workers: Processes to split the structures across. Bit-identical either way.

    Note:
        ``anchor="observed"`` is not a Hamming-ball constraint and deliberately so. An
        unconstrained chain started at the observed map is still sampling :math:`P(\sigma)`;
        restricting it to a ball would sample something else and would need its own justification.
        ``hamming_obs`` is the diagnostic such a constraint would have been reached for.
    """
    # ponytail: `anchor="observed"` biases nothing — it only shortens burn-in. If a genuinely
    # constrained ensemble is ever wanted (sample only within Hamming radius r of the observed
    # map), that is a different distribution and needs its own normalisation, not this flag.
    if anchor not in ("observed", "random"):
        raise ValueError(f"anchor must be 'observed' or 'random', got {anchor!r}")
    pot = potential or tcren2()
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    phi, undefined = site_energies(q, pot)
    payload = _payload(upid, starts, e, sigma, E, eoff, kv)
    payload["phi"], payload["undefined"] = phi, undefined
    rows = _map_structures("energy_ensemble", len(upid), payload,
                           dict(beta=beta, anchor=anchor, chains=chains, burn=burn, draws=draws,
                                thin=thin, particles=particles, steps=steps, seed=seed),
                           workers)
    return pl.DataFrame(rows).with_columns(pl.lit(pot.name).alias("potential"))


# ---------------------------------------------------------------------------- the configurations

def _body_map_ensemble(s: int, *, top_k: int, chains: int, burn: int, draws: int, thin: int,
                       seed: int) -> dict:
    from .score import _SHARED
    pid, lo, hi, et, sg, A, cols = _unpack(s)
    phi = _SHARED_PHI()[lo:hi]
    groups = _SHARED["groups"][lo:hi]
    labels = _SHARED["group_labels"][s]
    n_groups = len(labels)
    _, _, phi_vals, e_vals, _, packed = _accumulate(
        et, A, cols, _rng_for(seed, pid), phi, sg, chains=chains, burn=burn, draws=draws,
        thin=thin, init=sg, groups=groups, n_groups=n_groups)
    if not packed.size:
        return {"rows": []}
    uniq, inv, counts = np.unique(packed, axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(-counts)[:top_k]
    obs_g = np.zeros(n_groups)
    np.add.at(obs_g, groups, sg)
    obs_packed = np.packbits((obs_g > 0).astype(np.uint8))
    rows = []
    for rank, k in enumerate(order):
        member = inv == k
        bits = np.unpackbits(uniq[k])[:n_groups].astype(bool)
        rows.append({"pdb.id": pid, "rank": rank,
                     "p_hat": float(counts[k]) / len(phi_vals),
                     "n_groups_engaged": int(bits.sum()),
                     "phi_mean": float(phi_vals[member].mean()),
                     "e_potts_mean": float(e_vals[member].mean()),
                     "is_observed": bool(np.array_equal(uniq[k], obs_packed)),
                     "engaged": ",".join(str(labels[i]) for i in np.flatnonzero(bits)),
                     "pattern": "".join("1" if b else "0" for b in bits)})
    return {"rows": rows, "n_distinct": int(len(uniq)), "n_draws": int(len(phi_vals))}


_BODIES["map_ensemble"] = _body_map_ensemble

#: What each ``by=`` coarse-grains onto, using :func:`tcren.potts.contact_map`'s own vocabulary.
#: ``"pair"`` is the raw site set and is only resolvable on a very small interface.
_GROUP_KEYS = {"pair": ("sid",),
               "loop": ("loop", "pos.par"),
               "position": ("pos.par",)}


def map_ensemble(sites: pl.DataFrame, model: PottsModel, potential: Potential | None = None, *,
                 by: str = "position", top_k: int = 20, chains: int = 64, burn: int = 100,
                 draws: int = 100, thin: int = 3, seed: int = 0,
                 workers: int | None = None) -> pl.DataFrame:
    r"""The plausible contact configurations, each with its probability.

    One row per ``(structure, rank)``, ranked by how often the configuration was drawn: which maps
    the interface is choosing between, how much probability the most likely one carries, whether
    the crystal's own map is among them (``is_observed``), and the mean interaction energy
    *conditional* on each.

    **A configuration has to be coarse-grained to carry a probability at all, and this is
    measured, not assumed.** On 1ao7's TCR:peptide interface (174 available sites, 6,400 draws)
    the raw site-level map recurs never — 6,400 distinct configurations out of 6,400 draws, every
    one at :math:`\hat p = 1.6	imes10^{-4}`, because the state space is :math:`2^{174}`.
    Coarse-grained to which **peptide positions** are engaged at all, the same draws collapse to
    **41 distinct patterns with the most likely at** :math:`\hat p = 0.3219`. So:

    ``by="position"`` (default)
        which peptide positions the receptor engages. 9 groups on a 9-mer; resolvable.
    ``by="loop"``
        which ``(CDR loop, peptide position)`` cells engage — the granularity the model's
        molecular-dynamics validation works at. 39 groups on 1ao7, 3,308 distinct patterns, top
        :math:`\hat p = 0.0216`: at the edge of what 6,400 draws resolve.
    ``by="pair"``
        the raw sites. Only meaningful for an interface small enough to enumerate, where
        :func:`tcren.potts.exact_log_z` gives the exact distribution instead.

    ``pattern`` is the bit string over that structure's groups in ``engaged``-label order;
    ``engaged`` names the groups that are on.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`.
        model: A fitted :class:`~tcren.potts.PottsModel`.
        potential: Read ``phi`` from. Default TCRen2.
        by: Coarse-graining, one of ``"position"``, ``"loop"``, ``"pair"``.
        top_k: Configurations kept per structure, most probable first.
        chains, burn, draws, thin, seed, workers: as :func:`energy_ensemble`.
    """
    if by not in _GROUP_KEYS:
        raise ValueError(f"by must be one of {sorted(_GROUP_KEYS)}, got {by!r}")
    pot = potential or tcren2()
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    phi, undefined = site_energies(q, pot)
    keys = _GROUP_KEYS[by]
    # Ordered on the key columns' own dtypes, not on their string forms. `np.unique` over a
    # concatenated tag sorts lexicographically, which puts peptide position 10 between 1 and 2 and
    # silently scrambles the bit order of every emitted pattern.
    groups = np.zeros(q.height, dtype=np.int64)
    labels = []
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        sl = q[lo:hi].select(keys)
        order = sl.unique(maintain_order=False).sort(list(keys)).with_row_index("__g")
        joined = sl.with_row_index("__r").join(order, on=list(keys), how="left").sort("__r")
        groups[lo:hi] = joined["__g"].to_numpy()
        labels.append(["|".join(str(v) for v in row) for row in order.drop("__g").iter_rows()])
    payload = _payload(upid, starts, e, sigma, E, eoff, kv)
    payload.update(phi=phi, undefined=undefined, groups=groups, group_labels=labels)
    rows = _map_structures("map_ensemble", len(upid), payload,
                           dict(top_k=top_k, chains=chains, burn=burn, draws=draws, thin=thin,
                                seed=seed), workers)
    flat = [r for part in rows for r in part["rows"]]
    if not flat:
        return pl.DataFrame()
    n_dist = {p["rows"][0]["pdb.id"]: p["n_distinct"] for p in rows if p["rows"]}
    return (pl.DataFrame(flat)
            .with_columns(pl.col("pdb.id").replace_strict(n_dist, default=None)
                          .alias("n_distinct_total"),
                          pl.lit(by).alias("by")))


# ------------------------------------------------------------------------ threading a candidate

def _eta_parts(q: pl.DataFrame, model: PottsModel, e: np.ndarray):
    r"""Split ``eta`` into the part the partner sequence does not touch, and the gather that does.

    A partner residue enters ``eta`` exactly twice, through ``h_par`` and through ``J``, so

    .. math::  \eta_a(x) = \underbrace{\eta_a - h^{\mathrm{par}}_{b_a} - J_{a_a b_a}}_{\text{fixed}}
               \;+\; h^{\mathrm{par}}_{b_a(x)} + J_{a_a\, b_a(x)}

    and threading a candidate is a gather, not a refit. Returns
    ``(eta_rest, rec_codes, h_par, J, aa_index, pos_par)``.
    """
    codes, _, _ = site_codes(q.drop(*[c for c in _ADDED if c in q.columns]), model)
    rec, par = codes[0], codes[1]
    h_par = np.asarray(model.h_par)
    J = model.coupling_array()
    eta_rest = e - h_par[par] - J[rec, par]
    aa_index = {a: i for i, a in enumerate(model.alphabet)}
    return eta_rest, rec, h_par, J, aa_index, q["pos.par"].to_numpy()


def threaded_ensemble(sites: pl.DataFrame, model: PottsModel, candidates: pl.DataFrame,
                      potential: Potential | None = None, *, beta: float = 1.0,
                      coupled: bool = True, chains: int = 32, burn: int = 100,
                      candidate_burn: int = 10, draws: int = 50, thin: int = 2,
                      seed: int = 0) -> pl.DataFrame:
    r"""Ensemble read-outs with the contact map re-equilibrated for each candidate sequence.

    The shipped soft-contact path — ``ddg(weights=...)``, ``cpl.predict_response(tcr_weights=...)``
    — takes **one** marginal vector computed on the template and applies it to every candidate,
    giving :math:`\sum_a p_a(\text{template})\,\varphi_a(x)`. That is not the ensemble mean: a
    candidate residue enters :math:`\eta_a` through :math:`h^{\mathrm{par}}` and :math:`J`, so the
    contact ensemble itself moves under substitution. This computes
    :math:`\sum_a p_a(x)\,\varphi_a(x)` and the non-linear read-outs beside it.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`.
        model: A fitted :class:`~tcren.potts.PottsModel`.
        candidates: Frame with ``pdb.id`` and ``peptide``. Every candidate must match the
            structure's own ``par.len``.
        potential: Read ``phi`` from. Default TCRen2.
        beta: Inverse temperature in ``phi_free``.
        coupled: ``True`` re-runs Gibbs per candidate. ``False`` holds the marginals at the
            template and only re-gathers :math:`\varphi`, which is the shipped approximation and is
            here so the two can be measured against each other on identical rows.
        chains, burn, draws, thin: Gibbs settings. ``burn`` applies to the first candidate of a
            structure; later candidates warm-start from the previous one and use
            ``candidate_burn``.
        candidate_burn: Sweeps between candidates once the chain is warm.
        seed: Seed for the sampler.

    Returns:
        One row per ``(pdb.id, peptide)`` with ``phi_mean``, ``phi_var``, ``phi_free``,
        ``n_contacts_model`` and ``n_sites``.
    """
    # ponytail: serial over structures, warm-started over candidates within one. The pool in
    # `_map_structures` splits by structure and would work here too, but the candidate table has to
    # ride along in the payload; add it when a benchmark run actually shows this as the bottleneck.
    pot = potential or tcren2()
    dense, pot_index = pot.as_matrix()
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    eta_rest, rec, h_par, J, aa_index, pos_par = _eta_parts(q, model, e)
    rec_pot = np.array([pot_index.get(a, -1) for a in q["aa.rec"]], dtype=np.int64)
    by_id = {k: v["peptide"].to_list()
             for k, v in candidates.partition_by("pdb.id", as_dict=True).items()}
    by_id = {(k[0] if isinstance(k, tuple) else k): v for k, v in by_id.items()}

    out = []
    for s, pid in enumerate(upid):
        peps = by_id.get(pid)
        if not peps:
            continue
        lo, hi = starts[s], starts[s + 1]
        et0, sg, A, cols = _one_slice(lo, hi, e, sigma, E, eoff, kv, s)
        rng = _rng_for(seed, pid)
        sl = slice(lo, hi)
        pp, rc, rp = pos_par[sl], rec[sl], rec_pot[sl]
        state, first = sg, True
        for pep in peps:
            idx = np.asarray([aa_index.get(a, -1) for a in pep], dtype=np.int64)
            pidx = np.asarray([pot_index.get(a, -1) for a in pep], dtype=np.int64)
            if pp.max(initial=-1) >= len(pep):
                raise ValueError(f"candidate {pep!r} (length {len(pep)}) is shorter than the "
                                 f"{pp.max() + 1} partner positions {pid} carries")
            b = idx[pp]
            if (b < 0).any():
                raise ValueError(f"candidate {pep!r} carries a residue outside the model alphabet")
            et = eta_rest[sl] + h_par[b] + J[rc, b]
            pb = pidx[pp]
            ok = (rp >= 0) & (pb >= 0)
            phi = np.zeros(hi - lo)
            phi[ok] = dense[rp[ok], pb[ok]]
            phi = np.nan_to_num(phi, nan=0.0)
            if coupled:
                occ, _, phi_vals, _, _, _ = _accumulate(
                    et, A, cols, rng, phi, sg, chains=chains,
                    burn=burn if first else candidate_burn, draws=draws, thin=thin, init=state)
                state, first = None, False
                out.append({"pdb.id": pid, "peptide": pep, "n_sites": int(hi - lo),
                            "n_contacts_model": float(occ.sum()),
                            "phi_mean": float(phi_vals.mean()),
                            "phi_var": float(phi_vals.var()),
                            "phi_free": _free_energy(phi_vals, beta)})
            else:
                occ, _, phi_vals, _, _, _ = _accumulate(
                    et0, A, cols, rng, phi, sg, chains=chains,
                    burn=burn if first else candidate_burn, draws=draws, thin=thin, init=state)
                first = False
                out.append({"pdb.id": pid, "peptide": pep, "n_sites": int(hi - lo),
                            "n_contacts_model": float(occ.sum()),
                            "phi_mean": float(occ @ phi),
                            "phi_var": float(phi_vals.var()),
                            "phi_free": _free_energy(phi_vals, beta)})
    return pl.DataFrame(out) if out else pl.DataFrame()
