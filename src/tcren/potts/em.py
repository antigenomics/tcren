r"""Fitting the contact-map model when the contact indicator is not fully observed.

:func:`tcren.potts.fit_potts` treats :math:`\sigma` as data. In a crystal it mostly is — but not
everywhere. A contact is called by thresholding the closest heavy-atom distance at
``cutoff`` (5 Å), and a pair sitting within half an ångström of that line is a **coin toss
recorded as a fact**: the same complex re-refined, or relaxed, would call it the other way. Those
sites are latent, and Monte-Carlo EM is what a latent indicator asks for.

**E-step.** Sample :math:`\sigma_{\mathcal L}` from :math:`P(\sigma_{\mathcal L} \mid
\sigma_{-\mathcal L})` under the current parameters, by block Gibbs with every confident site
**clamped** to what the structure did. The couplings are what make this worth doing: a borderline
pair next to three made contacts is imputed differently from one standing alone, which is exactly
the linkage the model was fitted to carry.

**M-step.** Run the existing penalised pseudolikelihood on the completed data — ``draws``
completions stacked, each at weight :math:`1/\text{draws}`, so the effective sample size and
therefore the ridge are unchanged. Each completion carries its own ``pdb.id`` suffix because every
coupling edge in :mod:`tcren.potts.kernel` joins on that column, and un-suffixed copies would link
sites across draws that never coexisted.

This **adds** a fit; it does not replace one. The shipped pseudolikelihood model is the point of
comparison, and :func:`consistency` is the read-out: which parameter blocks move, measured in units
of the shipped fit's own cluster-robust standard errors.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .fit import DEFAULT_RIDGE, fit_potts
from .model import PottsModel
from .score import _one, _prepare, _rng_for
from .sample import gibbs

__all__ = ["LATENT_MARGIN", "latent_mask", "impute", "fit_potts_em", "consistency"]

#: Half-width, Å, of the band around ``cutoff`` in which the contact call is not trusted. Half the
#: ``DBIN`` the distance profile is binned at: a call the model cannot resolve to better than one
#: bin is not a call.
LATENT_MARGIN: float = 0.5


def latent_mask(frame: pl.DataFrame, *, cutoff: float = 5.0,
                margin: float = LATENT_MARGIN) -> np.ndarray:
    """Boolean mask of the sites whose contact call is inside its own resolution.

    ``|d_heavy - cutoff| <= margin``. A pair at 4.9 Å and one at 5.1 Å are recorded as contact and
    non-contact, and nothing in the coordinates justifies the difference.

    **The band must reach both ways.** ``available_pairs`` measures ``d_heavy`` only out to its
    ``cutoff`` by default, so every pair that missed has a null distance and the band collapses to
    its inward half — an EM run on it could turn contacts off and never on, which is a systematic
    bias and not a small one. Build the pair table with ``available_pairs(..., probe=cutoff +
    margin)`` and this function will see the near misses. It raises rather than proceeding
    one-sided.

    Note:
        A second class of genuinely-latent site — a residue with no resolved side-chain atoms, so
        ``d_heavy`` is a backbone-only lower bound — is not detected here, because the pair table
        does not carry atom counts. Such a site reads as a confident non-contact and the fit
        believes it.
    """
    d = frame["d_heavy"].fill_null(np.inf).to_numpy()
    mask = np.abs(d - cutoff) <= margin
    if mask.any():
        sigma = frame["sigma"].to_numpy()
        if not (sigma[mask] == 0).any():
            raise ValueError(
                f"the latent band [{cutoff - margin}, {cutoff + margin}] A contains "
                f"{int(mask.sum())} sites and every one of them is a contact, so it only reaches "
                "inwards. Rebuild the pair table with "
                f"available_pairs(..., probe={cutoff + margin}) so a pair that missed carries how "
                "far it missed by; otherwise EM can only ever remove contacts.")
    return mask


def impute(sites: pl.DataFrame, model: PottsModel, *, margin: float = LATENT_MARGIN,
           draws: int = 5, chains: int = 1, burn: int = 50, thin: int = 5,
           seed: int = 0) -> tuple[pl.DataFrame, np.ndarray, float]:
    r"""The E-step: ``draws`` completions of the contact map, stacked into one frame.

    Every site outside the latent set keeps the value the structure gave it; the latent ones are
    resampled from their conditional under ``model``. Returns ``(stacked frame, latent mask, mean
    imputed sigma)``; the frame's ``pdb.id`` carries a ``#m`` draw suffix.
    """
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    mask = latent_mask(q, cutoff=model.cutoff, margin=margin)
    completions = np.repeat(sigma[None, :], draws, axis=0)
    for s in range(len(upid)):
        lo, hi = starts[s], starts[s + 1]
        et, sg, A, cols = _one(lo, hi, e, sigma, E, eoff, kv, s)
        m = mask[lo:hi]
        if not m.any():
            continue
        # Restricting each colour class to its latent members IS the conditional: the clamped
        # sites stay at their observed value in `sig` and still enter every field through `A`.
        free = [c[m[c]] for c in cols if m[c].any()]
        kept = []
        gibbs(et, A, free, _rng_for(seed, pid := upid[s]), chains=chains, burn=burn,
              draws=draws, thin=thin, init=sg, observer=lambda x: kept.append(x[0].copy()))
        assert pid == upid[s]
        for d in range(min(draws, len(kept))):
            completions[d, lo:hi] = kept[d]
    base = q.drop(*[c for c in ("sid", "loop", "pchain") if c in q.columns])
    stacked = pl.concat([
        base.with_columns(pl.Series("sigma", completions[d]),
                          (pl.col("pdb.id") + f"#{d}").alias("pdb.id"))
        for d in range(draws)])
    imputed = float(completions[:, mask].mean()) if mask.any() else float("nan")
    return stacked, mask, imputed


def _eta_vector(model: PottsModel) -> np.ndarray:
    """Every free parameter of a model as one vector, for comparing two fits."""
    return np.concatenate([[model.alpha], model.h_rec, model.h_par, model.g_dist,
                           model.g_region, model.g_role, model.g_class,
                           model.coupling_array().ravel(), model.kernel])


def fit_potts_em(sites: pl.DataFrame, *, margin: float = LATENT_MARGIN, draws: int = 5,
                 steps: int = 10, init: PottsModel | None = None, chains: int = 1,
                 burn: int = 50, thin: int = 5, tol: float = 1e-3, seed: int = 0,
                 ridge: float = DEFAULT_RIDGE, notes: str = "",
                 **fit_kwargs) -> tuple[PottsModel, pl.DataFrame]:
    r"""Monte-Carlo EM over the borderline contact calls, with the existing M-step.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`.
        margin: Half-width of the latent band around ``cutoff``, Å. ``0`` makes every site
            observed, in which case this reduces to :func:`tcren.potts.fit_potts` exactly.
        draws: Completions per E-step. Each enters the M-step at weight ``1/draws``.
        steps: Maximum EM iterations.
        init: Starting parameters. ``None`` runs the plain pseudolikelihood fit first, which is
            both the natural warm start and the model the result is compared against.
        chains, burn, thin: Gibbs settings for the conditional sampler.
        tol: Stop when no parameter moves by more than this between iterations.
        seed: Seed for the sampler.
        ridge, notes, fit_kwargs: Passed to :func:`tcren.potts.fit_potts`. A ``weights`` entry is
            per-structure redundancy weighting keyed on the **unsuffixed** ``pdb.id``; it is
            multiplied by the ``1/draws`` each completion enters at.

    Returns:
        ``(model, trace)``. The trace has one row per iteration: ``step``, ``n_latent``,
        ``mean_imputed_sigma``, ``pseudo_loglik``, ``max_shift``, ``max_shift_vs_init`` and
        ``converged``.

    Warning:
        Monte-Carlo EM with a fixed ``draws`` does not converge to a point — it random-walks in a
        ball around the optimum whose radius is set by the Monte-Carlo noise, so ``max_shift`` has
        a floor that ``tol`` may sit below. Read ``max_shift_vs_init`` beside it: a run that is
        jittering has a flat one, a run that is still moving does not. Raise ``draws`` to lower the
        floor. Everything is seeded, so a run is exactly reproducible regardless.

    Note:
        ``pseudo_loglik`` on the returned model is the M-step's value on **completed** data and is
        not comparable to the pseudolikelihood fit's value on observed data — different response.
        Compare the two fits through :func:`consistency`, not through their objectives.
    """
    base = fit_kwargs.pop("weights", None) or {}
    model = init or fit_potts(sites, ridge=ridge, notes=notes, weights=base or None, **fit_kwargs)
    n_struct, n_sites_obs = model.n_structures, model.n_sites
    prev = start = _eta_vector(model)
    latent_contacts = int(sites["sigma"].to_numpy()[
        latent_mask(sites, cutoff=model.cutoff, margin=margin)].sum()) if margin > 0 else 0
    rows, converged = [], False
    for step in range(steps):
        stacked, mask, imputed = impute(sites, model, margin=margin, draws=draws, chains=chains,
                                        burn=burn, thin=thin, seed=seed + step)
        # Each completion enters at 1/draws so the effective sample size, and therefore what the
        # ridge means, is unchanged. A redundancy weighting the caller supplied multiplies it, and
        # is looked up under the UNSUFFIXED id it was built for.
        w = {pid: base.get(pid.rsplit("#", 1)[0], 1.0) / draws
             for pid in stacked["pdb.id"].unique().to_list()}
        model = fit_potts(stacked, weights=w, ridge=ridge, notes=notes, **fit_kwargs)
        cur = _eta_vector(model)
        shift = float(np.max(np.abs(cur - prev)))
        prev = cur
        rows.append({"step": step, "n_latent": int(mask.sum()),
                     "n_latent_contacts": int(latent_contacts),
                     "mean_imputed_sigma": imputed, "pseudo_loglik": model.pseudo_loglik,
                     "max_shift": shift,
                     # Separates "drifting somewhere" from "jittering in place": Monte-Carlo EM
                     # with a fixed `draws` random-walks in a ball around the optimum rather than
                     # settling, so `max_shift` has a noise floor and this is the signal beside it.
                     "max_shift_vs_init": float(np.max(np.abs(cur - start))),
                     "converged": shift < tol})
        if shift < tol:
            converged = True
            break
    # The M-step saw `draws` copies of everything; report what the fit is actually about.
    model.n_structures, model.n_sites = n_struct, n_sites_obs
    model.n_contacts = int(round(model.n_contacts / draws))
    model.notes = (notes or model.notes) + f" | EM: margin={margin}, draws={draws}, " \
                                           f"steps={len(rows)}, converged={converged}"
    return model, pl.DataFrame(rows)


#: Parameter blocks, in the order :func:`_eta_vector` lays them out.
_BLOCKS = ("alpha", "h_rec", "h_par", "g_dist", "g_region", "g_role", "g_class", "J", "kernel")


def consistency(reference: PottsModel, other: PottsModel) -> pl.DataFrame:
    r"""How far a refit moved, per parameter block, in units of the reference's own scale.

    The yardstick is the reference fit's cluster-robust standard error where it has one — the
    coupling kernel, the only block :func:`tcren.potts.fit_potts` stores ``kernel_se`` for. For
    every other block the yardstick is that block's own standard deviation in the reference, which
    answers "did this move by a lot compared with how much it varies at all".

    Returns one row per block: ``n``, ``max_abs``, ``mean_abs``, ``scale``, ``max_scaled`` and
    ``n_beyond_1`` — how many coefficients moved by more than one unit of the yardstick.
    """
    def blocks(m: PottsModel):
        return {"alpha": np.array([m.alpha]), "h_rec": np.asarray(m.h_rec),
                "h_par": np.asarray(m.h_par), "g_dist": np.asarray(m.g_dist),
                "g_region": np.asarray(m.g_region), "g_role": np.asarray(m.g_role),
                "g_class": np.asarray(m.g_class), "J": m.coupling_array().ravel(),
                "kernel": np.asarray(m.kernel)}

    a, b = blocks(reference), blocks(other)
    se = np.asarray(reference.kernel_se) if reference.kernel_se else None
    rows = []
    for name in _BLOCKS:
        d = b[name] - a[name]
        if name == "kernel" and se is not None and len(se) == len(d):
            scale, label = se, "cluster-robust s.e."
        else:
            sd = float(np.std(a[name]))
            # A block with no spread (one level, or a constant covariate) has no yardstick: any
            # shift divided by ~0 reads as enormous and means nothing.
            scale = np.full(len(d), sd if sd > 1e-12 else np.nan)
            label = "block s.d." if sd > 1e-12 else "degenerate"
        scaled = np.abs(d) / scale
        finite = np.isfinite(scaled)
        rows.append({"block": name, "n": len(d), "max_abs": float(np.max(np.abs(d))),
                     "mean_abs": float(np.mean(np.abs(d))),
                     # `scale` is one number per block except for `kernel`, where it is the
                     # per-coefficient standard error; report the mean of the usable entries.
                     "scale": float(scale[finite].mean()) if finite.any() else float("nan"),
                     "yardstick": label,
                     "max_scaled": float(scaled[finite].max()) if finite.any() else float("nan"),
                     "n_beyond_1": int((scaled[finite] > 1.0).sum())})
    return pl.DataFrame(rows)
