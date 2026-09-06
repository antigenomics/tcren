r"""A self-consistent TCRen: the potential re-derived from the contacts the model predicts.

The shipped potential and the shipped contact model run in one direction. `derive_tcren` tallies
the amino-acid pairs of the contacts a crystal **shows**, and `fit_potts` reads the finished matrix
back as the pair term in its field. Nothing carries the model's own view of which contacts formed
back into the tally. This module closes that loop.

.. math::

    \eta_a^{\mathrm{post}} &= \eta_a(\theta^{(t)})
        + \operatorname{clip}\!\left(\frac{c - d_a}{m},\ \pm K\right) \\
    p_a^{(t)} &= \langle \sigma_a \rangle_{P(\sigma \mid \mathrm{obs},\, \theta^{(t)})} \\
    N^{(t)}(a,b) &= \sum_s w_s \sum_{a\, :\, \mathrm{pair}(a) = (a,b)} p_a^{(t)} \\
    e^{(t+1)}(a,b) &= -\ln \frac{N^{(t)}(a,b)\, N_{\mathrm{tot}}}
                                {N_{\mathrm{row}}(a)\, N_{\mathrm{col}}(b)} \\
    \theta^{(t+1)} &= \arg\max \ \mathrm{PL}(\theta \mid e^{(t+1)})

with :math:`c` the 5 Å contact cutoff, :math:`d_a` the site's closest-heavy-atom distance,
:math:`m` the margin, :math:`K` the logit clip, :math:`w_s` the structure's redundancy weight and
:math:`\sigma_a \in \{0,1\}` the contact indicator at site :math:`a` (a receptor residue paired
with a partner residue inside the Cα radius). Higher :math:`e` is less favourable.

**Every conformation is weighted by its probability, and that costs nothing.** The tally the author
asked for is :math:`\sum_\delta P(\delta \mid \mathrm{obs})\, n_{ab}(\delta)`, summed over whole
contact configurations :math:`\delta`. By linearity of expectation that **equals**
:math:`\sum_a p_a`, so the enumeration is never needed to get the counts —
:func:`tcren.potts.map_ensemble` lists the configurations as a read-out, not as a dependency.

**The estimator does not change; only the counts it is fed.** `derive_tcren` already takes a
per-contact fractional weight (`weight_col`), documented there as the way to give a partially
formed contact a fractional vote. So the log-odds, the contact-fraction reference state and the
pseudocount are the shipped ones, unmodified.

**The observation enters as one additive field, which is the whole design.** Conditioning on the
crystal through :math:`P(\sigma_a = 1 \mid d_a) = \mathrm{sigmoid}((c - d_a)/m)` puts a per-site
constant into :math:`\eta`, so the E-step is the existing block-Gibbs sampler called with a
modified vector — no clamping machinery and no second code path. One knob spans the family:

===================  =============================================================================
``margin`` → 0       the crystal map exactly, and therefore **TCRen2 exactly** (the reduction gate)
``margin`` = 0.5 Å   the threshold band only — half of ``DBIN``, matching :data:`tcren.potts.em.LATENT_MARGIN`
``margin`` large     the interface is released and the crystal enters only through the fit
===================  =============================================================================

Pass ``probe = cutoff + 3 * margin`` to :func:`tcren.potts.available_pairs`, or ``d_heavy`` is null
for every non-contact and the field can only ever remove contacts — the same one-sidedness
:func:`tcren.potts.em.latent_mask` refuses to fit.

**This is EM-shaped but it is not EM, and the convergence criterion follows from which.** The
E-step is a genuine conditional expectation. The potential update is not the maximiser of any one
objective — it is the shipped log-odds re-derivation applied to expected counts — so the
pseudo-log-likelihood is **not** guaranteed to rise, and measured on the 17 bundled crystals it does
not (-337.94, -335.05, -335.20, -335.53, -335.68, -335.74 over six steps at ``damping=0.5``: up
once, then drifting down while the matrix is still moving). Convergence is therefore declared on
``max_abs_de``, the largest cell movement in the potential's own natural-log units, which does decay
geometrically at the damping rate. ``pseudo_loglik`` rides along as a diagnostic, not as a
guarantee; treating it as one would be reading an EM property off a fixed-point iteration that does
not have it.

**The data is partly the loop's own output, so the mass is reported rather than assumed.**
``contact_mass`` is the total weighted count against the crystal's own contact total; a loop
inventing contacts shows up there, and ``mass_from_uncertain`` says how much of the tally came from
sites the crystal did not settle.

The framing is Murugan et al., *PNAS* 109:16161 (2012), SI §3.2, where a read's recombination
scenario is latent and each model parameter is updated to the probability-weighted count of the
feature it describes. Their §3.1 keeps the single-best-scenario estimate — "standard practice in the
literature" — explicitly for contrast. **TCRen2 is that deterministic estimate**; what this module
returns is the probabilistic one, and it ships beside rather than instead.

Example:
    >>> import polars as pl
    >>> from tcren.potts import available_pairs, fit_selfconsistent   # doctest: +SKIP
    >>> pairs = pl.concat([available_pairs(s, probe=6.5) for s in structures])  # doctest: +SKIP
    >>> pot, model, trace = fit_selfconsistent(pairs, margin=0.5)     # doctest: +SKIP
    >>> pot.to_csv("TCRen2_sc.csv")                                   # doctest: +SKIP
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..potential import Potential, derive_tcren, tcren2
from .fit import DEFAULT_RIDGE, fit_potts
from .model import PottsModel, kernel_names
from .score import _map_structures, _payload, _prepare

__all__ = ["LOGIT_CLIP", "observation_field", "contact_posterior", "weighted_contacts",
           "fit_selfconsistent"]

#: Where the observation log-odds saturates. ``sigmoid(-40)`` is 4e-18, so a site this far outside
#: the cutoff is a non-contact to well past float64's reach on any realistic contact total, and the
#: clip is what keeps ``margin -> 0`` from being a division by zero rather than a limit.
LOGIT_CLIP: float = 40.0


def observation_field(frame: pl.DataFrame, *, cutoff: float = 5.0, margin: float = 0.5,
                      clip: float = LOGIT_CLIP) -> np.ndarray:
    r"""Per-site log-odds that the crystal's distance read contributes to :math:`\eta`.

    ``(cutoff - d_heavy) / margin``, clipped. Positive inside the cutoff, negative outside, and
    steep in proportion to how finely the margin says the distance was resolved.

    A null ``d_heavy`` is a **confident non-contact**: :func:`tcren.potts.available_pairs` measures
    the heavy-atom layer only out to ``probe``, so a null means the pair missed by more than that,
    which is the reading — not missing data. Pass ``probe = cutoff + 3 * margin`` and the boundary
    sits where the field is already below −3.

    Args:
        frame: Sites carrying ``d_heavy`` (and ``sigma`` when ``margin <= 0``).
        cutoff: Contact definition, Å — the model's own ``cutoff``.
        margin: Resolution of the distance call, Å. ``0`` is the hard crystal map.
        clip: Saturation of the returned log-odds.

    Returns:
        One float per row of ``frame``.

    Example:
        >>> import numpy as np, polars as pl
        >>> f = pl.DataFrame({"d_heavy": [4.0, 5.0, 6.0, None]})
        >>> np.round(observation_field(f, margin=0.5), 3)
        array([  2.,   0.,  -2., -40.])
    """
    d = frame["d_heavy"].to_numpy().astype(float)
    absent = ~np.isfinite(d)
    if margin <= 0.0:
        # The limit, not a shortcut: an infinitely sharp read pins every site, so the field is
        # +-clip on the observed indicator and the posterior below is the crystal map itself.
        return np.where(absent | (d > cutoff), -clip, clip)
    return np.where(absent, -clip, np.clip((cutoff - d) / margin, -clip, clip))


def contact_posterior(sites: pl.DataFrame, model: PottsModel, *, margin: float = 0.5,
                      clip: float = LOGIT_CLIP, chains: int = 64, burn: int = 100,
                      draws: int = 100, thin: int = 3, seed: int = 0,
                      workers: int | None = None) -> pl.DataFrame:
    r"""The E-step: :math:`P(\sigma_a = 1 \mid \text{the crystal's distances}, \theta)` per site.

    :func:`tcren.potts.contact_probabilities` computed against the observation-shifted field, so
    the couplings still propagate a resampled site to its neighbours while the confidently read
    sites hold their value. Same sampler, same worker split, same
    ``blake2b(seed, pdb.id)`` streams — so the result is a function of the structure and not of its
    position in the table.

    At ``margin <= 0`` no sampling happens at all: the posterior is a point mass on the observed
    map, and ``p_post`` is ``sigma`` exactly. That exactness is what makes the reduction to TCRen2
    a test rather than an approximation.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`, ideally with
            ``probe = cutoff + 3 * margin``.
        model: The current :class:`PottsModel`.
        margin: Resolution of the distance call, Å; see :func:`observation_field`.
        clip: Saturation of the observation log-odds.
        chains, burn, draws, thin: Gibbs settings; see :func:`tcren.potts.gibbs`.
        seed: Seed for the sampler.
        workers: Processes to split structures across. ``None`` uses every core.

    Returns:
        The prepared site frame with ``eta`` (post-observation), ``eta_prior`` (the model's own)
        and ``p_post`` appended, one row per site.
    """
    upid, starts, e, sigma, E, eoff, kv, q = _prepare(sites, model)
    field = observation_field(q, cutoff=model.cutoff, margin=margin, clip=clip)
    out = q.drop("sid", "loop", "pchain").with_columns(
        pl.Series("eta_prior", e), pl.Series("eta", e + field))
    if margin <= 0.0:
        return out.with_columns(pl.col("sigma").cast(pl.Float64).alias("p_post"))
    rows = _map_structures("contact_probabilities", len(upid),
                           _payload(upid, starts, e + field, sigma, E, eoff, kv),
                           dict(chains=chains, burn=burn, draws=draws, thin=thin, seed=seed),
                           workers)
    p = np.zeros(len(e))
    for r in rows:
        p[r["lo"]:r["hi"]] = r["p_model"]
    return out.with_columns(pl.Series("p_post", p))


def weighted_contacts(post: pl.DataFrame, *, weight_col: str = "p_post") -> pl.DataFrame:
    """Rename :func:`contact_posterior`'s output into the schema :func:`derive_tcren` reads.

    The two halves of this loop were written years apart and name the same thing differently —
    ``aa.rec``/``aa.par`` against ``residue.aa.from``/``residue.aa.to``. This is the whole bridge.

    Every available pair is emitted, not only the contacts: a site at ``p_post`` near zero
    contributes near-zero mass and dropping it would be the same tally with more code. That the
    universe is the 15 Å available set rather than the observed contacts is exactly what lets a
    near-miss vote.

    Args:
        post: Output of :func:`contact_posterior`.
        weight_col: Which probability column becomes the fractional count.

    Returns:
        ``pdb.id``, ``residue.aa.from``, ``residue.aa.to`` and the weight column.
    """
    return post.select("pdb.id",
                       pl.col("aa.rec").alias("residue.aa.from"),
                       pl.col("aa.par").alias("residue.aa.to"),
                       pl.col(weight_col))


def _dense(pot: Potential) -> pl.DataFrame:
    return pot.matrix.select("residue.aa.from", "residue.aa.to", "value")


def _blend(old: Potential, new: Potential, damping: float) -> tuple[Potential, float]:
    """``(1 - damping) * old + damping * new``, and the largest cell movement it made."""
    j = (_dense(old).rename({"value": "v_old"})
         .join(_dense(new).rename({"value": "v_new"}),
               on=["residue.aa.from", "residue.aa.to"], how="inner"))
    if j.height != old.matrix.height:
        raise ValueError(f"the derived matrix covers {j.height} of the reference's "
                         f"{old.matrix.height} cells; the two derivations disagree on the alphabet")
    a, b = j["v_old"].to_numpy(), j["v_new"].to_numpy()
    mixed = (1.0 - damping) * a + damping * b
    return (Potential(name=new.name,
                      matrix=j.select("residue.aa.from", "residue.aa.to",
                                      pl.Series("value", mixed)),
                      alphabet=new.alphabet),
            float(np.max(np.abs(mixed - a))))


def fit_selfconsistent(sites: pl.DataFrame, *, potential: Potential | None = None,
                       margin: float = 0.5, damping: float = 0.5, steps: int = 15,
                       tol: float = 5e-3, pseudocount: int = 1, pin_centred: bool = False,
                       weights: dict[str, float] | None = None, ridge: float = DEFAULT_RIDGE,
                       chains: int = 64, burn: int = 100, draws: int = 100, thin: int = 3,
                       seed: int = 0, workers: int | None = None, notes: str = "",
                       **fit_kwargs) -> tuple[Potential, PottsModel, pl.DataFrame]:
    r"""Iterate potential :math:`\leftrightarrow` contact model to a fixed point.

    Each step recomputes the contacts under the current model (E-step), re-derives the potential
    from those contacts weighted by their probability, damps it against the previous matrix, and
    refits the model with the new matrix pinned as its pair term.

    **Two things are estimated each step, not one.** The 20x20 energy matrix is *re-derived* from the
    probability-weighted counts; the coupling kernel — the 18 classes over sequence offsets that
    carry the "crosses", a made contact pulling its own sequence neighbours onto the same partner
    residue — is *refitted* by pseudolikelihood, jointly with every field and with the scale on the
    matrix. They are not independent: the E-step's posterior for a borderline site depends on
    whether its neighbours contacted, so the couplings shape the counts the matrix is derived from.
    Both land in the trace.

    **The pair term is pinned, not free.** ``J`` carries one coefficient — the scale
    ``beta_matrix`` on the current matrix — instead of 400. That is what closes the loop: a free
    ``J`` would absorb the pair structure and leave the derived matrix with no route back into
    :math:`P(\sigma)`. It also makes the fixed point well posed and ``beta_matrix`` reportable, as
    the inverse temperature relating TCRen units to contact log-odds.

    ``pin_centred=False`` by default here, against ``fit_potts``'s own default: double-centring
    strips the one-body terms, which is right for *ranking* matrices against each other and wrong
    for a loop whose output is a matrix meant to be scored raw.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs` with ``probe`` past the cutoff.
        potential: Starting matrix. ``None`` (default) is TCRen2, the shipped deterministic
            estimate — the loop's :math:`e^{(0)}`.
        margin: Resolution of the distance call, Å. ``0`` returns after one step having reproduced
            the starting derivation.
        damping: Mixing weight :math:`\lambda` on the newly derived matrix. Fixed-point iteration
            on a statistic of its own output overshoots undamped.
        steps: Maximum iterations.
        tol: Convergence on ``max_abs_de``, in the potential's own (natural-log) units. The
            default sits **above** the sampler's noise floor deliberately — see the note below.
        pseudocount: Passed to :func:`derive_tcren`; ``1`` is the shipped value, and dropping it
            sends 17 of 380 cells to infinity.
        pin_centred: Double-centre the pinned matrix. See above.
        weights: Per-structure redundancy weights, applied to **both** the derivation and the fit —
            e.g. :func:`tcren.potential.balanced_weights`, which is what ships TCRen2.
        ridge: Ridge on the pseudolikelihood fit.
        chains, burn, draws, thin: Gibbs settings for the E-step.
        seed: Base seed; step ``t`` samples at ``seed + t`` so the steps are not correlated draws.
        workers: Processes for the E-step.
        notes: Free text recorded on the returned model.
        **fit_kwargs: Forwarded to :func:`fit_potts` (``radius``, ``cutoff``, ``joint``, …).

    Returns:
        ``(potential, model, trace)`` — the converged matrix, the model pinned to it, and one trace
        row per step carrying ``step``, ``max_abs_de``, ``r_vs_initial``, ``beta_matrix``,
        ``contact_mass``, ``mass_from_uncertain``, ``pseudo_loglik``, ``n_sites``, ``converged``,
        and one column per coupling class of :func:`tcren.potts.kernel_names` — the loop
        re-estimates the matrix and the couplings together, so the trace carries both.

    Raises:
        ValueError: If ``steps`` is below 1, or a derived matrix's alphabet stops matching.

    Note:
        **``max_abs_de`` has a Monte-Carlo floor and ``tol`` has to sit above it**, the same way
        :func:`tcren.potts.fit_potts_em`'s ``max_shift`` does. The E-step is sampled, so once the
        matrix stops moving systematically the trace oscillates at the sampler's own noise rather
        than reaching zero. Measured on Native2026 at ``--balance both``, 359 structures, the
        default Gibbs settings (``chains=64``, ``draws=100``, ``thin=3``): ``max_abs_de`` falls
        0.3499 -> 0.0033 over ten steps and then sits between 0.0033 and 0.0038 for the next five.
        Hence a default ``tol`` of 5e-3. Lower it only together with more ``draws`` or ``chains``,
        which is what lowers the floor; lowering ``tol`` alone buys a run that cannot converge.
    """
    if steps < 1:
        raise ValueError(f"steps must be at least 1, got {steps}")
    pot = potential if potential is not None else tcren2()
    initial = _dense(pot)["value"].to_numpy().copy()
    model = fit_potts(sites, coupling_matrix=pot, pin_centred=pin_centred, weights=weights,
                      ridge=ridge, notes=notes, **fit_kwargs)
    rows = []
    for step in range(1, steps + 1):
        post = contact_posterior(sites, model, margin=margin, chains=chains, burn=burn,
                                 draws=draws, thin=thin, seed=seed + step, workers=workers)
        p = post["p_post"].to_numpy()
        derived = derive_tcren(weighted_contacts(post), weights=weights, weight_col="p_post",
                               pseudocount=pseudocount)
        pot, de = _blend(pot, derived, damping)
        cur = _dense(pot)["value"].to_numpy()
        model = fit_potts(sites, coupling_matrix=pot, pin_centred=pin_centred, weights=weights,
                          ridge=ridge, notes=notes, **fit_kwargs)
        mass = float(p.sum())
        row = {
            "step": step,
            "max_abs_de": de,
            "r_vs_initial": float(np.corrcoef(initial, cur)[0, 1]),
            "beta_matrix": float(model.beta_matrix),
            "contact_mass": mass,
            "mass_from_uncertain": (float(p[(p > 0.05) & (p < 0.95)].sum() / mass)
                                    if mass > 0 else float("nan")),
            "pseudo_loglik": model.pseudo_loglik,
            "n_sites": int(post.height),
            "converged": bool(de <= tol),
        }
        # The couplings ride in the trace beside the scale, because they are the OTHER half of
        # what each step re-estimates: the matrix says which pairs are favourable, the kernel says
        # how a made contact pulls its sequence neighbours in. Reading only `beta_matrix` would
        # hide whether the crosses moved while the matrix did.
        row.update(dict(zip(kernel_names(model.joint), (float(k) for k in model.kernel))))
        rows.append(row)
        if de <= tol:
            break
    pot = Potential(name="TCRen-SC", matrix=pot.matrix, alphabet=pot.alphabet)
    model.notes = (f"{notes} | self-consistent: margin={margin}, damping={damping}, "
                   f"steps={len(rows)}, converged={rows[-1]['converged']}").strip(" |")
    return pot, model, pl.DataFrame(rows)
