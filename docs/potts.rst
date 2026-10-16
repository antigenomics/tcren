Historical contact-map Potts model
==================================

These fitted models and their empirical outputs are retained for historical and
methodological reproduction. They are excluded from final fit-free scoring.
The local benchmark contact GLM is distinct from this coupled model.

Every scoring path elsewhere in ``tcren`` reads a contact map that a structure *has*. This one
models the map itself: which residue pairs **could** have contacted, which did, and how the
answers at neighbouring cells depend on each other.

A **site** :math:`a = (i, j)` is an available pair — a receptor residue :math:`i` and a partner
residue :math:`j` whose C\ :sub:`α` atoms lie within ``radius`` (15 Å by default) — and
:math:`\sigma_a = 1` iff a heavy-atom contact formed within ``cutoff`` (5 Å, the TCRen contact
definition, unchanged). A whole contact map is the configuration :math:`\sigma`, and the model is
a Boltzmann distribution over it:

.. math::

   E(\sigma) = -\sum_a \eta_a \sigma_a - \tfrac{1}{2}\sum_{a,b} A_{ab}\,\sigma_a \sigma_b,
   \qquad P(\sigma) = \frac{e^{-E(\sigma)}}{Z},
   \qquad Z = \sum_{\sigma \in \{0,1\}^{n}} e^{-E(\sigma)}

Here :math:`n` is the number of available sites, :math:`a,b` index sites,
:math:`\eta_a` is the contact log-odds field, and :math:`A_{ab}` is a symmetric
zero-diagonal coupling. :math:`E` is a dimensionless statistical energy;
:math:`\log Z` and entropy use natural-log units (nats). These are not measured
binding free energies.

Site tables retain physical ``chain.id.rec`` and ``chain.id.par`` plus heavy-atom
measurement ``probe`` (Å). Partner offsets and lengths belong to a chain and
region together, so equal class-II groove-region names cannot merge chains.
Orientation lookups index the same compact Cα rows used to compute their matrix.
Unavailable orientations remain NaN; a model fitted with orientation rejects
missing geometry explicitly instead of assigning it a bin.
Legacy peptide-only tables remain readable. Legacy MHC tables without partner
chain identity must be rebuilt from their structures. Frozen coefficient arrays
are unchanged; class-II predictions require fresh assessment with corrected sites.

The one-body term is additive over categorical blocks,

.. math::

   \eta_a = \alpha + h^{\mathrm{rec}}(a_i) + h^{\mathrm{par}}(b_j) + J(a_i, b_j)
            + g_{\mathrm{dist}} + g_{\mathrm{region}} + g_{\mathrm{role}} + g_{\mathrm{class}}

so the **fields** carry single-residue propensity — "the backbone put you in reach; did your side
chain engage?" — ``J`` carries the pair chemistry, and ``g`` the geometry and annotation the fields
must be adjusted for. Every block is in the zero-sum (Ising) gauge, which makes ``J`` directly
comparable with a double-centred :class:`~tcren.potential.Potential`.

Why a reference state of *available* pairs
------------------------------------------

A TCRen potential is a Boltzmann inversion **conditioned on a contact existing**, so a residue that
could have reached the peptide and declined contributes nothing at all. Here that non-event is the
observable, which is what lets the fields separate reach from chemistry.

Why the couplings
-----------------

Contacts on a real interface are not independent. On the 362 αβ ``Native2026`` crystals, at fixed
C\ :sub:`α` distance, the chance of a contact runs 0.155 → 0.401 when the same receptor residue's
neighbouring peptide residue contacts (odds ratio 5.37 over 48,921 pairs). Three coupling families
capture it, all defined on sequence offsets so no extra coordinates are needed:

``K(di, dj)``
   **within-loop**, :math:`|di| \le 2`, :math:`|dj| \le 2`, both sites in one receptor loop and one
   partner chain — 12 classes after the :math:`K(d) = K(-d)` symmetry;
``L(|dj|, same chain?)``
   **cross-loop**, the two sites in different hypervariable loops — 6 classes;
``M``
   **cross-class**, the same receptor residue against both the peptide and the MHC groove — 1
   class, joint models only;
``X(same receptor residue, |dj| > 2)``
   **long-offset exclusion**, the complement of ``K``'s reach on one receptor residue — 1 class,
   ``far=True`` only. New in 3.8.0.

Everything else is asserted uncoupled: the model's statement that linkage falls to zero across
loops and beyond two residues.

.. note::

   ``X`` exists because that statement is not neutral at long offset. ``K(+0,+1)`` and ``K(+0,+2)``
   are strongly positive — a receptor residue gripping peptide position *j* also grips *j+1* and
   *j+2* — and nothing in the model then says the run stops, though a side chain gripping one end of
   the groove is plainly not also at the other. Fitted on 367 crystals with ``J`` pinned to
   ``tcren2``, the class comes out at **K = −0.457, s.e. 0.049, z = −9.34**, worth **+81.8 nats** of
   pseudo-log-likelihood for one parameter.

   It is **off by default** and every shipped model carries ``far=False``, so nothing already
   written changes. Where it bites is a *predicted* Cα geometry: on the crystal's own coordinates
   the field is already sharp and the couplings have little left to move (Jensen–Shannon similarity
   0.7133 → 0.7162 over 367 crystals), while under a geometry averaged over reference crystals it
   improves both placement and level on three of the four ``Native2026`` receptors with held-out
   priors — js 0.4846 → 0.4957, 0.5226 → 0.5350, 0.4834 → 0.4924, with the expected contact count
   moving 36.80 → 34.11 against 32 observed, 33.81 → 31.40 against 33, and 29.97 → 27.82 against 25.
   On the fourth (``3qdg``) it is flat, and that receptor is already under-counted.

On the crystals every **axial** class comes out positive and every **off-axis** class negative — a
made contact recruits its own sequence neighbours onto the *same* partner residue and suppresses
the diagonal one. That sign pattern is not visible in the raw data; it appears only once the axial
terms and the distance profile are held fixed.

Fitting needs no partition function
-----------------------------------

The conditional :math:`P(\sigma_a = 1 \mid \sigma_{-a})` is logistic in
:math:`\eta_a + \sum_k K_k n_k(a)` with :math:`n_k(a)` the count of contacting neighbours in class
``k``, so the coupled fit is an ordinary weighted-binomial GLM with a handful of extra integer
covariates, and it stays concave. That is Besag's pseudolikelihood; consistency for this model
class is Ravikumar, Wainwright & Lafferty (`arXiv:1010.0311 <https://arxiv.org/abs/1010.0311>`_),
and plmDCA (`arXiv:1211.1281 <https://arxiv.org/abs/1211.1281>`_) is the same recipe on Potts
sequence variables. The design is over-parametrised and identified by an :math:`\ell_2` ridge, then
projected to the zero-sum gauge — **penalise then project**, because an :math:`\ell_2` penalty
otherwise picks its own gauge.

``PottsModel.n_parameters()`` counts stored fitted coefficients, including
orientation and redundant or unobserved one-body levels. It is not the number of
identifiable degrees of freedom. ``pseudo_loglik`` stores the unpenalized weighted
value at the penalized optimum. IRLS validates its numerical domain and backs off
Newton steps that decrease the penalized objective; failure to converge raises.

Scoring does, and it is exact at the reference
----------------------------------------------

``Z`` is estimated by annealed importance sampling (Neal, *Stat. Comput.* **11**:125–139, 2001),
annealing *only* the coupling term. At :math:`\beta = 0` the model is the uncoupled one, whose
partition function is exact and closed form,
:math:`\log Z_0 = \sum_a \log(1 + e^{\eta_a})` — so the reference is a verified model rather than
an approximation, and the estimator is unbiased in ``Z``. Transitions are block Gibbs on a greedy
colouring of the actual coupling graph; same-colour sites are conditionally independent, which is
asserted against the real edge lists rather than argued.

Always read ``ais_ess`` in the output. It is the effective sample size of the AIS weights out of
``--particles``; close to ``--particles`` means the annealing schedule was long enough, and a small
value means it was not.

Bound versus unbound, for the whole interface
---------------------------------------------

A single site has two states and :math:`\eta_a` is the free energy between them. The same contrast
for the whole interface needs a macrostate, and the contact count :math:`N(\sigma) = \sum_a \sigma_a`
defines one. Because :math:`E(\varnothing) = 0` exactly, the observed map's log-odds against the
empty configuration is just its negated energy, and the scores already emitted by
:func:`~tcren.potts.score_sites` decompose:

.. math::

   \underbrace{-E(\sigma^{\mathrm{obs}})}_{\text{binding log-odds}}
   \;=\; \underbrace{\log Z}_{\text{capacity}}
   \;+\; \underbrace{\mathcal{L}(\sigma^{\mathrm{obs}})}_{\text{typicality}}

:func:`~tcren.potts.bound_unbound` gives three readings of that contrast, all from one Gibbs pass,
because every tilt in :math:`N` is an exponential family:

``df_empty``
    :math:`\log[P(N \ge 1)/P(N = 0)] = \log(Z - 1)`, exact, from the AIS ``log Z``.
``df_threshold``
    :math:`\log[P(N \ge x)/P(N < x)]` from the sampled histogram — ``Z`` cancels, so no AIS is
    needed, but it is only finite while ``x`` sits inside the sampled range.
``mu_star``
    the chemical potential at which :math:`\langle N\rangle_\mu` equals the observed count. Positive
    means the structure made more contacts than its fields and couplings warrant. ``nan`` outside
    the sampled support, where reweighting would be extrapolating.

The two are not competing estimates of one number. The unbound basin of a *docked* pose is
astronomically improbable — the model is conditioned on an available set that already holds the
receptor against the peptide — so no sampler reaches :math:`N = 0`, and only the ``log Z`` route
gives it. :func:`~tcren.potts.count_profile` returns the pooled :math:`F(N) = -\log p(N)` profile so
the landscape can be inspected before a threshold is chosen.

Constraining a statistic of the whole configuration
---------------------------------------------------

Jaynes' construction lifts a constraint from a single site to any statistic of the configuration:
fix :math:`\langle O_k\rangle` and its multiplier enters the Hamiltonian as :math:`-\lambda_k O_k(\sigma)`.
Tkačik *et al.* (*PLoS Comput. Biol.* **10**\ (1):e1003408, 2014) is this with ``O`` the total
activity — the "K-pairwise" model — and it is why :func:`~tcren.potts.gibbs` takes an ``observer``
callback: a statistic of whole configurations can be accumulated during sampling without
materialising every draw.

Because such a statistic depends on :math:`\sigma` only through a low-dimensional summary, the
tilted expectation is an importance-weighted average over draws taken at :math:`\lambda = 0`, so one
sampling pass serves every step of a moment-matching fit. :func:`~tcren.potts.tilt_mean` and
:func:`~tcren.potts.mu_star` implement that reweighting for the contact count. A **linear** tilt is
exactly a constant added to every field, :math:`E - \mu N = -(\eta + \mu)\cdot\sigma`, which is what
makes the reweighting identity checkable against direct simulation rather than assumed.

Shipped models
--------------

.. list-table::
   :header-rows: 1
   :widths: 22 12 12 54

   * - Key
     - Sites
     - Contacts
     - Provenance
   * - ``potts_tcr_peptide``
     - 64,622
     - 7,865
     - **The default.** TCR:peptide over the 362 αβ ``Native2026`` crystals, redundancy-balanced
       on both axes. Reproduce with
       ``tcren potts fit -s data/Native2026 -o potts_tcr_peptide.json --balance both``.
   * - ``potts_tcr_peptide_sc``
     - 64,064
     - 7,799
     - The same interface with the pair term **pinned on** :func:`tcren.potential.tcren2_em` rather
       than on TCRen2, over the 359 crystals that survive the withdrawn-entry rule;
       ``beta_matrix`` +0.7629. Reproduce with ``tcren potts selfcon`` (see :doc:`potentials`) —
       it is the model that loop leaves behind, not a separate fit.
   * - ``potts_tcr_mhc``
     - 239,093
     - 15,451
     - TCR:MHC groove over the same crystals — twice as many contacts as the peptide interface.
       Reproduce with ``tcren potts fit -s data/Native2026 --partner mhc --balance both -o …``
       (needs the allele reference: ``tcren build-mhc-ref``).
   * - ``potts_pep_mhc``
     - 195,674
     - 23,492
     - **The presentation arm**: the groove's grip on the peptide, over the same crystals. Here the
       groove is the *receptor*, so ``region.rec`` is a groove region and the model carries
       :data:`tcren.potts.MHC_RECEPTOR_REGIONS` rather than the TCR loop set. Reproduce with
       ``tcren potts fit -s data/Native2026 --receptor mhc --partner peptide --balance both -o …``.

Load any of them with :meth:`tcren.potts.PottsModel.bundled`, or by name on the command line —
``-m/--model`` on ``potts score``, ``contacts``, ``map``, ``ensemble`` and ``scan`` takes a
bundled name as well as a JSON path::

   tcren potts contacts -s complex.pdb.gz -m potts_tcr_peptide_sc -o contacts.tsv


The three arms are three separate Hamiltonians, not one field read three ways, and the fit says so:
the peptide:MHC coupling block correlates with the shipped TCRen2 at **r = +0.017** and with
Miyazawa-Jernigan at **r = +0.042**. The groove's chemistry is its own. Its per-position contact
profile is the textbook one — on 1ao7, peaks of 10, 9 and 10 contacts at P1, P2 and P9 against 1 at
the TCR-facing bulge.

Which potential belongs on which interface
------------------------------------------

Fixing ``J`` to one scale on a bundled potential (``--coupling-matrix``) gives every candidate an
identical parameter count and an identical design, so their pseudo-log-likelihoods compare
directly. Because the matrix is double-centred it contributes nothing to the one-body marginals —
the fields carry those, refitted freely — so the comparison is about pair structure alone.

On the ``Native2026`` crystals the ranking **inverts** between the two interfaces: TCRen2 beats MJ
by 103.3 nats on TCR:peptide, and MJ beats TCRen2 by 35.5 nats on the TCR:MHC groove. TCRen2's
fitted scale falls 5.4-fold across that move (:math:`\beta` = +1.131 → +0.209) while MJ's barely
changes (+0.803 → +0.974). This is the measurement behind ``tcren``'s long-standing default of
scoring ``Phi_tcr_mhc`` with Miyazawa–Jernigan and reserving TCRen for TCR:peptide.

Three limits of one free energy
-------------------------------

The interface free energy is closed form, :math:`\Phi = \sum_a \operatorname{softplus}(\eta_a)`, and
its derivative is the contact marginal, :math:`\partial\Phi/\partial\eta_a = p_a`. So :math:`\Phi`
is an interaction sum weighted by how often each pair is *in* contact, and the fixed contact map
every other scoring path in ``tcren`` reads is the special case :math:`p_a \in \{0, 1\}`. Three
limits follow, and they are not competing models — each is the one the interface's own physics
selects.

.. list-table::
   :header-rows: 1
   :widths: 18 34 48

   * - limit
     - weight on each pair
     - the interface it fits
   * - hard contact
     - :math:`\mathbf{1}[d \le 5\,\text{Å}]`
     - reproduces :func:`tcren.ddg.reference_delta` exactly, when ``J`` is pinned with
       ``pin_centred=False``
   * - smoothed
     - :math:`p_a \in (0, 1)`
     - a *plastic* interface, where an ensemble average is a different number from one snapshot
   * - saturated
     - :math:`p_a \to \{0, 1\}`
     - a *frozen* interface, where the two above coincide and a fixed map is already exact

Which limit an interface is in is measurable, not a matter of taste. Over four deposited 100 ns
trajectories, collapsed onto the peptide and restricted to the pairs inside the availability radius,
**36 of the 38 engaged peptide positions (94.7%) reach a maximum groove-pair contact frequency above
0.98, against 14 of 38 (36.8%) on the receptor side** — at the groove a peptide residue is in
contact in essentially every frame or not at all, while at the receptor loops the same maximum runs
the whole range from 0.00 to 1.00. Among pairs that touch at all the groove mean is 0.468 against
0.322, and 28.6% of groove pairs sit above 0.98 against 17.4% of receptor pairs.

That is why the groove takes a fixed map and Miyazawa–Jernigan: it is already in the saturated
limit, so smoothing it can only add noise. The receptor loops are not, which is what the smoothed
limit is for.

What the marginals settle, and what they do not
-----------------------------------------------

The interaction energy of a contact map is linear in the map,
:math:`\Phi(\sigma) = \sum_a \sigma_a \varphi_a`, so by linearity of expectation

.. math:: \langle\Phi\rangle \;=\; \sum_a \langle\sigma_a\rangle\, \varphi_a

for **any** joint distribution over :math:`\sigma`, coupled or not. Dependence between sites does
not break that identity. What it does is fix its two conditions, and both are easy to violate:

1. :math:`\langle\sigma_a\rangle` must be the marginal of the **coupled** model — ``p_model``, by
   block Gibbs — and not ``p_independent``, the factorised :math:`\sigma(\eta_a)`.
2. It must be read **at the sequence being scored**. ``ddg(weights=...)`` and
   ``cpl.predict_response(tcr_weights=...)`` take one weight vector computed on the template and
   apply it to every threaded candidate, giving
   :math:`\sum_a p_a(\text{template})\,\varphi_a(x)`. A candidate residue enters :math:`\eta_a`
   twice, through :math:`h^{\mathrm{par}}` and through :math:`J`, so the contact ensemble itself
   moves under substitution and the two are different numbers.

:func:`tcren.potts.threaded_ensemble` recomputes the marginals per candidate.
It accepts peptide-partner sites only and validates each whole candidate against
``par.len``, including unrepresented tail residues. Later candidates start from
the previous retained matrix of chain configurations and use ``candidate_burn``;
candidate ordering remains part of the sampling protocol.
:func:`tcren.potts.energy_ensemble` carries the read-outs no set of marginals can determine:

.. list-table::
   :header-rows: 1
   :widths: 20 22 58

   * - read-out
     - units
     - what it is
   * - ``phi_var``
     - (potential units)²
     - :math:`\operatorname{Var}\Phi`, which needs the two-point function
       :math:`\langle\sigma_a\sigma_b\rangle`. How much the energy depends on *which* of the maps
       the interface admits is the one realised.
   * - ``phi_free``
     - potential units
     - :math:`-\beta^{-1}\log\langle e^{-\beta\Phi}\rangle`, a soft minimum over configurations.
       :math:`\beta \to 0` recovers :math:`\langle\Phi\rangle`; :math:`\beta \to \infty` the best
       map the geometry admits.
   * - ``map_entropy``
     - nats
     - :math:`\langle E\rangle + \log Z`, how many contact maps the interface is choosing between.
   * - ``hamming_obs``
     - count
     - mean distance from the map the structure itself made.

On the 14 αβ complexes in ``tests/assets/pdb`` the spread is not a rounding effect:
:math:`\operatorname{Var}\Phi` runs **0.341 to 1.286 squared TCRen units**, a 3.8-fold range, so
how plastic an interface is varies as much between complexes as its mean energy does.

A configuration needs coarse-graining before it has a probability
-----------------------------------------------------------------

:func:`tcren.potts.map_ensemble` returns the contact configurations themselves, ranked, each with
its probability. The granularity is not a free choice, and it is measured rather than assumed. On
1ao7's TCR:peptide interface — 174 available sites, 6,400 draws — the raw site-level map recurs
**never**: 6,400 distinct configurations out of 6,400 draws, every one at
:math:`\hat p = 1.6\times10^{-4}`, because the state space is :math:`2^{174}`. Coarse-grained to
which **peptide positions** are engaged at all, the same draws collapse to **41 distinct patterns
with the most likely at** :math:`\hat p = 0.3219` — and that modal pattern is the one the crystal
itself made.

.. list-table::
   :header-rows: 1
   :widths: 16 20 28 36

   * - ``by``
     - groups on 1ao7
     - distinct patterns / 6,400
     - reading
   * - ``"pair"``
     - 174
     - 6,400, top :math:`\hat p` 0.0002
     - unresolvable; use :func:`tcren.potts.exact_log_z` on an interface small enough to enumerate
   * - ``"loop"``
     - 39
     - 3,308, top :math:`\hat p` 0.0216
     - the granularity the molecular-dynamics validation works at; at the edge of resolvable
   * - ``"position"``
     - 9
     - 41, top :math:`\hat p` 0.3219
     - which peptide positions the receptor engages — the default

Across those same 14 complexes the modal engagement pattern is the observed one in **6 of 14**,
and ``n_distinct_total`` runs 37 to 556.

When the contact indicator is not data
--------------------------------------

:func:`tcren.potts.fit_potts` treats :math:`\sigma` as observed. In a crystal it mostly is — but a
contact is called by thresholding the closest heavy-atom distance at 5 Å, and a pair sitting within
half an ångström of that line is a coin toss recorded as a fact. :mod:`tcren.potts.em` treats those
sites as latent and runs a Gibbs-imputation/pseudolikelihood hybrid around the
existing fit. The full-posterior imputation and conditional-pseudolikelihood update
do not share the ordinary likelihood-EM monotonicity theorem.

**E-step.** Sample :math:`\sigma_{\mathcal L}` from :math:`P(\sigma_{\mathcal L} \mid
\sigma_{-\mathcal L})` by block Gibbs with every confident site **clamped** to what the structure
did. The couplings are what make this worth doing: a borderline pair beside three made contacts is
imputed differently from one standing alone.

**M-step.** The existing penalised pseudolikelihood on the completed data, ``draws`` completions
stacked at weight :math:`1/\text{draws}` so the effective sample size — and therefore what the ridge
means — is unchanged. Each completion carries a ``#m`` suffix on its ``pdb.id``, because every
coupling edge in :mod:`tcren.potts.kernel` joins on that column and unsuffixed copies would link
sites across draws that never coexisted.

.. warning::

   **The latent band has to reach both ways.** ``available_pairs`` measures ``d_heavy`` out to
   a requested ``probe=cutoff``, so every pair that missed carries a null distance.
   The current default probe is the availability radius. A band built on that table
   is one-sided: on the 17 bundled test structures it selects **44 sites, every one of them a
   contact**, and an EM run on it could remove contacts and never add them. Build the table with
   ``available_pairs(..., probe=cutoff + margin)`` — which populates ``d_heavy`` past the cutoff
   without moving ``sigma`` — and the same band selects **115 sites, 44 contacts and 71 near
   misses**. :func:`tcren.potts.latent_mask` checks the recorded ``probe`` against
   ``cutoff + margin``; an all-contact band alone does not establish truncation.
   Legacy tables without probe provenance need regeneration for publication.
   ``margin=0`` explicitly makes every site observed.

This **adds** a fit; it does not replace one. The shipped pseudolikelihood models are unchanged and
are the point of comparison. :func:`tcren.potts.consistency` reports how far a refit moved per
parameter block, scaled by the reference fit's own cluster-robust standard errors where it has them
(the coupling kernel) and by the block's own standard deviation where it does not. A block with no
spread — ``g_class`` on a single-partner model — is reported as ``degenerate`` rather than divided
by something near zero.

.. note::

   Monte-Carlo EM with a fixed ``draws`` does not converge to a point. It random-walks in a ball
   around the optimum whose radius is set by the Monte-Carlo noise, so ``max_shift`` has a floor
   that ``tol`` may sit below. Read ``max_shift_vs_init`` beside it: a run that is merely jittering
   has a flat one. Raise ``draws`` to lower the floor. Every sampler is seeded, so a run is exactly
   reproducible either way.

.. code-block:: text

   tcren potts fit -s <structures> -o potts_em.json --em --em-draws 5 --em-steps 10 \
       --em-trace-out em_trace.tsv

Retired geometry experiments
----------------------------

The geometry prior, coordinate-chain model and cross-product experiment are archived
under ``attic/2026-10-audit`` in the source repository. They are not shipped APIs.
The contact model on supplied coordinates and its fitting recipes remain supported.

Scoring the agreement: placement and level
------------------------------------------

:func:`tcren.potts.map_agreement` reads any scored site table --- from
:func:`contact_probabilities` or ``tcren potts contacts`` --- and reports **two** numbers per map, because no single one answers
both questions:

``ratio``
    The **level**, :math:`\sum_a p_a / \sum_a \sigma_a`. One is the contact count reproduced.

``js``
    The **placement**, Jensen--Shannon similarity of the two maps read as distributions over sites,
    :math:`q_a = \sigma_a/\sum\sigma` against :math:`r_a = p_a/\sum p`, in bits:

    .. math::  \mathrm{js} = 1 - \tfrac12 D_{\mathrm{KL}}(q\|m) - \tfrac12 D_{\mathrm{KL}}(r\|m),
               \qquad m = \tfrac12(q+r)

    Bounded in :math:`[0, 1]`, symmetric, and **scale-free in the level** --- normalising divides
    the total out, which is why ``ratio`` sits beside it rather than being folded in.

.. code-block:: console

   tcren potts agreement -i predicted.tsv
   tcren potts agreement -i contacts.tsv --column p_independent -o agreement.tsv

**Read ``js`` against the same model on the structure's own geometry, not against 1.**
:math:`\mathrm{js} = 1` requires the model to put all its mass on exactly the contacting cells in
equal proportion. A probabilistic map :math:`p=c\sigma`, :math:`0<c<1`, also
reaches 1 after normalization; ``ratio`` records its different level. On the four receptors of the manuscript's Figure 2c the crystal's own
C\ :math:`\alpha` geometry reads 0.75--0.81 and the predicted geometry 0.42--0.58, so the first is
the ceiling the second is measured against.

**R² is deliberately not offered.** On a binary outcome ``1 - SS_res/SS_tot`` *is* the Brier skill
score --- the two agree to three decimals on every map checked here --- so reporting both is
reporting one number twice. It is also dominated by the many easy zeros and scales with prevalence,
which makes two maps of different site counts incomparable under it; that is exactly the comparison
the predicted and measured geometries require, since the predicted site set is roughly twice the
size.

Reproducing a referenced score
------------------------------

``centred_potential`` double-centres by default, which is right for *ranking potentials against
each other* — with the one-body terms removed the comparison is about pair structure alone. It is
wrong for *reproducing a referenced contact-map score*. The peptide-referenced energy is a
difference of one-body sums, so double-centring re-injects a burial-scaled composition term
:math:`n_i\,c(a)`: the position's contact count times the potential's partner-residue column mean.
On TCRen2 that mean has s.d. 0.0668 and runs −0.212 to +0.027 over the 19 residues the potential
observes, and :math:`n_i` runs 1 to 54.

Pass ``pin_centred=False`` (or ``centre=False`` to ``centred_potential``) and the coupling *is* the
potential, so any linear read-out of the field reduces to the potential's own score up to the
fitted scale :math:`\beta`.

Command line
------------

.. code-block:: bash

   # fit a model (the alpha-beta TCR:pMHC HARD RULE applies, as in derive-potential)
   tcren potts fit -s structures/ -o potts.json --balance both

   # energy, log Z and the likelihood of each structure's observed contact map
   tcren potts score -s structures/ -o scores.tsv          # bundled model by default

   # per-residue-pair contact probability
   tcren potts contacts -s complex.pdb -o contacts.tsv

   # close those onto the grids an experiment measures
   tcren potts map -s complex.pdb --by loop     -o map.tsv    # CDR loop x peptide position
   tcren potts map -s complex.pdb --by position -o import.tsv # peptide residue importance

   # historical log-Z0 substitution effects (nats; not calibrated binding energy)
   tcren potts scan -s complex.pdb -o scan.tsv

``score`` emits one row per structure: ``n_sites`` and ``n_contacts`` (the available pairs and
how many of them engaged), ``energy`` and ``neg_energy`` (the Hamiltonian of the observed map and
its negation, lower energy being more favourable), ``log_z`` and ``log_z0``, ``log_lik`` and ``psi`` (the log-likelihood, and the
same per available pair so it compares across interfaces of different size),
``pseudo_log_lik``/``psi_pseudo`` as the MCMC-free cross-check, and ``ais_ess``.

``contacts`` emits one row per site with three probabilities, whose differences *are* the
couplings: ``p_independent`` from the one-body model alone, ``p_model`` the marginal of the full
coupled model by block Gibbs (the one to use), and ``p_conditional`` =
:math:`P(\sigma_a = 1 \mid \text{the observed rest})`.

``map`` groups pairs onto a loop/position grid. Its ``p_any`` uses the independence approximation

.. math:: P_{\mathrm{ind}}(N \ge 1) = 1 - \prod_j (1-p_j),

where :math:`p_j` is the coupled marginal of pair :math:`j`. Distinct pairs need not be independent:
coupled marginals alone do not determine the probability of their union. This approximation can
be compared empirically with a trajectory contact-frequency map, but is not its exact coupled
model expectation. No numerical definition changed in the audit.

``p_expected`` is :math:`\sum_j p_j`, the expected number of contacts; linearity makes this valid
without independence, and it need not lie in [0, 1]. Other columns are ``n_pairs``, ``n_observed``
and ``observed``. ``p_any`` lies in [0, 1]; the product is accumulated in log space. Neither this
probability approximation nor the expected count is a binding free energy.

``score``, ``contacts`` and ``map`` all take ``--workers`` (default: every core). The per-structure
numbers are functions of ``(seed, pdb.id)`` alone, so splitting the work changes nothing.

Historical substitution log-partition and linear response
----------------------------------------------------------

``map`` reads how engaged a position is expected to be *before any residue identity is scored*.
:func:`~tcren.potts.peptide_free_energy` reads what happens when the identity changes. The partner
residue enters :math:`\eta` in exactly two places, so for a site :math:`s` and a candidate
residue :math:`a`

.. math:: \eta_s(a) = r_s + h^{\mathrm{par}}(a) + J(\mathrm{rec}_s, a)

with :math:`r_s` — the intercept, the receptor field, the distance profile, region, role and class
— independent of :math:`a`. Threading :math:`a` through partner position :math:`i` moves every site
carrying that position, and the uncoupled log-partition changes with it:

.. math::

   \Phi^{\mathrm{Potts}}(x) = \log Z_0\big(\eta(x)\big) = \sum_s \log\!\big(1 + e^{\eta_s(x)}\big),
   \qquad
   \Delta Phi_i(a) = \Phi^{\mathrm{Potts}}(x_{i \to a}) - \tfrac{1}{20}\sum_b \Phi^{\mathrm{Potts}}(x_{i \to b})

Higher is more favourable, and the reference is the **equimolar** one — the mean over the twenty
residues at that position rather than the residue the structure carries — which is the null a
positional-scanning library holds its other positions at. ``coupled=True`` takes the linear
response about the observed sequence, :math:`\Delta \log Z \approx \sum_s p_s\,\Delta\eta_s`, since
:math:`\partial \log Z/\partial \eta_s = \langle\sigma_s\rangle`: one Gibbs pass, then a dot
product per cell. Because :math:`\log Z_0` is a sum over independent sites the result is additive
over positions, so one :math:`L \times 20` table scores a single substitution and any whole partner
sequence alike. Only ``aa.par`` changes — the backbone, the Cα distances, the receptor residues and
the partner roles stay the structure's own, the same fixed-backbone approximation every threading
score in the package makes. Custom radii/cutoffs must also agree with the model's fitted categorical blocks; changing the radius alone does not create missing field levels. A partner position carrying two different residues has no sequence to
substitute into and is rejected rather than averaged.

The positive :math:`\log Z_0` is dimensionless and measured in nats. The uncoupled
statistical free energy is :math:`-\log Z_0`; neither is a calibrated binding free
energy. With ``coupled=True``, the preserved ``log_z0`` column contains the frozen-
marginal linear surrogate, not coupled :math:`\log Z`. Finite substitutions do not
make the first-order approximation exact. Existing names and values are unchanged.

Only the full-field centred ``F_full`` quantity overlaps the uncoupled scan.
The centred ``F_full_p`` quantity overlaps its fixed-native-marginal path when
identical sites, parameters and marginal values are supplied. Benchmark
composition-removed, pair-only, observed-support, normalized-contact and
within-clone standardized combinations have no complete equivalent CLI route.
Do not replace them with ``potts scan``. Their equations, original results and
losses remain historical evidence; no new fitted replacement is authorized.
See the issue-32 disposition record in ``audit/2026-10-16/``.

The two terms it separates are different quantities, and which one a task needs is an empirical
question. :math:`h^{\mathrm{par}}` is **composition** — how much a residue engages an available
partner at all, wherever it is put — while :math:`J` is **complementarity**, the pair chemistry, and
is the block :math:`\beta_\Phi` ties to TCRen2 above.

.. math:: \text{eta} \;=\; \underbrace{h^{\mathrm{rec}} + h^{\mathrm{par}}}_{\text{composition}}
          \;+\; \underbrace{J}_{\text{complementarity}} \;+\; \underbrace{g}_{\text{geometry}}

Closing the loop: the potential re-derived from the contacts the model predicts
-------------------------------------------------------------------------------

Everything above reads a **finished** potential. :func:`tcren.derive_tcren` tallies the amino-acid
pairs of the contacts a crystal *shows*, and :func:`tcren.potts.fit_potts` reads the resulting
matrix back as the pair term in its field. Nothing carries the model's own view of which contacts
formed into the tally. :mod:`tcren.potts.selfcon` closes that loop.

.. math::

    \eta_a^{\mathrm{post}} &= \eta_a(\theta^{(t)})
        + \operatorname{clip}\!\left(\frac{c - d_a}{m},\ \pm K\right) \\
    p_a^{(t)} &= \langle \sigma_a \rangle_{P(\sigma \mid \mathrm{obs},\, \theta^{(t)})} \\
    N^{(t)}(a,b) &= \sum_s w_s \sum_{a\,:\,\mathrm{pair}(a)=(a,b)} p_a^{(t)} \\
    e^{(t+1)}(a,b) &= -\ln \frac{N^{(t)}(a,b)\,N_{\mathrm{tot}}}
                                {N_{\mathrm{row}}(a)\,N_{\mathrm{col}}(b)}, \qquad
    \theta^{(t+1)} = \arg\max\ \mathrm{PL}\big(\theta \mid e^{(t+1)}\big)

with :math:`c` the 5 Å cutoff, :math:`d_a` the site's closest-heavy-atom distance, :math:`m` the
margin, :math:`K` the logit clip and :math:`w_s` the structure's redundancy weight.

**Weighting each conformation by its probability costs nothing.** The tally over whole contact
configurations, :math:`\sum_\delta P(\delta \mid \mathrm{obs})\, n_{ab}(\delta)`, **equals**
:math:`\sum_a p_a` by linearity of expectation — the same identity that governs
:func:`tcren.potts.energy_ensemble`. So the configurations never have to be enumerated to get the
counts; :func:`tcren.potts.map_ensemble` lists them as a read-out, not as a dependency.

**The estimator does not change; only the counts it is fed.** :func:`tcren.derive_tcren` already
accepts a per-contact fractional weight, documented there as the way to let a pair that "will only
sample an interacting geometry some of the time" cast a fractional vote. The log-odds, the
contact-fraction reference state and the pseudocount are the shipped ones, untouched.

**The observation enters as one additive field.** The distance sigmoid
:math:`\operatorname{sigmoid}((c-d_a)/m)` defines a multiplicative observation field.
A Bayesian posterior interpretation requires its logit to be the observation
likelihood ratio, or the sigmoid to use an equal-odds reference prior; an arbitrary
posterior probability given distance does not suffice. The field adds a constant to
:math:`\eta`, so the E-step is the existing block-Gibbs sampler called with a modified vector — no
clamping machinery, no second code path. One knob spans the family: ``margin`` → 0 pins the crystal
map, ``margin`` = 0.5 Å releases the threshold band alone, and a large ``margin`` lets the crystal
enter only through the fit.

**The pair term is pinned, not free.** ``J`` carries one coefficient — ``beta_matrix``, the scale on
the current matrix — rather than 400 cells. That is what closes the loop: a free ``J`` would absorb
the pair structure and leave the re-derived matrix no route back into :math:`P(\sigma)`. It also
makes ``beta_matrix`` reportable as a statistical scale relating TCRen
natural-log units to contact log-odds. It is not a calibrated physical temperature. :func:`tcren.potts.fit_potts` accepts a :class:`~tcren.potential.Potential` object for
``coupling_matrix``, not only one of the five bundled names, and a model pinned that way stores the
matrix itself so it round-trips.

.. warning::

   **The heavy-atom layer has to reach past the cutoff**, exactly as for
   :mod:`tcren.potts.em`. Build the table with ``available_pairs(..., probe=cutoff + 3*margin)``, or
   ``d_heavy`` is null for every non-contact and the observation field can only ever *remove*
   contacts. ``tcren potts pairs --probe`` and ``tcren potts selfcon`` both do this.

.. note::

   **This is EM-shaped but it is not EM.** The E-step is a genuine conditional expectation, but the
   potential update is the shipped log-odds re-derivation applied to expected counts, not the
   maximiser of any one objective — so the pseudo-log-likelihood need not rise, and on the 17
   bundled crystals it does not (:math:`-337.94` then :math:`-335.05`, drifting to :math:`-335.74`
   over six steps at ``damping=0.5``). Convergence is declared on ``max_abs_de``, the largest cell
   movement in natural-log units, which contracts at the damping rate. ``pseudo_loglik`` rides along
   as a diagnostic.

**The reduction gate.** At ``margin = 0`` the crystal's distance read pins every site, the posterior
is the observed map exactly, and the probability-weighted tally reproduces the plain unweighted one
**bit for bit** — on the bundled crystals, ``max |soft − hard| = 0`` over all 380 cells and a
posterior mass of 309.000000 against 309 observed contacts. ``contact_mass`` and
``mass_from_uncertain`` sum unweighted site probabilities; structure weights enter
the fit and potential derivation separately. That is the check that the two halves of
the loop are wired to the same contacts, and it is a test rather than an approximation.

.. note::

   That gate is exact **on one site table**, and the Potts availability set is not quite the
   derivation's contact set. Over Native2026 the 15 Å Cα radius holds **7,799 of the deposited
   derivation's 7,816 TCR:peptide contacts (99.78 %)**; the 17 it misses are spread over 17
   structures at one or two each, side chains touching while their Cα sit further apart than the
   radius. So ``margin = 0`` reproduces the derivation *from the availability set*, not the shipped
   ``TCRen2_potential.csv`` — which is derived over the full contact map.

.. code-block:: text

   tcren potts selfcon -s <structures> --balance both --margin 0.5 \
       -o TCRen_sc.csv --model-out potts_sc.json --trace-out sc_trace.tsv

API
---

.. automodule:: tcren.potts
   :members:
   :undoc-members:
   :show-inheritance:

Constants
^^^^^^^^^

The alphabet, the region and role vocabularies, the coupling offsets and the fitting defaults.
``automodule`` does not reach re-exported data members, so they are named here explicitly.

.. currentmodule:: tcren.potts

.. autodata:: AA
.. autodata:: REGIONS
.. autodata:: ROLES
.. autodata:: CLASSES
.. autodata:: CDR_LOOPS
.. autodata:: GROOVE_REGIONS
.. autodata:: MHC_PARTNER
.. autodata:: OFFSETS
.. autodata:: CROSS_DJ
.. autodata:: DBIN
.. autodata:: DEFAULT_RIDGE
