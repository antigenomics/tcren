Contact-map Potts model
=======================

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
   class, joint models only.

Everything else is asserted uncoupled: the model's statement that linkage falls to zero across
loops and beyond two residues.

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
sites as latent and runs Monte-Carlo EM around the existing fit.

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

   **The latent band has to reach both ways.** ``available_pairs`` measures ``d_heavy`` only out to
   its own ``cutoff``, so every pair that missed carries a null distance. A band built on that table
   is one-sided: on the 17 bundled test structures it selects **44 sites, every one of them a
   contact**, and an EM run on it could remove contacts and never add them. Build the table with
   ``available_pairs(..., probe=cutoff + margin)`` — which populates ``d_heavy`` past the cutoff
   without moving ``sigma`` — and the same band selects **115 sites, 44 contacts and 71 near
   misses**. :func:`tcren.potts.latent_mask` raises rather than proceeding one-sided.

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

Placing a pair without a template
---------------------------------

:func:`tcren.potts.available_pairs` reads :math:`d_{\mathrm{C}\alpha}` off a structure, so every
read-out above is conditional on holding a template of the right shape — and lengths are not evenly
covered. :mod:`tcren.potts.geom` models the distance matrix instead,

.. math::  P\big(d_{\mathrm{C}\alpha} \mid \text{loop},\, i,\, L,\ j,\, M,\ \text{MHC class}\big)

and emits the same schema, so :func:`tcren.potts.site_codes`, :func:`tcren.potts.eta` and every
scoring function accept it unchanged.

**Length generalisation is the normalised coordinate, not the parameterisation.** A loop of length
:math:`L` against a peptide of length :math:`M` is resampled onto :math:`(u, v) = (i/(L+1),
j/(M+1))` — the graphon coordinates :func:`tcren.contactmap.registered_map` already uses — so
:math:`L = 10` and :math:`L = 17` land on one surface that reads back at any pair of lengths,
including a pair the fit never saw.

**Three loops, not one.** On the bundled crystals CDR3 carries **69.3% of TCR:peptide contacts**
and CDR1 a further **26.5%**; CDR1 + CDR2 + CDR3 is **98.4%**. A CDR3-only prior would cap its own
recall below 0.70 by construction, so all three carry a surface. CDR1 and CDR2 are germline and
barely vary in length — they need the surface, not the length generalisation.

**Feasibility is checked against the loop physics, not against the fit.**
:func:`tcren.geometry.reach_max` gives, in closed form from the discrete-Frenet virtual-bond chain,
the maximum apex reach of a length-``L`` Ω-loop pinned at neck separation :math:`R`. A predicted
distance beyond it is one no loop of that length can adopt, and the emitted table flags it in
``feasible`` rather than scoring it silently.

``sigma`` on a predicted table is a **placeholder**, not data: nothing here observed a contact, so
it is zero throughout and ``d_heavy`` is null. ``d_ca_sd`` carries the surface's own spread at that
cell, which is what says how much to trust the placement.

**The fitted prior ships.** ``geom_tcr_peptide`` is averaged over the same 359 alpha-beta
``Native2026`` crystals :func:`tcren.potential.tcren2` is derived from, twelve surfaces at grid 12
(TRA/TRB x CDR1/CDR2/CDR3 x MHCI/MHCII), and :meth:`tcren.potts.GeomPrior.bundled` loads it. Fit
your own with

.. code-block:: console

   tcren potts geom -s <structures> -o geom.json --grid 12 --exclude 3qib,3qiu

``--exclude`` is how a held-out prediction is built; the design belongs to whoever is evaluating,
so the command only offers the exclusion. The MHC groove is annotated whether or not it is asked
for, because the surface is keyed on MHC class and without it every crystal reports class I.

Scoring a map with no structure of the complex is then one command:

.. code-block:: console

   tcren potts predict -s model.pdb.gz -m potts_tcr_peptide_sc -o predicted.tsv
   tcren potts predict --sequences loops.tsv -m potts_tcr_peptide_sc -o predicted.tsv

The first opens the structure for its loop sequences, its peptide and its MHC class **only**, and
joins the crystal's own ``sigma`` and ``d_ca`` back on so the table can be scored against what the
crystal did. The second reads no structure at all: a TSV naming ``id``, ``peptide``, optionally
``mhc_class``, and any of the six loop columns ``TRA:CDR1`` ... ``TRB:CDR3``. Both emit the
``potts contacts`` columns plus ``d_ca_sd`` and ``feasible``.

The difference from :func:`tcren.potts.contact_probabilities` on ``available_pairs`` is not the
model — it is the same Hamiltonian and the same block Gibbs — but **what the 15 A cut is applied
to**: there, the Ca distance the complex was solved at; here, the distance the surface predicts.

Average over geometries, not over surfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Everything above scores **one** map, read off the mean surface. That is the wrong order of
operations, and it costs most of the map. Two reasons, and neither is a matter of degree:

* A contacting pair sits at the **near end** of the distribution being averaged. Over every pair
  the measured and predicted geometries share, the predicted :math:`d_{\mathrm{C}\alpha}` runs
  0.2–2.3 A long; over the pairs that actually contact it runs **2.6–4.5 A** long. The model's
  distance term swings from :math:`+5.77` at 4.5–5.0 A to :math:`-2.33` at 11.0–11.5 A, so those
  angstroms are worth about four natural-log units of one-body field.
* **The distances in one map are not independent.** Consecutive C\ :math:`\alpha` atoms are a
  fixed virtual bond apart and a loop cannot bend arbitrarily, so no set of per-pair marginals
  recombines into a geometry a chain could adopt. Averaging them lets every pair be close at once.

:func:`tcren.potts.predict_contact_map` averages the other way round. Every crystal in
:attr:`tcren.potts.GeomPrior.surfaces` is read back at the target lengths, the map is scored under
*that* geometry, and the contact probabilities are averaged:

.. math::  \bar p_a = \frac{1}{S}\sum_{s} p_a\big(d^{(s)}\big)\,
           \mathbb{1}\big[d^{(s)}_a \le r\big]

Every member is a chain that was actually solved, so the virtual-bond constraint holds by
construction. Held out by structure over **364** ``Native2026`` crystals:

.. list-table::
   :header-rows: 1

   * - estimator
     - mean :math:`\sum_a p_a`
     - ratio to observed
     - median contact ROC-AUC
   * - mean surface
     - 9.80
     - 0.466
     - 0.9242
   * - **geometry ensemble**
     - **24.62**
     - **1.171**
     - **0.9325**

against **21.02** contacts observed per structure on CDR1–CDR3 (median 21, s.d. 6.08). Both the
level and the ordering improve; the mean surface loses more than half the contacts.

.. code-block:: console

   tcren potts geom -s <structures> --grid 32 --keep-surfaces -o prior.npz
   tcren potts predict --sequences loops.tsv --ensemble -m potts_tcr_peptide_sc -o predicted.tsv

``--keep-surfaces`` needs an ``.npz`` output, because a stack of surfaces does not fit in JSON;
distances are stored to 0.01 A as ``uint16``, four decimal places finer than the 0.5 A bin the
model reads. The shipped ``geom_tcr_peptide_ens`` carries **2,154** surfaces from the same 359
crystals at grid 32, 3.5 MB, and ``--ensemble`` selects it by default. The scored table gains
``p_avail``, the fraction of reference geometries placing that pair within the radius at all, and
costs about 0.16 s per member.

What the ensemble does **not** buy is receptor specificity: the correlation between the predicted
and the observed contact count across the 364 structures is Pearson :math:`+0.23`
(Spearman :math:`+0.18`), and the per-structure mean absolute error, 5.95 contacts, is about the
spread of the population itself. It predicts the typical interface at the right level, which is
what loop lengths alone can carry.

Scoring the agreement: placement and level
------------------------------------------

:func:`tcren.potts.map_agreement` reads any scored site table --- from
:func:`contact_probabilities`, :func:`predict_contact_map`, ``tcren potts contacts`` or
``tcren potts predict`` --- and reports **two** numbers per map, because no single one answers
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
equal proportion, which is a deterministic map; a calibrated probabilistic model cannot reach it
against a single realisation. On the four receptors of the manuscript's Figure 2c the crystal's own
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

   # the free-energy effect of every substitution at every peptide position
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

``map`` closes those pairs onto a coarser grid, which is where they become comparable with an
experiment. The residues of one CDR loop are distinct pairs with different marginals, so the number
of simultaneous contacts is Poisson-binomially distributed and has no closed form; the event "at
least one" does,

.. math:: P(N \ge 1) \;=\; 1 - \prod_j \left(1 - p_j\right),

with :math:`p_j` the model marginal of pair :math:`j` in the group. That is the quantity a molecular
-dynamics trajectory reports as the fraction of frames in which a loop touches a peptide position,
so ``--by loop`` is directly comparable with a measured contact-frequency map. ``--by position``
collapses the loops and reads how engaged each peptide residue is expected to be, before any residue
identity is scored. Emitted columns are ``p_any`` (the frequency above), ``p_expected``
(:math:`\sum_j p_j`, the expected number of contacts), ``n_pairs``, ``n_observed`` and ``observed``.
The sum is accumulated in :math:`\log(1 - p)`, so a twelve-residue loop does not underflow and a
saturated pair returns 1 exactly rather than ``nan``.

These are contact **frequencies** — dimensionless, in :math:`[0, 1]`. They are not free energies and
carry no :math:`k_\mathrm{B}T`, so they belong to the diagnostic and importance side of the model
rather than to any energy block; ``score``'s ``neg_energy`` and :func:`~tcren.potts.peptide_free_energy`'s ``log Z0`` are the quantities with units.

``score``, ``contacts`` and ``map`` all take ``--workers`` (default: every core). The per-structure
numbers are functions of ``(seed, pdb.id)`` alone, so splitting the work changes nothing.

Substituting a residue: the free energy, not the frequency
----------------------------------------------------------

``map`` reads how engaged a position is expected to be *before any residue identity is scored*.
:func:`~tcren.potts.peptide_free_energy` reads what happens when the identity changes. The partner
residue enters :math:`\eta` in exactly two places, so for a site :math:`s` and a candidate
residue :math:`a`

.. math:: \eta_s(a) = r_s + h^{\mathrm{par}}(a) + J(\mathrm{rec}_s, a)

with :math:`r_s` — the intercept, the receptor field, the distance profile, region, role and class
— independent of :math:`a`. Threading :math:`a` through partner position :math:`i` moves every site
carrying that position, and the interface free energy moves with it:

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
score in the package makes. A partner position carrying two different residues has no sequence to
substitute into and is rejected rather than averaged.

Unlike ``map``'s frequencies this **is** an energy: :math:`\log Z_0` carries :math:`k_\mathrm{B}T`.

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

**The observation enters as one additive field.** Conditioning on the crystal through
:math:`P(\sigma_a{=}1 \mid d_a) = \operatorname{sigmoid}((c - d_a)/m)` puts a per-site constant into
:math:`\eta`, so the E-step is the existing block-Gibbs sampler called with a modified vector — no
clamping machinery, no second code path. One knob spans the family: ``margin`` → 0 pins the crystal
map, ``margin`` = 0.5 Å releases the threshold band alone, and a large ``margin`` lets the crystal
enter only through the fit.

**The pair term is pinned, not free.** ``J`` carries one coefficient — ``beta_matrix``, the scale on
the current matrix — rather than 400 cells. That is what closes the loop: a free ``J`` would absorb
the pair structure and leave the re-derived matrix no route back into :math:`P(\sigma)`. It also
makes ``beta_matrix`` reportable as the inverse temperature relating TCRen units to contact
log-odds. :func:`tcren.potts.fit_potts` accepts a :class:`~tcren.potential.Potential` object for
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
weighted mass of 309.000000 against 309 observed contacts. That is the check that the two halves of
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
