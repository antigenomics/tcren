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
   * - ``potts_tcr_mhc``
     - 239,093
     - 15,451
     - TCR:MHC groove over the same crystals — twice as many contacts as the peptide interface.
       Reproduce with ``tcren potts fit -s data/Native2026 --partner mhc --balance both -o …``
       (needs the allele reference: ``tcren build-mhc-ref``).

Load either with :meth:`tcren.potts.PottsModel.bundled`.

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
scoring ``F_tcr_mhc`` with Miyazawa–Jernigan and reserving TCRen for TCR:peptide.

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
   \Delta F_i(a) = \Phi^{\mathrm{Potts}}(x_{i \to a}) - \tfrac{1}{20}\sum_b \Phi^{\mathrm{Potts}}(x_{i \to b})

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
