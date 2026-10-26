Screen a cohort
===============

``tcren screen`` ranks separately modelled candidate complexes for one task.
It reports each model's components, a frozen task-matched combination where the
cohort supports one, and the reasons a score is unavailable.

.. code-block:: console

   $ tcren screen -s models.tar.gz -o scores.tsv
   $ tcren screen -s models.tar.gz --confidence confidence.tsv -o scores.tsv

Choose the task and cohort
--------------------------

A cohort holds one side fixed: one receptor with varying peptides is epitope
ranking; one pMHC with varying receptors is TCR ranking. Automatic detection uses
the annotated peptide sequence and CDR3 pair. Keep the MHC background consistent;
this identity check does not establish experimental restriction or binding.
A mixed set, or a set in which neither side varies, is refused with a diagnostic.
``--task epitope`` or ``--task tcr`` declares a task but does not override that check.

.. image:: _static/screening-workflow.svg
   :alt: Two ranking tasks select different structural references, followed by component reporting and a scorable-cohort check.
   :width: 100%

The four exposed components
---------------------------

**Confidence:** ``iptm`` is the generator's interface predicted TM score, supplied
from the input set's metadata or a unique ``complex.id``/``iptm`` table passed with
``--confidence``. tcren does not compute it. Missing confidence is reported.

**Energy:** a directed statistical contact potential sums residue-pair entries on
the observed TCR:peptide contact map:

.. math::

   \Phi_{\mathrm{TP}}(x)=\sum_{(u,v)\in\mathcal C_{\mathrm{TP}}(x)}\phi_{a_u b_v},
   \qquad
   \Delta_s\Phi_{\mathrm{TP}}(x)=\Phi_{\mathrm{TP}}(x)-\Phi_{\mathrm{TP}}(A_s x).

Here :math:`x` is one annotated structure, :math:`\mathcal C_{\mathrm{TP}}(x)` its
residue contacts at a closest-heavy-atom cutoff of 5 Å, :math:`u` a TCR residue,
:math:`v` a peptide residue, and :math:`a_u,b_v` their amino-acid identities.
:math:`\phi_{ab}` is the directed TCRen2 table entry. :math:`A_s` truncates the
selected side :math:`s` to alanine in three dimensions, retaining backbone and
Cβ atoms, and rebuilds the map at the same cutoff. The energy is a dimensionless
statistical contact quantity; lower differences are the declared favourable direction.

.. image:: _static/structural-alanine.svg
   :alt: Truncate side chains in coordinates, rebuild the contact map, and subtract reference energy from real energy.
   :width: 100%

The substituted side follows the side that varies. Epitope ranking uses
``dPhi_tcr_pep`` from peptide substitution; TCR ranking uses
``dFT_tcr_pep_cdr`` from CDR substitution. A fixed-map virtual identity substitution
is a different quantity. ``screen`` always computes the structural reference.

**Geometry:** ``sc_shape`` is the Pearson correlation between finite paired height
samples of the pMHC face and receptor underside in the same frame:

.. math::

   r_h=\frac{\sum_{k\in\mathcal K}(h_k-\bar h)(t_k-\bar t)}
   {\sqrt{\sum_{k\in\mathcal K}(h_k-\bar h)^2}\sqrt{\sum_{k\in\mathcal K}(t_k-\bar t)^2}}.

:math:`\mathcal K` contains the paired interface cells, :math:`h_k,t_k` are their
heights in Å and the bars denote means on those cells. The correlation is
dimensionless; degenerate or insufficient samples are unavailable, not neutral zero.
See :mod:`tcren.topology.surface` for grid and contact-plane definitions.

**Topology:** ``g_recognition_cyclo_frac`` is the recognition contact-graph cycle fraction:

.. math::

   g_{\mathrm{cyclo}}=\frac{E-V+C}{E}, \qquad E>0.

:math:`E` counts unique residue-contact edges, :math:`V` contacted vertices and
:math:`C` connected components. The graph contains CDR:peptide, peptide:MHC and CDR:MHC edges at 5 Angstrom.
CDR1, CDR2 and CDR3 are included; framework, backbone-sequence and internal
peptide edges are excluded. ``g_cyclo_frac`` retains the older CDR:pMHC definition.
This dimensionless graph quantity is carried in the TCR-ranking combination;
topology is excluded from the epitope-ranking combination.

Frozen standardization and combination
--------------------------------------

For each observed component :math:`i`, its frozen mean :math:`\mu_i`, positive
scale :math:`\sigma_i`, declared direction :math:`d_i\in\{-1,+1\}` and fixed weight
:math:`w_i\geq0` give

.. math::

   z_i=d_i\frac{x_i-\mu_i}{\sigma_i}, \qquad
   S=\frac{\sum_{i\in\mathcal O}w_i z_i}{\sqrt{\sum_{i\in\mathcal O}w_i^2}}.

:math:`x_i` is the raw value and :math:`\mathcal O` contains the observed components.
Missing values are omitted from both sums. The combination is unavailable when
no weighted contribution remains. Higher :math:`S` is the declared favourable
direction. It is dimensionless and unbounded, not a binding free energy or a
calibrated probability. Dependence among components is not a claim of independent
evidence or a calibrated Stouffer significance test.

The bundled TCR reference uses three bounded components:

.. math::

   S=(1-0.41q-0.22(1-q))I+0.41qB+0.22(1-q)M.

Here :math:`I` is ipTM, :math:`B` is recognition graph beta1/E and :math:`M`
is empirical whole-mock contact enrichment. The weights are nonnegative and
sum to one. Energy, geometry and frozen z-scores remain diagnostics.

The gate :math:`q` is the inclusive empirical CDF of the unlabelled cohort mean
pLDDT/100 against the five frozen means 0.9179320394646805,
0.92612607434279, 0.9299804672775766, 0.9312745381105264 and
0.9497418534263323. Comparisons round normalized confidence to twelve decimals
only to stabilize numerical ties. All supplied feature rows contribute to the
mean, including rows without binding labels. Supply ``plddt`` on the 0–100 scale
in metadata or ``--confidence``. Missing components and confidence use the frozen
training means and are flagged. The default cohort floor still applies.

Local residue probabilities use conditional mock upper tails for added contact
cycles (atom contacts minus distinct partner residues). Fisher evidence sums
``-2*ln(p)`` over all six CDR loops and peptide residues. Whole-control empirical
survival, with inclusive ties and a plus-one correction, converts it to
``M = 1 - p_mock``. All 2,056 eligible mock graphs have unit weight; each control's
own model was excluded from its local calibration ranks. The packaged histogram
and gzipped TSV calibration table are frozen estimated reference dependencies.
See :doc:`motif-confidence` and :doc:`motif-mixture` for the definitions.

Nothing is fitted when scoring. Cohort composition does affect the measured mean
and therefore its weights; an isolated model is a different context. The two
coefficients were selected using TCRvdb labels. The accepted development macro AP
is 0.768183 over five epitopes and 1,782 classified models; template-absent VDJdb
macro AP is 0.561450 across 16 cohorts and 738 models, versus 0.572973 for the
preceding prepared baseline. These are recorded development/transfer results,
not an untouched test. All per-cohort gains and losses remain in the development
record. Scores are not calibrated binding probabilities.

Epitope ranking retains Stouffer confidence and energy with equal weights.
Custom references without ``score_method`` retain Stouffer combination. Historical
TCR resources ``screen_reference_tcr_v43.json`` and
``screen_reference_tcr_v44.json`` preserve their original formulas.

Template availability and effective weight
------------------------------------------

``--template-status present|absent`` or ``--template-metadata inventory.tsv``
records explicit context. Under the new default, this does not change the shared
confidence rule. Status and inventory cannot be combined. ``--weight`` overrides
only the confidence coefficient while retaining the structural Betti:motif ratio.
Outputs record ``confidence_weight``, ``topology_weight``, ``motif_weight`` and
``gate_q``. Availability does not establish predicted-pose accuracy.

The historical two-weight reference retains 0.55 present/unspecified and 0.90
absent confidence weights; older custom references retain their original gate.
Epitope ranking does not accept template status. Zero-contact discard, cohort
floor and reporting-only flags retain their existing policies.

Energy, orientation and target relief
-------------------------------------

Keep the structural alanine energy alongside ipTM and the aggregate: it measures
CDR:peptide chemistry for receptor ranking and peptide chemistry for epitope
ranking. A gate that omits energy from one aggregate does not invalidate its
separate readout or establish that confidence is always stronger.

The geometry component ``sc_shape`` measures receptor--pMHC height-field
complementarity. TCR:MHC orientation (crossing angle, docking torsion and height)
and receptor-free pMHC topography (surface relief, accessible peptide and surface
chemistry) remain distinct toolkit outputs; they are not additional terms in
this aggregate. ``pitch`` is excluded as generator-confidence leakage.
The graph topology and surface geometry can both describe packing and interface
size, so their presence does not establish independent evidence. Near-orthogonal
confidence and structural directions motivate a standardized sum, but do not
make its denominator a calibrated significance test. Inspect the individual
channels, their missingness, the target surface and docking geometry with flags
before interpreting a ranking.

Discard and cohort floor
------------------------

Models with zero TCR:peptide contacts remain in output with ``no_interface`` true
and a reason in ``score_unavailable``, but no score. They do not count toward
``--min-cohort`` (default 30). Below the floor every row retains its component
readouts without a combined score. A flag response cannot restore discarded or
withheld scores or erase their original unavailability reasons.

Generator confidence reports pose quality; a binder-labelled pair may still have
an off-groove prediction. Reversing such a score would reward modelling failure.
Inspect model quality and provenance before interpreting discordant labels.

Flags and explicit responses
----------------------------

Flags report by default. ``small_cohort`` describes the scorable floor;
``no_interface`` reports zero-contact models; ``low_confidence`` compares cohort
mean ipTM with the frozen reference distribution; ``confidence_disagrees`` reports
opposing confidence and structure rankings. ``forced_pose`` is unavailable unless
strain is standardized against the population used for its reference cut.

``cohort_confidence_z`` is the standardized cohort mean, using the frozen
per-epitope mean 0.74562 and scale 0.035717; ``low_confidence`` fires below -0.6.
These stored constants describe the reference, not query-derived parameters.
A flag is an applicability observation, not a binding label or a validated correction.

``--on-flag`` controls the response to low confidence:

* ``report`` (default): retain the score and report the flags.
* ``structure``: omit confidence from the combination.
* ``geometry``: report the standardized geometry contribution as the score.
* ``withhold``: withhold the combined score.
* ``reverse``: negate the score; opt-in and not recommended as a modelling correction.

Responses change only the combined score. Raw components and missingness remain
visible; ``flag_response`` and ``flag_response_reason`` record an applied response.
References lacking a geometry component cannot use the geometry response.

Output and interpretation
-------------------------

The TSV keeps ``complex.id``, raw components, ``z_*`` contributions, ``missing_*``
flags, ``n_components``, ``screen_score`` and ``score_unavailable``. Cohort flags,
confidence statistics and zero-contact counts accompany each row. The provenance
sidecar records the task, reference hash, confidence weight, floor and flag response.

For receptor ranking, keep ipTM primary and read the structural components beside
it. The combination is complementary evidence, not an established universal
improvement. Report template coverage separately. Peptide-intrinsic surface
properties and a high CPL ranking alone do not establish receptor specificity.
Scientific evaluation, figures and historical results belong to the manuscript
and benchmark records; this page specifies the software's current behavior.

Library API
-----------

.. automodule:: tcren.screen
   :members:

Exploratory contact-profile evaluation
--------------------------------------

``tcren contact-profile-evaluate --contacts contacts.parquet --atoms atoms.parquet
--baseline baseline.parquet -o evaluation -t 0`` evaluates residue involvement
from fresh ``export --quantity contact-atoms --interface tcr_peptide`` and
``export --quantity atom-inventory --interface all`` exports at 5 Angstrom.
The baseline requires unique ``id``, ``epitope``, binary ``y``, ``iptm``,
``prepared`` and ``matched`` columns. This experiment is separate from the
released/default screening score.

Flat atom-pair, residue-pair and zero-inclusive residue tables retain observed
TCR/CDR3 and peptide sequences. Counts are distinct observed heavy atoms,
including backbone; occupation divides contacted by observed available atoms.
It does not assert canonical completeness or polar chemical satisfaction.
Profiles use anchor-excluded CDR3 positions padded after residue 3 and peptide
positions, separately by chain. Missing positions have explicit indicators.
Separate alpha-only and beta-only CDR3 arms accompany the joint profiles.
The state-aware arm retains contacted atom count, any-contact and full-contact
indicators, and occupation restricted to partially contacted residues (zero
otherwise). This separates masses at zero and one without fitting a density.

Five outer folds and three inner folds keep transitive components sharing either
exact CDR3 chain together. Inner folds select ridge regularization, profile family
and nonnegative mixture weights using mean AP on the two priority epitopes named
in the experimental implementation. Imputation and scaling use training rows
only. Outputs retain all five epitope metrics, the no-contact mixture, fitted
parameters, folds and additive logit contributions. The convex mixture is bounded
in [0,1], not an experimentally calibrated binding probability. Descriptive
Mann--Whitney P-values are uncorrected and model-level; related receptors may be
dependent. These outputs are internal exploratory validation, not independent
external validation or automatic promotion of a fitted score.

Contact engagement additions
----------------------------

``tcren contact-engagement-evaluate`` compares engagement additions against a
supplied, fixed ``prepared`` ipTM/Betti baseline. It takes fresh ``--contacts``,
``--atoms`` and ``--baseline`` tables, plus corresponding ``--mock-contacts``,
``--mock-atoms`` and ``--mock-baseline`` inputs. Baselines additionally identify
``dataset`` and ``stratum``; ambiguous outcomes remain in the raw benchmark
inventory for reference exclusion but do not enter classification.

The mock reference excludes benchmark epitopes, either exact CDR3 chain and
coordinate duplicates (the caller checks coordinate hashes before export).
Reference residues receive equal-epitope/equal-CDR3 model weights. Inclusive
amino-acid-specific empirical tails retain zero contacts and use one tail
pseudocount. Zero-contact information is exactly zero. Contributions never
divide by residue atom capacity. Raw counts and endpoint profiles remain
separate candidates; signed categorical likelihood ratios preserve zero,
partial, full and missing-position states. These are composite ranking scores,
without independent-residue or calibrated-probability claims.

Five outer/three inner exact-chain-grouped folds choose regularization and
``S=(1-lambda)*prepared+lambda*contact``. ``--policy priority`` shares the
TTD/NLV-priority choice across epitopes. ``safeguard`` chooses by each epitope's
inner AP, retaining the fixed baseline as a candidate. ``ranking`` additionally
fits a smoothed AP objective anchored to the baseline log odds. Its aggregate
families include an explicit ipTM-minus-baseline contrast, allowing the
confidence/topology balance to change without treating the inputs as
independent evidence. Imputation, indicators and scaling use training rows.
``--policy stable`` retains those AP candidates but admits an addition for a
cohort only when its AP does not fall below the fixed baseline in any of the
three inner folds. This is a selection-stability rule, not a significance test
or a guarantee about outer-fold performance.
``--policy transport`` uses empirical quantiles from each fitted model's training
scores to map its contact output monotonically into the training baseline's
score units before mixture selection. Ties use midrank percentiles. This mapping
uses neither held-out labels nor a query-cohort distribution; its saved knots
are part of the fitted model. It changes score scale, not the physical
per-residue contact definition or residue atom capacity.
``--policy shared`` also considers a compact conditional model trained across
the TCRvdb training cohorts: absolute contacted-atom, engaged-residue and
full-contact totals for each chain/region, four separate type-tail sums and a
confidence contrast. It uses the same grouped splits and training-only mapping.
VDJdb never enters its training or selection.
``--policy bagged`` averages the three fitted member probabilities from the
existing inner training partitions instead of refitting one outer-training
model. Members retain their own imputation, scaling, score mapping and feature
contributions. All members exclude the outer test components; no new split or
bootstrap control is introduced.

Outputs retain all five TCRvdb comparisons, complete precision--recall curves,
fold parameters, feature logit contributions and reference-only VDJdb results
within each template stratum. Exact VDJdb chains are excluded from transfer
training. A positive AP delta does not establish precision dominance at every
recall. Repeated panel-guided development remains exploratory; no default
screening score is replaced by this command.

Exploratory additive contact search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``python -m tcren.contacts.additive --contacts contacts.parquet --atoms atoms.parquet
--topology topology.tsv --baseline baseline.parquet -o evaluation -t 0``
compares maximum TRA/TRB atom-pair counts, maximum distinct contacted TCR atoms,
distinct contacted peptide atoms, and both maxima divided by the peptide atom
count. Inputs are fresh ``export`` outputs at 5 Angstrom (``contact-atoms``,
``atom-inventory``, ``recognition-topology``); receptor contacts cover CDR1/2/3.
Zero-contact counts and ratios are zero. Betti means the existing recognition
``beta1/E`` on CDR:peptide + peptide:MHC + CDR:MHC, without internal edges.

The explicit exploratory score is ``iptm + a*betti + b*contact``. Shared and
separate nonnegative coefficients use zero and 1, 2, 5 times powers of ten
from 1e-5 to 50. Equal-epitope mean average precision selects coefficients;
separation p-values never select a model. Five folds keep receptors sharing
either exact chain together. Held-out scores and both baselines receive the
same training-only, within-epitope midrank transform before pooled evaluation.
Full-fit and held-out outputs remain separate. Two-sided Mann--Whitney tests
report equal-distribution separation with signed rank-biserial effects; these
are exploratory, not external validation or evidence of affinity calibration.

In-sample contact GLMs
~~~~~~~~~~~~~~~~~~~~~~

``python -m tcren.contacts.glm --contacts contacts.parquet --topology topology.tsv
--baseline baseline.parquet -o glm`` fits two pooled, unpenalized binomial GLMs
on all called models with finite ipTM and betti. Here ``betti`` is the previous
prepared score, not the isolated graph cycle fraction. The five predictors are
ipTM, betti and distinct contacted peptide, TRA and TRB heavy-atom counts.
The second fit divides each atom count by the total number of unique
TCR--peptide residue-contact edges; zero-contact densities are zero.
Standardization is numerical only; exported coefficients use input units.
These are deliberately in-sample fits, without held-out accuracy claims.
