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
A mixed set is refused with a diagnostic. Automatic mode also refuses a set in
which neither side varies. For a single model, explicitly use ``--task tcr`` with
``--confidence-context nearest-plddt --min-cohort 1``. The
``nearest-mean-plddt`` option uses the observed mean of a declared small same-pMHC
set; see :doc:`motif-confidence` for support and missing-value rules.

.. image:: _static/screening-workflow.svg
   :alt: Two ranking tasks select different structural references, followed by component reporting and a scorable-cohort check.
   :width: 100%

Components and diagnostics
--------------------------

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

Default TCR combination
-----------------------

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
coefficients were selected during labelled development; frozen application is
not an independent validation. Scores are not calibrated binding probabilities.

Epitope ranking retains Stouffer confidence and energy with equal weights.
Custom references without ``score_method`` retain Stouffer combination. Historical
TCR resources ``screen_reference_tcr_v43.json`` and
``screen_reference_tcr_v44.json`` preserve their original formulas.

Epitope and legacy standardization
----------------------------------

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
Epitope ranking does not accept template status. Zero-contact rejection, cohort
floor and reporting-only flags apply to both tasks.

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
size. Their association does not establish independent evidence; the default
TCR score is the bounded mixture defined above, not a significance test. Inspect the individual
channels, their missingness, the target surface and docking geometry with flags
before interpreting a ranking.

Zero-contact rejection and cohort floor
---------------------------------------

Models with zero TCR:peptide residue contacts at the 5 Å heavy-atom cutoff receive
``screen_score=0``: no peptide contact means no binding interface in the supplied
model. They remain in output with ``no_interface`` true and the structural
rejection reason in ``score_unavailable``. The CLI emits a warning with the exact
number of zero-contact models. A high ipTM cannot override this rejection;
downstream exports must retain zero rather than substitute raw ipTM.

Zero-contact models do not count toward ``--min-cohort`` (default 30). Below the
floor their scores remain zero; other rows retain components without a combined
score. Explicit flag responses preserve these zeros and their original reasons.
Zero is the lower bound of the default TCR and raw ipTM/Betti scores. With a
historical standardized reference it is a rejection marker, not necessarily the
lowest numerical score. This behavior applies from version 4.4.1 onward; older
versions withheld the zero-contact score.

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

For receptor ranking, sort by ``screen_score`` with higher values first and retain
ipTM and the structural components for interpretation. Performance depends on
the evaluation population; the combination is not a binding decision. Report template coverage separately. Peptide-intrinsic surface
properties and a high CPL ranking alone do not establish receptor specificity.
Scientific evaluation, figures and historical results belong to the manuscript
and benchmark records; this page specifies the software's current behavior.

Library API
-----------

.. automodule:: tcren.screen
   :members:
