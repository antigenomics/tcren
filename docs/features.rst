Interface feature reference
===========================

``tcren recognize`` turns a set of TCR–pMHC structures into one flat table — **one row per
structure**, one interface descriptor per column — for downstream classification and ranking. It is
the single source of truth for every feature used in the TCRen2 benchmarks.

.. code-block:: console

   # 34 core descriptors + P(real):
   tcren recognize -s structures/ -o table.tsv
   # + the 18 CDR3-frame descriptors (52 features):
   tcren recognize -s structures/ -o table.tsv --full
   # + the frozen "good-results" scores p_bind and p_forced:
   tcren recognize -s structures/ -o table.tsv --scores

Output is **TSV**. The first column ``complex.id`` is the structure-file stem (the SHA-256 ``TCR_hash``
for the modelled sets), which is the join key to labels and AlphaFold confidences. ``--features-only``
skips the models; ``--scores`` implies ``--full``. Degenerate or undefined terms are ``NaN``. Every
feature is also available programmatically from :func:`tcren.recognition.recognition_features` (pass
``full=True``); the column-name tuples are :data:`tcren.recognition.RECOGNITION_FEATURES`,
:data:`~tcren.recognition.CDR3_FRAME_FEATURES` and :data:`~tcren.recognition.FULL_FEATURES`.
:data:`tcren.recognition.DESCRIPTORS` gives every column's family and whether the receptor enters its
definition; :func:`tcren.recognition.descriptors` selects on both (see :ref:`descriptor-families`).

Conventions used below: **tp** = TCR↔peptide interface, **tm** = TCR↔MHC, **pm** = peptide↔MHC.
``F_*`` is a raw interface energy Φ and ``dF_*`` its poly-alanine reference delta ΔΦ — every energy
column is named ``F``, because the potential is fixed by the interface rather than chosen per column:
TCR↔peptide uses the **TCRen** potential, TCR↔MHC and peptide↔MHC the **Miyazawa–Jernigan (MJ)**
potential. Energies are in dimensionless statistical-potential units (more negative = more
favourable); they are log-odds ratios of contact frequencies and are *not* in kT.

.. _descriptor-families:

Families, and which descriptors involve the receptor
----------------------------------------------------

Every emitted column is catalogued in :data:`tcren.recognition.DESCRIPTORS` as
``name -> (family, involves_tcr)``, and selected with :func:`tcren.recognition.descriptors`.
The three families are the three physical channels the method reports:

``geometry``
    Coordinates, docking angles, and the contact topology and chemistry read off them — the kind of
    quantity Eq. Q is built from.

``physics``
    Statistical-potential interface energies ``F`` and their poly-alanine references ``dF``.

``kinetics``
    The interface as a network of breakable springs: stiffness, anisotropy, strain, rupture, and the
    residues coupling the pre-formed scaffold to the interface (``tcren mechanics``, plus the
    contact-fragility terms ``recognize`` emits).

A fourth group, ``score``, holds the fitted and cohort-relative composites (``p_real``,
``p_real_bn``, ``p_forced``, ``p_bind``, ``q_bind``, ``s_strain``). These are **outputs** built from
the descriptors above and must never be fed back in as inputs, so
:func:`~tcren.recognition.descriptors` omits them unless ``with_scores=True``.

``involves_tcr`` is ``False`` for five columns — ``F_pep_mhc``, ``dF_pep_mhc``, ``mhc_class_bin``,
``F_pep_int`` and ``n_pep_int`` — each computed from the peptide and MHC alone. Two structures of the same epitope
on the same allele share their values whatever the receptor, so such a column carries **cohort
identity** rather than interface physics, and a model given one can reach a cohort-level label
without learning anything about recognition. Any analysis whose question is about receptors should
select ``tcr_only=True``::

    from tcren.recognition import descriptors

    descriptors("physics", tcr_only=True)
    # ('F_tcr_pep', 'F_tcr_mhc', 'F_cdr12', 'F_cdr3a', 'F_cdr3b', 'dF_tcr_pep')

Core recognition descriptors (34)
---------------------------------

The base feature set the shipped real-vs-shuffled recognizers consume
(:data:`tcren.recognition.RECOGNITION_FEATURES`), emitted by ``tcren recognize`` with no extra flags.

Coverage & burial
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 16 44 20

   * - Column
     - Unit
     - Description
     - Source
   * - ``extent``
     - count
     - Distinct TCR interface residues contacting the pMHC (interface size).
     - contact map
   * - ``chain_balance``
     - ratio [0, 0.5]
     - ``min(a, b) / (a + b)`` over TCR:peptide contacts by TCR chain (TRA=a, TRB=b); 0.5 = both
       chains engage equally, 0 = one chain only (a degenerate/forced pose signature).
     - contact map
   * - ``n_contacts_tp``
     - count
     - Number of TCR↔peptide residue–residue contacts.
     - contact map
   * - ``n_pep_contacted``
     - count
     - Distinct peptide residues contacted by the TCR.
     - contact map
   * - ``n_contacts_tm``
     - count
     - Number of TCR↔MHC residue–residue contacts.
     - contact map
   * - ``burial``
     - Å²
     - Interface ΔSASA = SASA(TCR) + SASA(pMHC) − SASA(complex), biopython Shrake–Rupley.
     - biopython SASA

Docking geometry
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 16 44 20

   * - Column
     - Unit
     - Description
     - Source
   * - ``pitch``
     - degrees
     - Incident (tilt) angle of the TCR over the pMHC groove — a clean structural angle.
     - :func:`tcren.orient.docking_angles`
   * - ``crossing``
     - degrees
     - TCR crossing (scanning) angle relative to the groove long axis.
     - :func:`tcren.orient.docking_angles`
   * - ``dock_d``
     - Å
     - MHC-stub → TCR-stub rigid-body separation (native TCRdock geometry).
     - ``orient.tcrdock_geometry``
   * - ``dock_torsion``
     - radians
     - Rigid-body dihedral of the TCR about the MHC stub (circular; wraps at ±π).
     - ``orient.tcrdock_geometry``
   * - ``dock_tcr_uy``, ``dock_tcr_uz``
     - unit
     - y/z components of the TCR stub unit vector in the MHC frame.
     - ``orient.tcrdock_geometry``
   * - ``dock_mhc_uy``, ``dock_mhc_uz``
     - unit
     - y/z components of the MHC stub unit vector.
     - ``orient.tcrdock_geometry``

Interface energies
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 16 44 20

   * - Column
     - Unit
     - Description
     - Source
   * - ``F_tcr_pep``
     - TCRen
     - Raw TCR↔peptide interface energy (whole interface, all TCR regions).
     - :mod:`tcren.pipeline` energy
   * - ``F_tcr_mhc``
     - MJ
     - Raw TCR↔MHC interface energy.
     - :mod:`tcren.pipeline` energy
   * - ``F_pep_mhc``
     - MJ
     - Raw peptide↔MHC interface energy.
     - :mod:`tcren.pipeline` energy
   * - ``F_cdr12``
     - TCRen
     - TCR↔peptide energy over the CDR1+CDR2 loops only.
     - :mod:`tcren.pipeline` energy
   * - ``F_cdr3a``
     - TCRen
     - TCR↔peptide energy over the CDR3α loop only.
     - :mod:`tcren.pipeline` energy
   * - ``F_cdr3b``
     - TCRen
     - TCR↔peptide energy over the CDR3β loop only.
     - :mod:`tcren.pipeline` energy
   * - ``dF_tcr_pep``
     - TCRen
     - Poly-alanine reference delta of the TCR↔peptide energy (geometry-normalized ΔΦ).
     - :func:`tcren.ddg.reference_delta`
   * - ``dF_pep_mhc``
     - MJ
     - Poly-alanine reference delta of the peptide↔MHC energy.
     - :func:`tcren.ddg.reference_delta`

Contact-type tallies
~~~~~~~~~~~~~~~~~~~~~~

Per-interface counts of contacts of each chemical type (``tp`` = TCR↔peptide, ``tm`` = TCR↔MHC),
from :func:`tcren.contact_types.contact_type_counts`.

.. list-table::
   :header-rows: 1
   :widths: 34 16 50

   * - Columns
     - Unit
     - Description
   * - ``ct_tp_salt_bridge``, ``ct_tm_salt_bridge``
     - count
     - Salt-bridge contacts on each interface.
   * - ``ct_tm_hydrogen_bond``
     - count
     - Hydrogen-bond contacts on the TCR↔MHC interface. (The TCR↔peptide count is emitted once,
       as ``n_hbond`` — the name Eq. Q uses.)
   * - ``ct_tp_aromatic``, ``ct_tm_aromatic``
     - count
     - Aromatic (π-stacking) contacts on each interface.
   * - ``ct_tp_hydrophobic``, ``ct_tm_hydrophobic``
     - count
     - Hydrophobic contacts on each interface.
   * - ``ct_tp_other``, ``ct_tm_other``
     - count
     - Remaining contacts on each interface.
   * - ``n_hbond``
     - count
     - Hydrogen-bond contacts on the TCR↔peptide interface; a term of Eq. Q.

MHC class
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 16 64

   * - Column
     - Unit
     - Description
   * - ``mhc_class_bin``
     - 0/1
     - 1 if any MHC chain is MHC class II, else 0 (from MHC annotation).

Interface quality — clashes & contact stability
-----------------------------------------------

Coordinate-only reads of forced-pose quality, always emitted by ``recognize`` (not part of the 35
model features): a steric-clash burden (:mod:`tcren.clashes`) and TCR:peptide contact fragility
(:mod:`tcren.stability`). Both are computed natively (``_geom``).

.. list-table::
   :header-rows: 1
   :widths: 20 16 44 20

   * - Column
     - Unit
     - Description
     - Source
   * - ``n_clashes``
     - count
     - Peptide↔partner heavy-atom pairs overlapping by more than 0.4 Å (Bondi vdW radii); the
       signature of a non-physical / wrong-register pose.
     - clashes
   * - ``clash_score``
     - Å
     - Summed overlap of all clashing pairs — a total steric-burden measure.
     - clashes
   * - ``exp_lost``
     - count
     - Expected TCR:peptide contacts lost under a 1 Å isotropic shift, ``Σ clip((1 − margin)/2, 0, 1)``
       over contacts (``margin = 5 − dmin``).
     - stability
   * - ``mean_margin``
     - Å
     - Mean contact margin ``5 − dmin`` over TCR:peptide contacts; larger = contacts sit deeper below
       the 5 Å cutoff.
     - stability
   * - ``frac_robust``
     - ratio [0, 1]
     - Fraction of TCR:peptide contacts with margin ≥ 1 Å (robust to a 1 Å shift).
     - stability

CDR3-frame descriptors (18) — ``--full``
----------------------------------------

The FramePose layer the whole-TCR features miss: each CDR3 loop projected onto the **pMHC groove
frame** (basis ``u, w, n``; origin = peptide Cα centroid). Computed by
:func:`tcren.recognition.recognition_features` with ``full=True``
(:data:`~tcren.recognition.CDR3_FRAME_FEATURES`). The ``cdr3b_*`` strain terms are the load-bearing
signal for forced-pose / hallucination detection. Columns are prefixed by loop (``cdr3a_`` for TRA,
``cdr3b_`` for TRB):

.. list-table::
   :header-rows: 1
   :widths: 20 16 64

   * - Suffix
     - Unit
     - Description
   * - ``reach``
     - Å
     - Distance from loop Cα centroid to the groove origin (how far the loop reaches out).
   * - ``ou``, ``ow``, ``on``
     - unit
     - Projection of the (centroid − origin) direction onto ``u``, ``w``, ``n`` — where over the
       groove the loop sits.
   * - ``au``, ``aw``, ``an``
     - unit
     - Projection of the loop N→C axis onto ``u``, ``w``, ``n`` — the loop's orientation over the
       groove (FramePose orientation).
   * - ``topep``
     - Å
     - Minimum Cα–Cα distance from the loop to the peptide — engagement depth.
   * - ``ext``
     - Å
     - Loop end-to-end extension ``|Cα_C − Cα_N|``.

The intra-peptide term (``--full``)
-----------------------------------

The three interface energies above all sum over contacts between two **different** chains, so a
peptide held in its bound conformation by its own side chains costs the same as one that is not.
:func:`tcren.intra_peptide_energy` is that omitted term, and ``recognize --full`` emits it
(:data:`~tcren.recognition.PEPTIDE_INTERNAL_FEATURES`). It is computed over
:func:`tcren.peptide_internal_contacts` — heavy atoms within 5 Å, sequence separation ≥ 3 — under a
**symmetrised** potential, since an intra-chain pair has no ``from``/``to`` orientation to respect.

The 5 Å cutoff is the same contact definition the rest of the package uses, so an internal contact
and an interface contact mean the same thing. The separation floor is what does the work: it
excludes pairs that touch because they are bonded. Over the 17 deposited complexes in
``tests/assets/pdb`` the count is 18 contacts at ``|i−j| ≥ 3`` and 134 at ``|i−j| ≥ 2``, and that
sevenfold jump is the ``i``/``i+2`` pairs of an extended chain — geometry rather than folding.

That also sets expectations for the term's size: a canonical extended class-I 9-mer makes **zero to
two** internal contacts, so this separates candidates only where the peptide is genuinely bulged or
packed against itself. It is off everywhere by default.

.. list-table::
   :header-rows: 1
   :widths: 20 16 64

   * - Column
     - Unit
     - Description
   * - ``F_pep_int``
     - MJ
     - The peptide's contact energy with **itself**, symmetrised potential. Lower = more
       favourable, as everywhere in tcren.
   * - ``n_pep_int``
     - count
     - How many such contacts the peptide makes.

As a **scoring term** rather than a descriptor, it is opt-in at each layer, weighted by ``w`` and
added to the energy it accompanies (``w = 0``, the default, computes nothing and leaves every score
byte-identical):

.. code-block:: console

   $ tcren score -s c.pdb -c candidates.txt --intra-weight 0.5   # score = Φ_TP + w·E_intra
   $ tcren scoring -s c.pdb --intra-weight 0.5                   # reports F_pep_int, folds w·it into F_total

.. code-block:: python

   from tcren import ContactMap, intra_peptide_energy, score_peptides
   from tcren.potential import mj, tcren

   cm = ContactMap.from_structure(structure, peptide_internal=True)   # required for the term
   intra_peptide_energy(cm, mj())                                     # the native peptide's own energy
   intra_peptide_energy(cm, mj(), peptide="KQWLVWLFL")                # a candidate threaded onto its pose
   score_peptides(cm, candidates, tcren(), intra_weight=0.5, intra_potential=mj())

The term's potential defaults to MJ, not TCRen: TCRen is derived from TCR↔peptide contacts and says
nothing about the contacts a chain makes with itself.

Scores (``--scores``)
---------------------

``p_real`` / ``p_real_bn`` come from ``recognize`` by default. ``--scores`` adds the **recommended
fit-free** cohort scores ``q_bind`` / ``s_strain`` (see :mod:`tcren.cohort`) plus the fitted
``p_bind`` / ``p_forced`` (retained for reproducibility). The frozen models are probabilities in
[0, 1]; ``q_bind`` / ``s_strain`` are cohort z-scores (unbounded, centred on the input set).

.. list-table::
   :header-rows: 1
   :widths: 16 40 44

   * - Column
     - What it discriminates
     - Model
   * - ``p_real``
     - Genuine TCR–pMHC recognition interface vs a wrong-TCR shuffle.
     - Distribution-aware Bayesian logistic (:class:`tcren.recognition.BayesianLogisticRecognizer`),
       trained on Shuffled2026 decoys.
   * - ``p_real_bn``
     - Same, via a conditional-linear-Gaussian Bayes net.
     - :class:`tcren.recognition.GaussianBNClassifier`.
   * - ``q_bind``
     - Binder vs non-binder (screen many TCRs against one epitope). **Recommended.**
     - Fit-free :func:`tcren.cohort.q_score` (``Q``): equal-weight mean of 5 within-cohort z-scores.
       TCRvdb **raw-label** macro AUC ~0.80 vs AlphaFold ipTM 0.794; generalises across cohorts where
       ``p_bind`` does not (benchmark ledger C25).
   * - ``s_strain``
     - Crystal-natural vs AF-forced pose. **Recommended.**
     - Fit-free :func:`tcren.cohort.strain_z` (``S_strain``): signed z of the strain terms, grading
       crystal < AF-real < AF-decoy. Reproducible; no training set.
   * - ``p_bind``
     - Binder vs non-binder — the fitted counterpart of ``q_bind``.
     - Frozen 5-feature logistic (:func:`tcren.binder.binder_score`); TCRvdb **raw-label** macro AUC
       0.796 vs AlphaFold ipTM 0.794. Matches ``Q`` in-sample but does not transfer; prefer ``q_bind``.
   * - ``p_forced``
     - Crystal-natural vs AF-forced pose — the fitted counterpart of ``s_strain``.
     - Frozen 6-feature strain logistic (:func:`tcren.recognition.forced_pose_score`,
       :data:`~tcren.recognition.FORCED_POSE_MODEL`); crystal-vs-AF 5-fold AUC 0.762. Coefficients
       not re-derivable (ledger C23); prefer ``s_strain``.

.. note::

   **Cohort-relative combinations live in** :mod:`tcren.cohort`. Scores that z-standardize a
   feature over the *set being scored* — the no-fit :func:`~tcren.cohort.phi_bind`, the
   interface-quality :func:`~tcren.cohort.q_score`, and the crystal-calibrated
   :func:`~tcren.cohort.strain_z` gradient — are not per-structure frozen models, but they are
   computed by ``tcren``, not downstream. Pass the whole cohort you are ranking; pass the crystal
   cohort as ``reference=`` to calibrate strain against natural geometry.

   (Before v2.2.3 these were declared out of scope and lived in manuscript scratch, which left the
   published headline numbers un-regenerable from this package. That is fixed.)
