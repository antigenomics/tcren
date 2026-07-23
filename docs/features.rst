Interface feature reference
===========================

``tcren recognize`` turns a set of TCR–pMHC structures into one flat table — **one row per
structure**, one interface descriptor per column — for downstream classification and ranking. It is
the single source of truth for every feature used in the TCRen2 benchmarks.

.. code-block:: console

   # 35 core descriptors + P(real):
   tcren recognize -s structures/ -o table.tsv
   # + 18 CDR3-frame + 12 matrix-swap descriptors (65 features):
   tcren recognize -s structures/ -o table.tsv --full
   # + the frozen "good-results" scores p_bind and p_forced:
   tcren recognize -s structures/ -o table.tsv --scores

Output is **TSV**. The first column ``complex.id`` is the structure-file stem (the SHA-256 ``TCR_hash``
for the modelled sets), which is the join key to labels and AlphaFold confidences. ``--features-only``
skips the models; ``--scores`` implies ``--full``. Degenerate or undefined terms are ``NaN``. Every
feature is also available programmatically from :func:`tcren.recognition.recognition_features` (pass
``full=True``); the column-name tuples are :data:`tcren.recognition.RECOGNITION_FEATURES`,
:data:`~tcren.recognition.CDR3_FRAME_FEATURES`, :data:`~tcren.recognition.MATRIX_SWAP_FEATURES` and
:data:`~tcren.recognition.FULL_FEATURES`.

Conventions used below: **tp** = TCR↔peptide interface, **tm** = TCR↔MHC, **pm** = peptide↔MHC.
``F_*`` is a raw interface energy Φ; ``dF_*`` its poly-alanine reference delta ΔΦ; ``e_*`` a per-loop
or per-interface contact energy. Energies are dimensionless statistical-potential units (more negative
= more favorable); TCR↔peptide uses the **TCRen** potential, TCR↔MHC and peptide↔MHC the
**Miyazawa–Jernigan (MJ)** potential.

Core recognition descriptors (35)
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
     - :func:`tcren.pipeline` energy
   * - ``F_tcr_mhc`` / ``e_tcr_mhc``
     - MJ
     - Raw TCR↔MHC interface energy (identical columns, kept for schema stability).
     - :func:`tcren.pipeline` energy
   * - ``F_pep_mhc``
     - MJ
     - Raw peptide↔MHC interface energy.
     - :func:`tcren.pipeline` energy
   * - ``e_cdr12``
     - TCRen
     - TCR↔peptide energy over the CDR1+CDR2 loops only.
     - :func:`tcren.pipeline` energy
   * - ``e_cdr3a``
     - TCRen
     - TCR↔peptide energy over the CDR3α loop only.
     - :func:`tcren.pipeline` energy
   * - ``e_cdr3b``
     - TCRen
     - TCR↔peptide energy over the CDR3β loop only.
     - :func:`tcren.pipeline` energy
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
   * - ``ct_tp_hydrogen_bond``, ``ct_tm_hydrogen_bond``
     - count
     - Hydrogen-bond contacts on each interface.
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
     - Duplicate of ``ct_tp_hydrogen_bond`` (kept for schema stability).

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

Matrix-swap descriptors (12) — ``--full``
-----------------------------------------

The same TCR:peptide contacts scored under **TCRen vs the generic MJ** potential, per interface group;
the per-group difference isolates the recognition-specific component (generic packing cancels because
both potentials read the identical contacts). Groups ``g`` ∈ {``tp`` (whole TCR:peptide), ``cdr12``,
``cdr3a``, ``cdr3b``}. From :data:`tcren.recognition.MATRIX_SWAP_FEATURES`.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Column
     - Unit
     - Description
   * - ``tcren_{g}``
     - TCRen
     - Group energy under the TCRen potential. (``tcren_tp``/``tcren_cdr12``/``tcren_cdr3a``/
       ``tcren_cdr3b`` duplicate the core ``F_tcr_pep``/``e_cdr12``/``e_cdr3a``/``e_cdr3b`` by
       construction — kept for parity.)
   * - ``mj_{g}``
     - MJ
     - Group energy under the generic MJ potential.
   * - ``d_{g}``
     - Δ
     - ``tcren_{g} − mj_{g}`` — the recognition-specific contrast.

Scores (``--scores``)
---------------------

The frozen per-structure "good-results" scores. ``p_real`` / ``p_real_bn`` come from ``recognize`` by
default; ``p_bind`` and ``p_forced`` are added by ``--scores``. All are probabilities in [0, 1].

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
   * - ``p_bind``
     - Binder vs non-binder (screen many TCRs against one epitope).
     - Frozen 5-feature logistic (:func:`tcren.binder.binder_score`); TCRvdb denoised AUC 0.928 vs
       AlphaFold ipTM 0.867.
   * - ``p_forced``
     - Crystal-natural vs AF-forced pose ("too good to be true" hallucination).
     - Frozen 6-feature strain logistic (:func:`tcren.recognition.forced_pose_score`,
       :data:`~tcren.recognition.FORCED_POSE_MODEL`); crystal-vs-AF 5-fold AUC 0.762.

.. note::

   **Cohort-relative combinations live in** :mod:`tcren.cohort`. Scores that z-standardize a
   feature over the *set being scored* — the no-fit :func:`~tcren.cohort.phi_bind`, the
   interface-quality :func:`~tcren.cohort.q_score`, and the crystal-calibrated
   :func:`~tcren.cohort.strain_z` gradient — are not per-structure frozen models, but they are
   computed by ``tcren``, not downstream. Pass the whole cohort you are ranking; pass the crystal
   cohort as ``reference=`` to calibrate strain against natural geometry.

   (Before v2.2.3 these were declared out of scope and lived in manuscript scratch, which left the
   published headline numbers un-regenerable from this package. That is fixed.)
