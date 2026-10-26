Working with modelled complexes
===============================

A generated structure describes a proposed pose. Generator confidence, observed
contacts and a statistical contact potential describe different aspects of that
pose. Start with the task and retain those observations alongside the combined
ranking rather than treating one number as a binding decision.

Prepare one cohort
------------------

For peptide ranking, supply candidate complexes for a fixed receptor. For receptor
ranking, supply candidate receptors against one pMHC background. Keep MHC identity
consistent within the cohort. Automatic cohort detection reads the annotated
peptide and CDR3 pair; it is not an experimental identity or restriction check.
A CDR3 excludes the conserved junction anchors; use the supplied chain annotation.

.. code-block:: console

   $ tcren screen -s models.tar.gz --confidence confidence.tsv -o screen.tsv

The optional confidence table has unique ``complex.id``, ``iptm`` and ``plddt`` columns (pLDDT on 0–100).
Without it, ``screen`` looks for the set's ``metadata.tsv``. If ipTM is unavailable,
the default TCR output flags missing values and uses the frozen imputation policy.
See :doc:`screening` for equations, the structural reference and missingness rules.

Read the output before using the ranking
----------------------------------------

#. Check ``no_interface`` and ``score_unavailable``. A zero-contact model is retained
   but discarded from scoring; a cohort below the scorable floor reports components.
#. Inspect the raw component columns, ``z_*`` contributions and ``missing_*`` flags.
   ``n_components`` reports observed inputs rather than an imputed complete vector.
#. Read ``cohort_mean_iptm``, ``cohort_confidence_z`` and the cohort flags. Flags report
   by default; they do not automatically change the score.
#. For receptor ranking, sort finite ``screen_score`` values from highest to lowest.
   Retain ipTM as a comparator and inspect ``gate_q`` and the component weights.
#. Use the provenance sidecar to recover the task, frozen-reference hash, weight,
   cohort floor and requested flag response.

A label describes the candidate pair, not whether the generator successfully
modelled it. Inspect off-groove poses before interpreting a low confidence score
as a failure of ranking. ``--on-flag reverse`` is opt-in and is not recommended
as a correction for such modelling failures.

Inspect the components
----------------------

.. code-block:: console

   $ tcren features -s models.tar.gz -i interface,topology,energetics -o features.tsv
   $ tcren contacts -s models.tar.gz --interface tcr_peptide -o contacts.csv
   $ tcren scoring -s models.tar.gz --delta --reference-mode structural -o energies.csv
   $ tcren export -s models.tar.gz --quantity receptor-reference --reference-mode structural -o receptor.tsv

``features`` provides individual descriptors; its existing ``dPhi`` columns use
the historical fixed-map route. For the structural epitope alanine difference use
``scoring --delta --reference-mode structural``. For contacted-CDR substitution
use the receptor-reference export. ``screen`` computes the task-appropriate
structural reference itself. :doc:`screening` defines the difference explicitly.

``surface`` inspects pMHC relief, hydropathy and charge in the supplied bound
conformation, and the receptor underside in the same frame. Peptide-intrinsic
surface properties do not by themselves establish receptor specificity.
``inspect-models`` retains confidence and descriptor outputs; ``associate``
compares quantities on identical rows within the declared task and strata.
Report template availability rather than inferring it from a confidence value.

Other supported workflows
-------------------------

Use ``score`` for fixed-template peptide threading, ``ddg`` for explicitly defined
mutation readouts, and ``md``/``mechanics`` for their separate trajectory and
mechanical quantities. None makes a contact log-odds a physical binding free energy.
The historical ``recognize``, ``assess`` and fitted Potts APIs remain available
for reproduction with their estimated dependencies documented in :doc:`fit-free`.
They are not the current cohort screen. The withdrawn integral commands are no
longer supported. Scientific campaign results and publication figures are maintained
with the manuscript and benchmark records, separately from this API guide.
