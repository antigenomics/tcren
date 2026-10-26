Residue contact motifs above background
=======================================

This experimental CLI workflow uses all six CDR loops and the epitope. It does
not change released screening scores. A repertoire motif can be compared with
control rearrangements; here a local paratope--epitope contact graph is compared
with mock contact graphs. This is a local contact-subgraph analogy,
without an assumption of formal self-similarity.

Residue-side diagnostic
-----------------------

The diagnostic aggregate is H=sqrt(T*P): T averages the six CDR-loop means and
P averages epitope residues. A residue contributes either mock rarity,
``1-p_mock``, or positive enrichment, ``max(0,log2(p_cognate/p_mock))``.
These are ordinary means and log frequency ratios, with no arbitrary powers.
For the final separately weighted ipTM, whole-complex Betti and Fisher-enrichment
score, use :doc:`motif-mixture`. Earlier ``residue-surprise-evaluate`` mixtures
are diagnostic records rather than the fitted three-component score.

Local molecular graph
---------------------

The graph contains atoms of a focal residue and its contacting residues, with
covalent bonds within each amino acid and observed focal-to-neighbour atom
contacts. Neighbour-to-neighbour contacts and peptide bonds between residues
are outside this local graph. For internally connected residues, m atom contacts
to q distinct partner residues add exactly **m-q** cycles: beta1=E-V+C and the
component count changes from q+1 to one. Intrinsic amino-acid rings cancel.
This exact contraction retains repeated contacts from the same physical atom.
One atom contacting two atoms of one neighbour adds one cycle; contacting two
separate neighbours once each adds none. No molecular-graph traversal is needed.

Export ``occupied`` distinct contacted atoms, ``pairs`` physical atom contacts,
``partners`` contacting residues and ``cycles=pairs-partners``, including zeros.
``available`` counts observed heavy atoms. Canonical heavy-atom completeness
is explicit; incomplete focal/partner residues abstain without reconstructed
atoms or inferred geometry. ``--atom-scope sidechain`` counts only contacts
incident to focal side-chain atoms; all amino acids remain internally connected.
The partner atom may be backbone or side chain. Glycine then has zero available
side-chain atoms. The default scope includes all heavy atoms.

Controls and empirical tails
----------------------------

**Retain every available mock/negative model outside the benchmark.** No cap,
downsampling or paired-CDR3 deduplication is applied. Controls may share TCRs
or epitopes with queries. Exclude only benchmark model IDs, coordinate hashes
and exact full observed TRA+TRB+epitope pairings. Full chain sequences distinguish
different receptors with matching CDR3s. Original CLI accounting supplies hashes.
Source metadata must identify every reference model and declare label zero;
``--reference-label 1`` instead selects a cognate reference.

Reference cells condition on chain role, loop/peptide region, amino acid and
observed atom count. For statistic x and N mock residues in its cell,
``p=(1+count(mock>=x))/(N+1)``. Ties are inclusive; x=0 gives exactly p=1;
beyond all observed values p=1/(N+1). Raw integer residue/model counts stay
explicit. There is no query-cohort fitting or selected residue position.

Optional ``--balance-epitopes`` retains every model and uses inverse epitope
frequency weights, normalized to mean model weight one. Then
``p=(1+sum(weight*I(mock>=x)))/(1+sum(weight))``. Integer counts and weighted
masses are exported separately. This is composition adjustment, not an assay
weight or a fitted score coefficient. Receptor counts are not capped or collapsed.

Residue tables retain probabilities, bits, positive/negative reference counts,
weighted masses and support. Unsupported cells remain null there and contribute
zero evidence in the full aggregation denominator. All six loops must be present.
The scalar summaries are an initial local-motif representation, rather than an
exact chemical-subgraph likelihood or a calibrated hypothesis test.

CLI workflow
------------

Use fresh ``export --quantity contact-atoms --interface tcr_peptide`` and
``export --quantity atom-inventory --interface tcr_peptide`` outputs at 5 Å.
Keep the inventory's original ``.accounting.tsv`` beside it.

::

   tcren residue-surprise-counts --atoms atoms.parquet --contacts contacts.parquet -o residues.parquet
   tcren residue-surprise-counts --atoms mock-atoms.parquet --contacts mock-contacts.parquet -o mock-residues.parquet
   tcren residue-surprise-fit --features mock-residues.parquet --query-features residues.parquet --metadata mock-metadata.tsv -o mock.json
   tcren residue-surprise-score --features residues.parquet --reference mock.json --statistic cycles -o scores
   tcren residue-surprise-evaluate --baseline baseline.parquet --scores scores/scores.tsv -o comparison

The count command also writes ``residues.identities.tsv`` with full chain/peptide
sequences and coordinate hashes. Fit reads these sidecars for both inventories;
combining feature inventories requires combining their identity sidecars too.
An optional ``--positive-reference cognate.json`` in scoring supplies enrichment;
``--enrichment`` selects it during evaluation. Ordinary evaluation uses rarity.

The baseline requires id,dataset,stratum,epitope,y,iptm,prepared,matched. CLI tables
own AP/AUROC, exact binder/nonbinder counts and all cohort/template-stratum AP
changes. ``all_labelled`` and historical ``observed_baselines`` populations are
reported separately. Optional saved R5/R7 comparators are historical records,
not candidates or new fitted models. Inspected panels remain internal development.

Complexity and scope
--------------------

For A inventory atoms, M contact pairs and R residues, columnar grouping/joins
use O(A+M+R) table storage. Sorted grouped as-of joins avoid a query-by-reference
cross product. Existing Polars/NumPy kernels suffice; no new native kernel or
worker pool is required. Empirical background frequencies are estimated
reference dependencies. Dependent residues and source-labelled mismatches do
not establish calibrated binding probabilities or independent evidence factors.

The graph cycle-space construction follows the standard
`cycle-basis definition <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.cycles.cycle_basis.html>`_.
Inclusive plus-one empirical ranks follow the construction in the
`Berkeley conformal prediction notes <https://statgrades.berkeley.edu/~ryantibs/statlearn-s24/lectures/conformal.pdf>`_;
the required exchangeability guarantee is not asserted here.

.. automodule:: tcren.contacts.residue_surprise
   :members:
