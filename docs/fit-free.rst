Reference dependencies and historical models
============================================

The default TCR screen uses a frozen confidence-dependent mixture of ipTM,
recognition topology and mock-calibrated contact enrichment. Its two coefficients
were selected during labelled development, and its local histograms, whole-mock
calibration and confidence representatives are estimated reference dependencies.
Nothing is fitted to the query, but the default is not a fit-free model.
See :doc:`screening` for the current formula and missing-value policy.

Fixed empirical contact potentials retain their own contact counts, pseudocounts,
weighting and structural/virtual reference conventions. Historical Gaussian,
GLM, fitted Potts and transition readouts remain distinct APIs. Their availability
does not make them components of the current screen.

Searched annotation references
------------------------------

ARDA 2.36.2 or later validates the source reference and actual MMseqs database
payloads by content, preserving legacy indexes as evidence. Source-root discovery
recognizes ARDA itself rather than an unrelated enclosing project.
``scoring``, ``features``, ``recognize`` and ``export`` retain ``searched_indexes``
in provenance sidecars: resolved database prefixes and SHA256 identities of the
searched sequence, header and index files, including split payloads. These are
captured around annotation, not inferred from package versions after scoring.
Mid-search changes and mixed generations across input chunks fail explicitly.
An empty identity record does not certify annotation provenance.

Contact sums and references
---------------------------

For a set of observed residue contacts C and a selected pair-potential table V,
the contact sum is Phi = sum over (i,j) in C of V(a_i,a_j). The identities a_i and
a_j are the observed amino acids. The virtual reference difference is
DeltaPhi = Phi(observed) - Phi(reference), with C unchanged. Lower values follow
the energy direction. These are dimensionless statistical-potential sums, not a
calibrated binding free energy. The manuscript reference uses structural alanine truncation: preserve backbone
and C-beta coordinates, remove the remaining side-chain atoms, construct C-beta
for glycine using the existing substitution geometry, and rebuild contacts.
This differs from the historical virtual reference and can change rankings.

``scoring --delta`` references peptide identities. ``dPhi_peptide_sum`` is the
unweighted TCR:peptide plus peptide:MHC reference sum. Historical ``Phi_total`` and
``dPhi_total`` use estimated native-reference scale weights and are excluded from
final scoring; they must not substitute for the unweighted sum. ``export`` references receptor
identities, or exports existing exact decompositions. No fitted replacement is
introduced:

.. code-block:: console

   tcren scoring -s structures/ --delta --reference-mode structural --no-superimpose --tcr-peptide-potential tcren2 -o peptide.csv
   tcren export -s structures/ --potential tcren2 --quantity receptor-reference --reference-mode structural -o receptor.tsv
   tcren export -s structures/ --potential tcren2 --quantity position -o positions.tsv
   tcren export -s structures/ --potential tcren2 --quantity residue -o residues.tsv
   tcren export -s structures/ --potential tcren2 --quantity contacts -o contacts.tsv
   tcren export -s structures/ --potential tcren2 --quantity placement -o placement.tsv

Use ``--reference-mode structural`` explicitly for the manuscript. The default
``virtual`` mode remains available for exact historical reproduction. Output rows
and provenance identify the mode; do not mix modes in a comparison. Structural
peptide references require alanine and use the selected cutoff and contact weighting.
The peptide reference contact map is built once per structure and reused across
both peptide interfaces. Potential sums then operate on the contact tables;
SASA and coordinate-based descriptors still require their structural calculations.
Use the command's ``--threads`` CPU budget and batch all structures in one call.
Do not wrap internally threaded annotation calls in an outer worker pool.

The receptor export substitutes alanine in all, CDR, CDR1/2 or CDR3 receptor
positions, separately for TCR:peptide and TCR:MHC and their explicit sum.
Position exports use ``position_profile``; residue exports use
``residue_importance``. Leave-one-out ``explain`` is not an equivalent decomposition.
Each contact contributes to both residue endpoints, so the sum of endpoint
contributions is twice the contact sum. Distances and placement are in angstroms.

Identity and accounting
-----------------------

Exports retain full source stems, physical chain IDs, sequential residue indices,
author residue numbers and insertion codes. CDR regions exclude conserved junction
anchors. One batched annotation is used per species plus one MHC annotation batch.
Files, directories, glob/manifests and supported structure archives are accepted.

Every input receives an accounting row, including parsing, annotation and export
errors. Empty/all-failed runs exit nonzero and keep diagnostics. Mixed runs retain
all failures for explicit acceptance by their consumer; a benchmark must verify
its full key set rather than silently discard failed rows.

Sidecars retain residue assignments, chain annotations, installed versions,
reference hashes, the exact potential tables and their hashes. A custom potential
has no inferred derivation or unit calibration. Source paths and annotations may
be private; review them before sharing generated outputs.

Historical assessment availability
----------------------------------

``score_table`` exposes a finite-value ``<readout>_available`` flag per readout,
``assessment_available`` for binder-score availability, and a JSON list of
``missing_descriptors``. Absent columns are marginalized; missing values in
present required columns remain missing. No imputation or refit is performed.
``inspect-models`` omits fitted assessment unless ``--historical-assessment`` is
explicitly requested.

Documentation builds require the ``docs`` extra and the external Pandoc executable
for nbsphinx. IPython supplies the notebook syntax lexer.

Molecular publication consumers
-------------------------------

Scene producers read exact ``position``, ``residue``, ``contacts``, ``placement``
and ``pmhc`` exports. They do not substitute leave-one-out ``explain`` for a
position decomposition or calculate new contact scores during drawing. Physical
chain, author residue number and insertion code identify selections; region
membership comes from annotation. Preserve raw-coordinate, searched-reference,
potential, camera and scene-geometry identities and check the displayed atom
colours against their authoritative rows. Inspect the final PNGs and composed
PDF pages. Keep class II, epitopes and template strata separate. Historical
images and manifests remain archived when labels or pixels change.

Structural class-II registers
-----------------------------

Sequence-only ``anchor_indices`` remains a prediction for historical callers.
For an experimentally verified core reference, export structural correspondence
in the same annotation batch::

   tcren export -s native.txt --quantity contacts --potential tcren2 --cutoff 8 -t 2 \
       --register-reference-id 1fyt --register-p1 3 -o contacts.tsv

The reference must be included in ``native.txt``. ``contacts.tsv.registers.tsv``
contains every forward/reverse nine-residue core, its C-alpha RMSD after aligning
both MHC groove chains, physical P1/P4/P6/P9 identities, groove RMSD, the exact rigid transform and the gap
between the two lowest core RMSDs (all distances in angstroms). No peptide atoms
enter the rigid fit and no sequence motif, trained coefficient or acceptance
threshold enters core ranking. A minimum is a correspondence candidate, not a
certificate of structural correctness; inspect residuals and the source structure.
Incomplete cores fail explicitly, with errors retained in input accounting.
Register exports retain known HETATM amino acids as ``X`` with their chemical
names and geometry; unsupported modified peptides are not energy-scored.
Historical parsing defaults are unchanged.

``--peptide-selection verified.tsv`` supplies source-verified antigen intervals:
``pdb.id``, ``chain.id``, ``first_pdb``, ``last_pdb`` (inclusive author numbering).
This handles engineered fused antigens in the presence of auxiliary peptides.
Original residue indices and coordinates remain unchanged; insertion-code
boundaries and overlap with annotated receptor regions are rejected. The selection and register
exports are hashed in provenance. Ordinary exports keep their previous behavior.
The Python APIs are ``tcren.refine.structural_register.structural_register`` and
``select_peptide``; they require already annotated structures and perform no search.

Local topology and topography of contacts
-----------------------------------------

``export --quantity contact-context --potential tcren2`` retains every original
TCR:peptide residue-pair contact, author identities and region markup. It adds
endpoint degrees and the number of bipartite four-cycles through each edge,
computed as ``(B B^T B)[u,v] - degree(u) - degree(v) + 1`` for an existing edge.
Here B is the binary incidence matrix of all TCR:peptide contacts, u is a TCR
residue and v a peptide residue. Framework contacts remain identifiable; these
are graph cycles, not geometric loops or spatial Betti numbers.

The closest peptide atom is projected onto the existing receptor-free pMHC
surface grid. Rows retain its height, the projected cell's peptide/MHC ownership,
surface height and KD hydropathy, and height minus the existing 90th-percentile
helix rim (angstroms). Cell ownership identifies peptide versus MHC, not a
particular peptide residue. Missing rims/out-of-window cells remain undefined.
``pair_energy`` is the supplied potential's exact per-pair contribution; missing
matrix entries remain NaN and are omitted only when reproducing the established
sum. Neither a zero retained sum nor a missing value is a binding probability.

This export supports inspection of graph-supported or protrusion-associated
contacts without learning thresholds from labels. It does not apply a filter,
replace the original score or establish that contacts outside a mask are false.
Keep assay targets, classes and template strata separate during evaluation.

``contact-summary -i contacts.tsv -o contacts.summary.tsv`` reduces that exact
stamped export and its accounting sidecar, with no parsing, annotation or search.
Canonical amino-acid-pair count matrices multiply the fixed potential vector:
identical contact compositions remain tied regardless of input row order. It
exports counts and sums for all contacts, edges participating in four-cycles,
peptide-owned projected cells above the established ridge, and their intersection.
These are declared structural selections, not calibrated nonbinding filters.
Accounting retains empty models; undefined topography remains NaN for affected
surface selections. Missing potential entries are counted and omitted as in the
original sum. ``--accounting`` supplies a relocated ledger explicitly. Compare
new sums with exact recorded baselines rather than perturbing old score ties.

Exact atom observations
-----------------------

.. code-block:: console

   tcren export -s native.txt --quantity contact-atoms --interface all -t 1 -o atoms.parquet

This opt-in quantity requires no potential. It selects CDR1/2/3 contacts with
peptide or annotated MHC helices and peptide contacts with those helices, counting
each physical heavy-atom pair once. Groove-floor and B2M contacts are outside this
observation scope. Missing helix annotation fails explicitly. Choose an existing
individual interface name instead of ``all`` to restrict the export.
The ``.parquet`` suffix selects compressed columnar output for this quantity;
other suffixes retain TSV output. Its summary/accounting sidecars remain TSV.

The per-residue conformer policy chooses the lexicographically first nonblank
alternate label plus blank atoms, preferring the selected label on atom-name ties.
This deterministic convention is not occupancy selection or proof of correlated
conformer compatibility between residues. Legacy contacts/scoring are unchanged.
Glycine contributes observed backbone atoms only; no C-beta is synthesized.

Physical sides retain atom/source indices, author residue/insertion, alternate
labels, coordinates, chain roles and regions. ``chem.a``/``chem.b`` sort whole
amino-acid/atom-name tuples, independently of physical side order. Observed
N/CA/CB geometry yields trRosetta omega, endpoint-specific theta/phi in radians
and C-beta distance in angstroms. Missing atoms and degenerate vectors yield null,
not a favourable geometry value. Existing v2 chemical flags remain distance/type
heuristics; stacking uses the existing ring geometry. Hydrogen-bond angles,
line-of-sight, plausibility probabilities and penalties are not calculated here.

The summary sidecar separates atom-pair and residue-pair frequencies by MHC class,
interface and receptor region. Residue rows use the closest retained atom distance
and the union of type flags across all atom pairs. Chemical types overlap; their
counts need not sum to the total. Summary groups contain observed contacts only;
the accounting sidecar separately retains successful zero-contact structures.
No count is labelled an independent sample size or a probability of physical
validity. Provenance retains the summary hash, conformer policy, annotation and
searched-reference identities. See :mod:`tcren.contacts.observations`.

Exploratory recognition graphs
------------------------------

``export --quantity recognition-topology`` exports seven fixed graphs: the
released CDR:peptide + CDR:MHC graph; CDR:peptide + peptide:MHC; the latter
with peptide sequence edges; and finally that graph with CDR:MHC edges.
Residues use physical chain/index identities and contacts use the supplied
heavy-atom cutoff (default 5 Angstrom). Sequence edges connect consecutive
observed peptide residues only with a C--N distance in (0, 1.8] Angstrom;
missing atoms and breaks are counted. No other backbone edges are added.

The two additional variants replace sequence edges with spatial peptide
side-chain contacts: ``tp_pm_sidechain`` (G2a) and ``tp_pm_sidechain_tm``
(G3a). A unique residue pair qualifies when at least one pair of side-chain
heavy atoms is within the supplied cutoff. N, CA, C, O and OXT are backbone;
hydrogens are excluded. Adjacent residues may qualify spatially, but sequence
adjacency alone creates no edge. Glycine has no side-chain heavy atom; missing
side chains are not synthesized. ``peptide_sidechain_edges`` records the count.

The ``tp_pm_tm`` variant adds CDR:MHC contacts directly to G1 without
peptide sequence or side-chain edges.

For each graph, beta1 = E - V + C and cyclo_frac = beta1/E, where E counts
unique edges, V incident vertices and C connected components. Empty graphs
have undefined cyclo_frac. Receptor contact availability is reported separately.
Cycles express redundant paths, not contact energy or the minimum cut that
would separate receptor from ligand. Internal peptide:MHC cycles and covalent
sequence edges are not equivalent to receptor recognition contacts. Neither
scalar alone assigns peptide anchors or register. The released score is unchanged.

Simplicial contact homology
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The recognition export also constructs each graph's clique complex and reports
``simplicial_beta1`` and ``simplicial_beta2`` over GF(2), together with triangle
and tetrahedron counts and the two boundary ranks. These are distinct from the
unchanged graph ``beta1`` and ``cyclo_frac``. No score default changes.

For each sorted edge (a,b), intersect the forward neighbor sets of a and b to
enumerate every triangle (a,b,c) once. Intersect the three forward neighbor sets
to enumerate every tetrahedron (a,b,c,d) once. The transposed boundary matrices
have a row for each simplex and a column for each codimension-one face, with
entry one exactly when that face belongs to the simplex. Existing vectorized
GF(2) Gaussian elimination computes ranks r2 and r3. Over GF(2), each face of a
boundary occurs twice, so consecutive boundary maps compose to zero. Thus::

    simplicial_beta1 = E - V + C - r2
    simplicial_beta2 = F - r2 - r3

Here F counts triangular faces; E, V and C retain their graph definitions.
The 3-skeleton suffices: higher simplices do not change these two homology
groups. A filled tetrahedron therefore has beta2=0; retaining triangles but
omitting tetrahedra would incorrectly report a shell. This is exact algebra,
not a geometric cavity reconstruction. All inter-entity and peptide-edge
policies remain those of the selected graph.

Implementation: ``topology.graph._clique_homology`` reuses
``topology.footprint._gf2_rank``. Boundary matrices are allocated sequentially
as uint8 arrays, with a 64 MB per-matrix ceiling and an explicit failure above
it. For q rows and p columns, dense elimination takes O(min(q,p)*q*p) time and
O(q*p) memory in the worst case. Edge, triangle and tetrahedron enumeration
stores O(V+E+F+T) objects, where T counts tetrahedra; clique enumeration is
output-sensitive and can grow combinatorially. No cutoff sweep, approximation,
fitted direction or persistence filtration is performed.

Exploratory normalization uses beta1/E and both beta2/E and beta2/F, setting
beta2/F=0 when F=0 (then beta2 is necessarily zero). No-edge fractions remain
undefined. Unlike beta1/E and beta2/F, beta2/E is not generally bounded by one.
Any fixed ipTM combination must record its chosen denominator and retain the
unchanged released score as a comparator. Homology counts describe unfilled
cycles; they do not measure binding energy, force or recognition strength.

The prepared 4.4.0 screen uses ``tp_pm_tm``. The ``released`` graph identifier
retains the historical 4.3 graph definition for reproducible comparisons.

Energy-weighted recognition graphs and persistence
--------------------------------------------------

``export --quantity weighted-recognition-topology`` uses the same seven graph
variants as ``recognition-topology``. It adds diagnostic descriptors only; the
screen score and frozen reference parameters are unchanged. Install the optional
``topology`` extra for ``export --quantity persistence``.

The graph vertices are incident physical residues, keyed by chain and observed
residue index. The contact cutoff remains 5 Å. Edge goodness is
``w = max(-energy, 0)``: TCRen2 is directed TCR→peptide and MHC→peptide;
MJ is directed TCR→MHC. No fitted interface coefficients are applied. Units are
the respective dimensionless table scales, which are not thermodynamic energies.
Peptide internal edges have no specified energy. An unsupported edge remains in
the unweighted graph and increments ``missing_weight``; all weighted descriptors
of that graph are undefined. Empty graphs and graphs lacking CDR contacts retain
explicit availability/count fields.

``w_total``, ``w_mean`` and ``w_max`` reduce edge goodness. Node strength is
``s_i = sum_j w_ij``; ``s_mean``, ``s_max`` and ``s_even`` summarize its distribution,
with Pielou entropy normalized by the number of incident vertices. Density is
``2E/(V(V-1))`` for E edges and V vertices. ``w_cyclo`` is goodness outside a
maximum-weight spanning forest divided by total goodness. It is undefined if
total goodness is zero. It measures redundant strength under that forest
convention, rather than assigning an edge uniquely to a cycle-space basis.

``giant_frac_q25/q50/q75`` retain edges at or above each model's own goodness
quantile and report the largest component's vertex count divided by V, retaining
isolated vertices. Quantiles are within one graph, never across a cohort.
``graph_energy`` is the sum of absolute eigenvalues of the symmetric weighted
adjacency. ``graph_energy_c`` subtracts the same-chain induced-subgraph energies
and divides by the number of inter-chain edges. ``graph_energy_bounded`` divides
by twice total goodness; the triangle inequality for matrix nuclear norms gives
the bound [0,1].

``clustering`` averages the fraction of neighbor pairs that are adjacent, with
zero for degree below two. ``transitivity`` divides summed per-vertex triangle
incidences by summed unordered neighbor-pair counts (zero when none exist).
``mincut_tcr_pmhc`` uses unit-capacity maximum flow with a super-source connected
to every CDR vertex and a super-sink to every pMHC vertex. Consequently, this
all-terminal definition removes every crossing edge; it is not a single-pair
bottleneck. ``mincut_frac`` divides by the smaller side's summed graph degrees.
The weighted variant rounds ``100*w`` to integer capacities, reports the cut in
original goodness units, and normalizes by the smaller summed rounded degree.
Zero denominators are undefined.

``alg_conn_w`` is the second eigenvalue of the normalized Laplacian on the
largest positive-weight component, requiring at least three vertices. Its range
is [0,2]; ``alg_conn_w_bounded`` divides it by two for mixtures requiring [0,1].
The existing unweighted connectivity API is unchanged. Betweenness uses
reciprocal-positive-goodness path lengths, excludes zero-weight edges and divides
ordered source/target dependencies by ``(V-1)*(V-2)``. Shortest-path dependencies use native triangular solves. Shortest-path ties use
relative/absolute comparison tolerance 1e-12 in the dependency reduction.
``btw_tcr`` and ``btw_pmhc`` average over the respective incident sides, including
zero-strength vertices. Weighted participation replaces degrees with strengths
in ``1-sum_s (strength_to_module_s/strength_total)**2``. Modules follow the
existing literature API: peptide/MHC for TCR vertices and six CDR loops for
pMHC vertices. Unsupported zero-strength participation is excluded from its mean.

Persistence uses the Cα coordinates of endpoints of the ``tp_pm_tm`` graph and
a companion Cβ cloud (Cα fallback for Gly only). Missing atoms and fallback counts
are explicit. GUDHI builds the exact Vietoris–Rips 2-skeleton through 12 Å over
GF(2). For H0 and H1, the export reports finite-bar counts, total/max/mean lengths
in Å and normalized entropy, plus censored infinite-bar counts. Zero/one-bar
entropy and absent max/mean are undefined. Infinite intervals never receive
fabricated death times. Total length divided by ``V*12 Å`` (V is the resolved point count for that cloud)
and maximum length
divided by 12 Å are exported separately; the former is not a universal H1 bound
and must pass a [0,1] eligibility check before entering any bounded mixture.
No H2 score is exported. A missing GUDHI import raises explicitly.

The implementation reuses native contact searches, SciPy graph primitives and
NumPy dense linear algebra; no new native kernel is needed. With V vertices,
E edges and F triangles, adjacency/eigenvalue work uses O(V²) memory and O(V³)
time. Persistence can require O(V³) simplices in the 2-skeleton. These are
small interface graphs, rather than whole-protein atom clouds. These descriptors
have no calibrated affinity, strain or nonbinding interpretation.

.. autofunction:: tcren.topology.weighted.weighted_graph_features

.. autofunction:: tcren.topology.persistence.persistence_features

Explicit MHC matrix weights and per-residue normalization
---------------------------------------------------------

For ``export --quantity weighted-recognition-topology``, ``--potential`` sets
TCR-to-peptide weights, ``--mhc-potential`` sets TCR-to-MHC weights, and the optional
``--pmhc-potential`` sets peptide-to-MHC weights. Each supplied asymmetric matrix
is looked up in its biological role direction before assigning an undirected
graph edge weight ``max(-energy,0)``. When no peptide-MHC matrix is supplied, the
historical MHC-to-peptide lookup in the TCR-peptide matrix is retained. Every
used matrix is exported and hashed in provenance; explicit inputs are honored.

Let E be the number of physical residue-contact edges, V the number of distinct
incident residues, and W the sum of their nonnegative goodness weights.
``edges_per_vertex=E/V`` is half the mean residue degree;
``w_per_vertex=W/V`` is half the mean weighted strength. Their fixed bounded
forms are ``E/(E+V)`` and ``W/(W+V)``, respectively. These require no all-pairs
contact assumption or fitted maximum degree. V includes endpoints of zero-weight
edges. Empty graphs remain undefined; missing matrix entries retain physical
edges and withhold energy-dependent descriptors. Historical ``density`` remains
available for compatibility and has a separate all-pairs definition.

When E denotes the sum of edge weights, the requested E/V statistic is
``w_per_vertex`` (W/V here), rather than ``edges_per_vertex``.

Boltzmann graph weights
-----------------------

``export --quantity weighted-recognition-topology --graph-weighting boltzmann
--graph-beta 1`` uses ``exp(-beta*E)`` for each finite supplied matrix entry E.
Beta must be finite and positive; beta=1 uses the dimensionless matrix scale
without temperature calibration or fitting. Neutral entries have weight 1;
unfavorable entries retain positive weight below 1. Missing entries remain
unsupported. Overflow or underflow raises an error rather than clipping.
The historical positive transform remains the default. Graph definitions,
W/V normalization and score mixtures are unchanged. Both weighting and beta
are recorded in export provenance.
