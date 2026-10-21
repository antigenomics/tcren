Experimental contact penalties
==============================

The four declared penalty constructions share one frozen geometry compatibility
weight. They are experimental development readouts, not validated probabilities
of natural contacts, nonbinding or affinity. Existing scoring defaults are unchanged.

Raw reference and freezing
--------------------------

``contact-prior-geometry --source DEPOSIT --exclude native_ids.txt -o NEW_DIR``
checks the deposited antibody contacts against current compressed raw mmCIFs,
one structure at a time. It requires exact amino-acid/atom identity and float32
coordinate agreement with a unique selected atom. IMGT and author residue numbers
are not equated. First model and deterministic single-conformer selection apply;
ambiguous or unmatched atoms retain missing recovered frames. The output retains
the deposited identities, current author identities and per-structure accounting.
``--limit`` explicitly produces a diagnostic subset.

``export --quantity contact-atoms --interface all`` supplies the fresh native
observations. Exclude every native PDB ID from the external reference first.
``contact-prior-fit --external external_geometry.tsv.gz --native native.parquet
-o NEW_RESOURCE.json`` freezes both reference sufficient statistics. No benchmark
labels, predictor confidence, fitted classifier or cohort normalization enters.
The source deposit and current-coordinate correspondence have different provenance;
their hashes and overlap exclusions must accompany a derived resource.

The distance transform is :math:`x=\log[d/(5-d)]`, where :math:`d` is the observed
heavy-atom separation in angstroms and :math:`0<d<5`. Orientation marginals use
:math:`x=\operatorname{atanh}(\cos a)` for the observed omega, theta and phi angles
:math:`a` in radians. Exact transform boundaries and missing frames stay unavailable;
there is no clipping floor or invented C-beta. Torsion sign is discarded. A normal
model for these transformed marginals is an approximation, not a joint angular
density or a demonstration of normal residuals.

Distances condition on unordered whole amino-acid:atom tuples; pooled chemistry
uses element plus backbone/sidechain identity. Orientations retain ordered binder-to-target residue pairs. Theta and phi
also identify binder or target endpoints, including separate fallback pools.
Active resources use contact-marginals-2; version-1 angular resources are rejected. Each physical
residue edge supplies orientation observations once, rather than once per atom pair.
An external exact group with fewer than two observations or zero variance uses
the declared external pool. A scoring category unseen during freezing uses that
same external pool, with no scoring-time update.

For each marginal, the external-data posterior under the reference prior
:math:`p(\mu,\sigma^2)\propto1/\sigma^2` becomes the native-update prior. Its
effective observation count is the actual external count, not a tuned strength.
The resource records external and native counts, means and squared-deviation sums
separately. Contacts within a complex are correlated; these counts are not numbers
of independent complexes. Pools overlap their exact groups and are alternatives,
not extra observations added to those groups.

Let :math:`n>1` be the resulting observation count, :math:`\bar x` its mean and
:math:`B>0` its sum of squared deviations. The dimensionless predictive scale is
:math:`h=\sqrt{B(1+1/n)/(n-1)}`. Compatibility is
:math:`s=2\Pr[T_{n-1}\geq |x-\bar x|/h]`, where :math:`T_{n-1}` is a standard
Student-t variable. This tail is dimensionless and lies in :math:`[0,1]`;
it is not a probability of physical naturalness. No product of marginal tails or
marginal-independence assumption is made.

The categorical Dirichlet update has external integer counts as its alpha values
on the observed external support. Native categories outside that support remain
explicitly unavailable; no smoothing mass is invented. These conditional frequency
summaries do not determine the geometry weight or label rare chemistry unnatural.
Reference consistency tables are descriptive in-sample checks, not held-out
calibration or independent posterior-predictive validation.

Residue contributions and comparison
------------------------------------

``contact-penalties -s structures.txt --resource RESOURCE.json --chunk 512 -t 0
-o NEW_DIR`` preserves the original all-TCR peptide contact sum for both TCRen2
and the bundled classic TCRen table. Its ``karnaukhov2022`` alias retains the
May 2026 update; it is distinct from the recovered original-paper table (see
``SOURCES``). Coefficients and defaults are unchanged. One bounded structure batch is parsed and
annotated at a time; searches stay in the parent, followed by coarse contiguous
spawn-worker scoring and ordered reduction. The default is all allocated CPUs
(``-t 0``) and a 512-structure batch; ``-t 1`` requests serial execution and
positive N selects an explicit bounded CPU budget. Numerical kernels remain
serial inside workers, and pool failures propagate without serial fallback.
Residue identities, searched-index lineage, potential tables and all failures are
exported. CDR1/2/3:peptide contacts receive the experimental increment; framework
contacts retain their original contribution and are explicitly unassessed.
The separate pMHC observations inform the reference but do not silently add another
interface energy to the legacy TCR:peptide baseline.

For physical residue edge :math:`r`, set :math:`u_r=1-\min s`, taking the minimum
over its supported atom-level and orientation marginals. Thus :math:`0\leq u_r\leq1`.
Unavailable marginals are counted; an edge with no support remains unavailable.
This worst-supported reduction gives large residues more opportunities for a
penalty. It was declared before accuracy inspection and is not scanned against labels.

Let :math:`e_r` be the original dimensionless potential entry and
:math:`M=\max_{a,b}|V_{ab}|` over finite entries of the full frozen table :math:`V`.
The four resulting contributions are:

* Fixed: :math:`e_r+u_r`, one natural-log potential unit at full weight.
* Matrix-scaled: :math:`e_r+u_rM`.
* Literal flip: :math:`(1-2u_r)e_r`, a signed diagnostic that can reward a contact
  whose original entry was positive.
* Favourable-only flip: :math:`e_r+2u_r\max(-e_r,0)`.

The other three increments are nonnegative. A zero table has :math:`M=0`; the
matrix-scaled arm is then degenerate. Missing table entries remain null. The
legacy scorer's omission of unknown entries is retained and counted explicitly;
such totals cover only defined entries. No absent entry is assigned zero or used
to define :math:`M`. Every physical edge exports its original contribution,
weight, increment and resulting contribution. The sum counts a residue edge once.

``contact-compare --labels labels.tsv --scores NEW_DIR/scores.tsv -o COMPARISON``
reports ROC AUC, average precision and the existing precision-at-recall metrics.
All four arms, original energy and native ipTM use one common finite population
within each potential, cohort, MHC class and template-coverage stratum. It writes
row membership, excluded counts and every delta/loss against both baselines.
Labels need unique ``id``, ``dataset``, ``task``, ``cohort``, ``mhc_class``,
``template_coverage``, binary ``y`` and native ``iptm``. Lower energies and higher
ipTM remain favourable. No significance tests, fitted coefficients, threshold scan
or label-selected direction is introduced.

API
---

.. automodule:: tcren.contacts.plausibility
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: tcren.contacts.prior_geometry
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: tcren.contacts.penalty_export
   :members:
   :undoc-members:
   :show-inheritance:


Mechanism audit without changing scores
---------------------------------------

``contact-geometry-audit`` accepts the same raw-source/resource/chunk options as
``contact-penalties`` but exports only CDR:peptide chemistry and geometry lineage.
It records any-atom chemical support, closest-atom support, overlapping v2 flags,
the minimum-tail family/key/parameter source, distance range and distinct marginal
opportunities. A winning orientation belongs to the residue frame, not a unique
atom pair. Existing v2 hydrogen-bond/cation–pi labels are distance/type heuristics,
not angular energetic validation. No scores or reference parameters are updated.

``contact-mechanism --run ACCEPTED_RUN --geometry NEW_DIAGNOSTICS -o NEW_AUDIT``
reads the accepted final decompositions and unchanged labels/scores. It requires
complete row membership and verifies physical keys, frozen weights, four signed
sum identities and every recorded AUROC delta by exact pair-concordance accounting
(ties carry half credit). It exports model/label profiles, chemistry/marginal
classes, motif sums, matrix-cell selection counts and illustrative changed-rank
pairs, separately per cohort/class/template stratum and potential.

Chemistry-supported and proximity-only counts are heuristic labels, not physical
truth. Continuous geometry equivalents are sum(1-u) and sum(u) across assessed
physical edges, not integer real/virtual classifications. Missing potential entries
are retained in contact counts and omitted from the original energy arithmetic,
with a separate known-entry geometry sum. Framework contacts remain unassessed.
The original potential coefficients are unchanged; filtering or weighting changes
the score and generally changes the contact population's log-odds interpretation.
No new metric, binary cutoff, model selection, fitting or default change occurs.

Multiplicity-conditioned weighted TCRen
---------------------------------------

``contact-weight-calibrate --external external_geometry.tsv.gz --native native.parquet
--resource marginal-resource.json -o calibration.json`` freezes one empirical
minimum-tail distribution per reference source and opportunity stratum. The native
calibration selects CDR:peptide contacts; antibody:antigen observations supply only
strata absent from native contacts. Other native interfaces do not enter these ranks.
The underlying native/external marginal resource is unchanged and its SHA256 is
bound into the new resource. Neither command accepts benchmark labels.

For each physical residue edge, let ``s`` be its worst supported marginal predictive
tail. Its stratum records the number of supported atom-distance rows and the number
of distinct ``(marginal key, transformed value)`` omega, theta and phi terms. Repeated
angular copies across atom rows contribute once. If the selected frozen stratum
contains ``N`` reference edges, the contact weight is
``w = (1 + count(reference minimum tail <= s)) / (N + 1)``. The added one is the
query's own inclusive rank contribution, with ties retained conservatively; it is
not a tuned floor or prior-strength parameter. Native strata take precedence over
external strata regardless of count. Unseen strata or unavailable geometry retain
``w=1`` with explicit status, rather than dropping the model. Every scored edge
exports reference-edge/structure counts and the selected source.

``contact-penalties -s manifest.txt --resource marginal-resource.json
--calibration calibration.json --chunk 32 -t 1 -o NEW_DIR`` adds weighted TCRen to
both unchanged matrix tables. Each assessed CDR entry ``e`` contributes ``w*e``;
framework contacts retain their original contribution. Unknown matrix entries keep
the legacy omission with explicit counts. Empty contact maps retain zero scores.
``contact-compare`` detects the weighted column and compares original, weighted
and native ipTM on the same finite model IDs, within each cohort/class/template
stratum. Score directions remain fixed. Prior four-penalty comparisons are unchanged
when no calibration is supplied.

``contact-weight-summary --run NEW_RUN --accepted ACCEPTED_RUN --mechanism
PREVIOUS_MECHANISM -o NEW_AUDIT`` validates unchanged raw identities, labels,
contact graphs, matrix entries and original scores. It exports signed magnitudes,
CDR-only retained/downweighted equivalents, native/external/uncalibrated support,
chemistry classes, motifs and exact positive-negative rank-credit accounting.
All model denominators include zero-contact models. The mechanism input supplies
previous final chemical classifications only; old geometry weights do not enter
the weighted score.

These are in-reference conditional empirical ranks, not contact-validity,
occupancy or binding probabilities. Conditioning accounts for the declared
opportunity counts, without asserting independent atoms or exchangeability across
chemistries. Reference edges helped estimate the underlying marginals, so neither
independent validation nor conformal coverage is claimed. Matrix coefficients are
numerically preserved; their weighting creates a different score, not an empirically
re-estimated conditional potential or a binding free energy. All gains and losses
remain development evidence on previously inspected benchmark labels.


Directed reference-partition correction
---------------------------------------

Residue counts, potentials and endpoint marginals retain binder-to-target
direction. No reverse residue counts or transpose averaging are used. Only
individual atom-pair chemistry is exchangeable. Internal contact-table endpoint
views are selected by biological role before counting or scoring.

Freeze a native class-I CDR-to-peptide reference partition before evaluation::

    contact-normalizer-fit --native native.parquet --resource prior-v2.json --mhc-class I --out normalizer.json --reference-checks checks.tsv

Exact opportunity strata contain supported distance and distinct omega/theta/phi
term counts. Let logZ=log(mean_reference(exp(-u))). The corrected increment
g=u+logZ satisfies mean_reference(exp(-g))=1 in the empirical group; a constant
reference penalty cancels exactly. No assay coefficient, threshold, floor or
cohort centering is fitted. The native partitions and external marginals remain
estimated dependencies. This relative factor may have either sign; it is not
an absolute nonnegative punishment, a contact-validity probability or a binding
free energy.

Add --normalizer normalizer.json to contact-penalties; contact-compare includes
the normalized arm alongside the original and four diagnostic arms on matched
finite IDs. Resource hash, marginal version and MHC class must match. Unseen
groups get zero increment with a flag; missing geometry stays unavailable.
Framework contributions and legacy omitted potential entries remain unchanged.

Neither composition nor geometry-partition normalization supplies local SASA,
pocket enclosure, water thermodynamics or all noncontact opportunities. Total
SASA is not a residue-resolved opportunity distribution.

Atom-resolved context exposure
------------------------------

``tcren export -s MANIFEST --quantity contact-exposure -o OUTPUT`` exports
heavy-atom SASA in square angstroms for isolated TRA/TRB, intact pMHC, and the
bound complex. The target keeps its MHC and beta-2-microglobulin occluders.
Each row preserves physical atom/residue identity and binder/target role,
including noncontact atoms. ``partner_burial_A2`` is isolated minus bound area;
residue exposure is the sum over its atoms. Peptide-only
``peptide_alone_sasa_A2`` and ``mhc_shielding_A2`` distinguish MHC shielding
from receptor burial. Other atoms have null values in these two fields.

The calculation reuses Bio.PDB Shrake-Rupley with 100 sphere points and a
1.4 angstrom probe in the supplied coordinate frame. Finite surface sampling
has orientation-dependent quadrature error. The coordinates are not relaxed;
removing a partner does not recover an experimental apo conformation. Solvent
waters are not modeled, and these areas are not hydration energies.
No score or competing-contact allocation is changed by this export. Process
bounded input chunks because the general exporter retains its output tables.

Conserved surface allocation and exposed reference
--------------------------------------------------

``tcren export --quantity surface-contacts`` exports directed atom pairs with
separate ``binder_area_A2`` and ``target_area_A2``. Its ``.atoms.parquet``
sidecar contains every selected heavy atom, including noncontacts. The native
SASA kernel uses its existing 100-point Fibonacci-y sphere, a 1.4 angstrom
probe and Bio.PDB atomic radii. This orientation differs from the Bio.PDB
reference export above; do not silently exchange the two finite grids.

An atom's patches buried within its own partner are excluded. Each remaining
patch covered by the opposite partner is assigned to the expanded sphere
with smallest signed surface gap (distance to centre minus expanded radius).
Exact ties share its area equally. Allocated area sums to isolated minus
bound SASA for each atom. The two endpoint surfaces are separate quantities,
not a symmetrized residue statistic. This is a geometrical allocation of a
finite surface budget, not a chemical-valence or hydrogen-bond assignment.

``contact-penalties --surface-context`` additionally computes a fixed
exposed-surface reference, retaining all original arms. Let ``S_B(a)`` be
isolated exposed CDR area of binder identity ``a`` and ``S_T(b)`` peptide
area of target identity ``b`` within intact pMHC. Normalize each role's area
composition separately. Their product, conditioned on finite ordered matrix
support, is ``Q_ab``. For the supplied matrix ``V_ab``, compute
``logZ = log(sum Q_ab exp(-V_ab))``. Then ``V_ab + logZ`` is the negative
log density ratio against that declared reference on its positive support;
``sum Q_ab exp(-V_ab-logZ) = 1``. A constant table offset cancels.

The new ``surface_context_sum`` is ``sum A_ab*(V_ab+logZ)``, where ``A_ab``
is binder-side area assigned from CDRs to peptide. MHC competes for patch
ownership but is excluded from this CDR-to-peptide chemistry score. Also
export ``surface_raw_sum=sum A_ab*V_ab`` and the supported area, so the
normalization effect is exactly ``A_supported*logZ``. Units are square
angstroms times dimensionless table scale: area-weighted composite log odds,
not binding free energies or independent-contact likelihoods. Lower is favorable.
No coefficient, label-dependent selection, pseudocount or floor is fitted.

Missing matrix cells and observed area outside positive exposed support are
omitted with explicit area accounting. Empty supported sums are zero; the
per-area diagnostic is null at zero area and is not used to shrink the matched
comparison population. The complete atom and pair products are bounded-batch
Parquet shards. Every source structure remains in accounting. Global exposed
compositions do not capture every local pocket motif or hydration effect;
independent predictive validation remains necessary.

For ``N`` atoms, ``P`` sphere points and ``K_i`` local neighbours of atom ``i``,
point testing costs ``O(P*sum K_i)`` after spatial indexing, plus sorting and
unique neighbour bookkeeping. Memory is ``O(N+P+max K_i+E)`` for ``E`` output
atom edges; no dense all-pair coordinate array is allocated. In dense adversarial
geometry the time can still reach ``O(P*N^2)``. The kernel releases the GIL and
uses one native thread. Bounded source batches limit retained table/coordinate
memory; they do not eliminate output-size cost. Finite sampling gives exact
conservation on the sampled grid, not an exact continuous-surface error bound.

Accessible-boundary ownership option
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``contact-penalties --surface-context --surface-ownership accessible_boundary``
uses the same native sphere grid and isolated/bound SASA as the default
``deepest_sphere`` method. Each lost patch is assigned to the nearest sampled
isolated boundary point of the opposite partner; exact distance ties divide area
equally among distinct atom owners. Peptide and MHC boundary points compete before
CDR→peptide selection. Endpoint areas remain separate and residue roles remain
directed. Every recipient atom has positive sampled isolated exposure, and each
endpoint's total allocated area still equals isolated minus bound SASA.

The implementation exports exposed points from the existing native sampler and
uses SciPy's batched ``cKDTree`` with one worker. Only exact ties require ragged
processing. It raises an error if lost patches have no sampled opposite boundary.
Nearest sampled boundary is a geometric approximation, not an exact molecular
surface, a chemical valence rule or evidence of predictive improvement. This
option changes ownership and therefore score support, but neither fits weights
nor changes the default scorer. Compare separately named runs on identical rows;
never overwrite a prior run. The method is recorded in export provenance.
Directed atom-count layer
-------------------------

``contact-atom-count-fit --native observations.parquet --external curated.parquet
--equivalence-resource marginal-resource.json -o atom-count.json`` derives one
frozen directed residue-conditional atom-contact enrichment table. Native input
is restricted to class-I CDR:peptide; external input to the existing curated
peptide-target rows. The optional resource supplies only its explicit chemical
equivalences, applied identically to references and queries. Physical atom pairs
remain separate observations. No assay labels or weights enter derivation.

For binder atom type ``(a,i)`` and target type ``(b,j)``, pooled native/external
integer count ``n`` and residue-pair total ``N_ab`` give
``D = -ln[(n/N_ab)/(qB(i|a)*qT(j|b))]``. Role-specific atom fractions ``qB`` and
``qT`` describe observed-contact endpoint composition across opposing residue
types. Counts and sources remain explicit. This enrichment retains geometry,
availability and context effects; it is not isolated chemical energy.

``contact-penalties ... --atom-count-resource atom-count.json`` adds the mean
``D`` across observed atoms once per known CDR residue edge, with fixed unit
coefficient, to each original potential. Framework and missing original entries
retain their existing policies. Unseen atom cells receive zero correction with
support counts and stay in the averaging denominator; zero does not establish
compatibility. Missing CDR observations fail explicitly. Empty contacts retain
zero. ``contact-compare`` includes the ``atom_count`` arm on matched finite rows.

The dimensionless edge average controls direct atom-count inflation; it is not a
joint likelihood, a binding free energy, or a calibrated validity probability.
The neutral unseen-cell policy does not complete a globally normalized density.
References are estimated dependencies; raw contacts are correlated within
structures. No score-time fitting or default change is introduced.

Symmetric chemical core with directed composition
-------------------------------------------------

``contact-chemical-core-fit --native observations.parquet --external curated.parquet
-o chemical-core.json -t 0`` freezes the experimental ``chemical-core-count-1``
resource from existing native class-I CDR:peptide and curated antibody
peptide-target atom contacts at 5 Å. Each physical atom pair supplies one integer
count. Native and external contacts are pooled with their actual counts; no assay
labels, source weights or fitted score coefficients enter the estimator.
``-t 0`` uses one allocation-sized Polars pool in a fresh table-only subprocess.
BLAS/OpenMP kernels stay serial, with no outer process pool. This lets the
command own its immutable Polars pool without changing other commands' runtime.

The potential-specific mapping retains the existing chemical signature, splitting
carbonyl carbon sites (backbone C and D:CG, E:CD, N:CG, Q:CD) from other carbons,
separating methionine thioether sulfur from unknown-state cysteine sulfur, and
assigning histidine ring nitrogen an unknown protonation/tautomer state. The
accepted reference has 19 types, hence 190 possible unordered cells. The saved
AA:atom-to-type map applies identically to references and queries. The general
contact classifier is unchanged. Remaining coarse chemical merges and terminal
state assumptions persist; these labels do not assign charge, hydrogen-bond
satisfaction or a force field.

Estimator and normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let ``n_cd`` count directed binder-type ``c`` to target-type ``d`` contacts,
``N = sum_cd n_cd > 0``, and ``b_c = sum_d n_cd/N`` and
``t_d = sum_c n_cd/N`` be the two normalized endpoint compositions. For unordered
pairs, define::

    m_cd = n_cd + n_dc,  Q_cd = b_c*t_d + b_d*t_c  for c < d,
    m_cc = n_cc,         Q_cc = b_c*t_c.

The sums over unordered cells are ``sum m = N`` and ``sum Q = 1``. Diagonal
contacts are counted once. On observed cells define the symmetric core
``W_cd = m_cd/(N*Q_cd)`` and ``u_cd = -ln W_cd``. Do not average directed log
energies. A positive observed count with zero reference mass is inconsistent
input. On zero-count cells use ``W_cd = 1`` and ``u_cd = 0`` explicitly, including
zero-reference-mass cells, without evaluating ``0/0``. This is a declared neutral
core completion, with no pseudocount or inferred compatibility.

For binder residue identity ``a`` and target residue identity ``b``, freeze the
separately normalized conditional endpoint frequencies ``qB(c|a)`` and
``qT(d|b)`` from the same directed contact population. Then::

    Z_ab = sum_c,d qB(c|a) * qT(d|b) * W_cd
    p(c,d|a,b) = qB(c|a) * qT(d|b) * W_cd / Z_ab
    Delta(c,d|a,b) = -ln W_cd + ln Z_ab.

The sum defining Z is over ordered type pairs. A complete positive finite W and
normalized conditional endpoints imply ``0 < Z_ab < infinity`` and
``sum_c,d p(c,d|a,b) = 1``. On positive conditional reference support, Delta is
exactly ``-ln[p/(qB*qT)]`` for this constructed composition model. Symmetric W
is an exchange-invariant interaction assumption; it does not fit the observed
directed frequencies exactly or exchange biological roles. Directed conditional
compositions can give ``Z_ab != Z_ba``.

For fixed residue identities, write ``q = qB*qT``. On their positive support,
``E_q[Delta] = KL(q || p) >= 0`` and ``E_p[Delta] = -KL(p || q) <= 0``.
Here ``E_r`` denotes expectation under r and
``KL(r || s) = sum r * ln(r/s)`` on their shared positive support.
These are identities for the constructed conditional distributions. The endpoint
product q is not a distribution of mock structures or unbound atom opportunities,
so these identities do not establish real/mock classification. The log partition
is constant within one residue pair; its sum over a structure varies with the
contacted residue pairs. Removing it changes the model rather than choosing an
energy zero.

If U is the set of zero-count unordered cells, neutral completion gives
``sum QW = 1 + sum_U Q``. A separate global normalization is unnecessary: scaling
every completed W, including unseen cells, by the same positive constant cancels
in Delta. Scaling only observed cells and resetting unseen W to one changes the
estimator. The declared enrichment/completion gauge is stored directly.

Application and missing support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``contact-penalties ... --chemical-core-resource chemical-core.json`` exports a
separate ``chemical_core`` arm; ``contact-compare`` retains it alongside original
scores and other requested arms on matched rows. On each physical CDR1/2/3 residue
edge with finite original potential support, add the arithmetic mean of Delta
over all observed atom pairs with coefficient one. Atoms sharing a chemical type
remain separate physical observations. Both original residue potentials receive
the same increment. Framework and missing-original-entry policies are unchanged;
an eligible scored CDR edge missing all atom observations is an error.

An unseen chemical pair with positive conditional endpoint support receives
``Delta = ln Z_ab``, which can be nonzero. An unmapped query label, zero conditional
endpoint support or unavailable residue normalizer instead abstains with final
zero, explicit support counts and an unchanged atom-pair denominator. These
fallbacks make no density-ratio claim outside reference support. Unseen chemistry
and unsupported queries have distinct counts; no infinite energy is substituted.
Outputs retain edge increments, physical-pair counts and resource identity.

Provenance, cost and interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The resource freezes its version, explicit mapping, actual source counts,
conditional endpoint frequencies and core values with input hashes and source
membership accounting. Numerical resources and scoring outputs have separate
identities; deriving a resource does not certify predictive accuracy. Recreate
reference inputs with the existing contact exporter and curation recipe, then
run the fit CLI into a new output. Benchmark consumers read CLI-owned energies.

The symmetric interaction table is only part of the estimated reference. Both
residue-conditioned endpoint tables also depend on the contact population. Their
row normalization constrains their entries, and all components share the same
observations; stored cell counts are not independent fitted-parameter counts.
The residue-pair partitions are deterministic functions of these quantities and
must not be counted again as estimated coefficients.

For R reference atom rows and k chemical types (k at most 19), vectorized
mapping/count reductions use O(R) expected grouping work and O(R + k squared)
working storage; no coordinate or dense atom-by-atom matrix is built. After
counts are available, core construction costs O(k squared); partitions for A
binder and B target residue identities cost O(A*B*k squared) directly, with
A and B at most 20. Query scoring maps and joins M physical observations to the
frozen map/core, then groups them by residue edge: O(M) expected hash work and
O(M) working/output storage. Hash-group bounds are expected rather than worst
case. Small dense type arrays and existing vectorized operations suffice; no
new native kernel, dependency, cache or nested worker layer is required.

All increments are dimensionless contact-composition log ratios. Contact endpoint
frequencies are not exposed-atom opportunity. Edge averaging and unit addition
are frozen conventions, with no assay calibration or derivation as a joint
structure likelihood or binding free energy. Original residue scores were
estimated under a different model; this normalization does not prove their
chemical signal has been removed. Reference contacts are correlated and actual
contact pooling can be antibody dominated. There is no radial/angular or site
fit, physical kcal/mol term, accuracy claim or change to the default scorer.

Opportunity-aware class audit
------------------------------

``export --quantity atom-inventory --interface all -o atoms.parquet`` needs no
potential. It retains selected observed heavy atoms across full TRA/TRB, peptide,
MHCa/MHCb and B2M components, including zero-contact atoms. Its conformer policy
and source-heavy-atom indexing match contact observations. Source indices are
not PDB serial numbers. Peptide positions are observed one-based ordinals.
The inventory includes full framework/MHC context and observed Gly extra atoms;
``contact_role``/``contact_eligible`` separately describe the narrower historical
observation scope. No missing atom, hydrogen, water, ion, bond, charge or
protonation state is inferred. Parsing selects one model. Standalone export
writes its complete accounting and residue sidecars, not an energy.

``contact-class-audit --manifest labels.tsv -o audit -t 0 --chunk 128`` rebuilds
label-free structural observations from each raw structure exactly once, with
one batched annotation phase followed by bounded coarse workers. The manifest
requires unique ``id``, binary ``label``, ``receptor_cluster``, ``dataset``,
``cohort``, ``epitope``, ``mhc``, ``mhc_class``, ``template_coverage``, raw-file
``source`` and expected ``sha256``. It creates a new output directory; no resume,
cache, score calculation, fitted potential or default promotion is performed.
Any failed input prevents statistical reduction and leaves complete accounting.
Reduction processes fixed feature blocks while retaining all models and receptor
clusters for each feature. Final Parquet paths may be partitioned datasets; read
them directly with ``polars.read_parquet(path)``.

Contacts join all selected TCR regions, kept separately, to peptide or annotated
MHC helices at the inclusive 5 Angstrom cutoff. Full inventory also preserves
non-target environment atoms. Atom contacts count every physical pair; residue
contacts count each physical residue pair once. Directed features retain TCR
chain/region and amino-acid identities; atom features also retain atom names.
The site view keeps peptide ordinal or MHC author residue/insertion position
within its exact allele/case. The chemistry view pools positions within a case.
Only chemistry features receive an across-case overall summary.

For complex m and feature f, C_mf is its observed contact count. N_B,mf and
N_T,mf count all available matching binder and target atoms or residues,
including uncontacted endpoints. O_mf = N_B,mf*N_T,mf is the number of available
endpoint combinations, and d_mf = C_mf/O_mf when O_mf > 0. These dimensionless
opportunities account for available amino-acid/atom composition; they are not
solvent-accessible, sterically allowed or independent bond opportunities.
No dense universe of all possible feature pairs is constructed. The observed
chemistry-feature union is formed across all cases before its sparse opportunity
joins, retaining eligible zero contacts in cases where a feature was never seen.

Observed binder-only, shared and nonbinder-only sets use class contact presence,
not a claim of chemical requirement or impossibility. Outputs report support,
all-model and eligible-model prevalence, unavailable-opportunity counts, C/O
sums and equal-complex mean densities for both classes. Mean-density differences
and unsmoothed log2 ratios remain explicit, including zero/undefined states.
Each exact epitope/MHC/coverage case has its own contrast. Overall contrasts use
equal weights across eligible cases, retaining coverage and MHC class separately;
TCRvdb epitopes remain separate. Receptor-cluster influence contributions are
summed across cases before variance estimation. CR1 clustered t inference is an
approximation assuming independent supplied receptor clusters; its formulas and
estimability checks are in ``contacts.class_statistics``. Unreplicated classes,
insufficient residual degrees of freedom and degenerate variance do not produce
spurious zero-error significance. No permutation, multiplicity correction or
significance-based feature selection is added.

The output retains source hashes, annotations, raw inventory/physical-pair
partitions, feature dictionaries, per-model C/O and every case/overall statistic.
These are observational diagnostics. Shared gene ancestry beyond the supplied
clusters, naming-equivalent atoms and generator/template effects remain relevant.
Class-exclusive contacts can reflect finite sampling. Discovering contacts from
these labels and rescoring the same models would leak evaluation information;
no filtering benefit or improved final potential follows from this audit alone.

Within-residue atom patterns
----------------------------

``contact-atom-pattern-fit --native native.parquet --external antibody.parquet
--equivalence-resource prior.json -t 0 -o pattern.json`` freezes
``within-residue-atom-pattern-1``. Each physical reference residue edge supplies
unit mass, divided over its physical atom pairs before chemical-equivalence
aggregation. Native class-I CDR:peptide and antibody peptide edges are pooled
without assay or source weights. Unlike integer atom-contact weighting, a dense
edge does not supply more reference mass. Source physical counts and fractional
masses are retained separately.

For directed residues a,b and atom labels i,j, m_abij is the sum of per-edge
atom-pair fractions. P(i,j|a,b) normalizes this mass within the residue pair;
qB(i|a) and qT(j|b) are its directed endpoint marginals. The frozen value is
D_abij = -ln[P(i,j|a,b)/(qB(i|a)qT(j|b))]. These are dimensionless empirical
contact-mark probabilities, not exposed-atom opportunities or independent
observations. All formulae and the evaluation contract are in
``ATOM_PATTERN_SPEC.md``.

``contact-penalties --atom-pattern-resource pattern.json`` adds the arithmetic
mean D of physical atom pairs once per known CDR residue edge to unchanged
original TCRen. Unknown cells abstain with zero and remain in the denominator;
framework and missing-table policies remain unchanged. Separate exports are
required for ``--atom-count-resource`` and ``--atom-pattern-resource`` so their
arms cannot be mislabeled or overwrite each other. ``contact-compare`` includes
``atom_pattern`` alongside original and other requested arms on identical rows.
No assay weights, fitted coefficient, distance regression, joint graph
likelihood, binding free energy or new package default is implied.

Inventory-aware whole patterns
------------------------------

``contact-inventory-pattern-fit`` accepts native and antibody Parquet count
products with directed amino-acid/atom labels, integer ``count`` and integer
``opportunity``. Opportunities include noncontacting cells of the complete
observed atom inventories on existing reference residue contacts. The resource
freezes Jeffreys-smoothed cell probabilities P and a pooled residue-pair density Q.
Both independent-Bernoulli distributions condition on at least one contact.

``contact-penalties --inventory-pattern-resource inventory.json`` adds the whole
binary-pattern negative log P/Q once per scored CDR edge to unchanged original
TCRen. Noncontact terms and the nonempty normalization remain included. Missing
opportunity support abstains for the whole edge. This dimensionless conditional
model assumes independent cells and does not describe binding free energy.
No assay weights, source weights, cohort normalization or default changes apply.
The repository's ``INVENTORY_PATTERN_SPEC.md`` defines the exact estimator,
physical multiplicities, equations and matched-cohort evaluation contract.

Output ``inventory_pattern_pairs`` counts Cartesian opportunity cells, while
``inventory_pattern_contact_pairs`` counts contacting physical pairs.
``inventory_pattern_unsupported_pairs`` and
``inventory_pattern_unsupported_edges`` expose whole-edge abstention. The saved
``inventory_pattern_contact_increment``,
``inventory_pattern_noncontact_increment`` and
``inventory_pattern_normalizer_increment`` sum to the total correction.
