# Contact plausibility and positive penalties — 18 October 2026

## Current author subgoal — 20 October 2026

The ordered continuation is
CONTACT_DISCRIMINATION_SUBROADMAP.md (private evaluation record; not distributed).
It reconciles generator and negative-label construction, native comparisons,
antibody atom-reference support and the completed accessible-boundary losses.
The immediate table-only dissection is complete; no structure rerun or new fit
was performed. Existing surface arms still use residue-only chemistry and do
not test an antibody-derived atom-specific matrix. Preserve the original score;
define the atomic observation unit, reference opportunity and transfer context
before estimating any next table. Classification restoration remains a measured
scientific objective, not a consequence of geometric conservation alone.

## First frozen comparison completed — 18 October 2026

The author approved worst-supported atom/marginal aggregation and transformed
Gaussian marginals with Student-t predictive tails, actual external prior counts
and declared pooled chemistry backoff. The resource was frozen before scoring.
[CONTACT_PENALTIES.md](audit/2026-10-18/CONTACT_PENALTIES.md) records the complete
3808-model Stratos comparison, exact original-energy reproduction, bounded memory,
all four arms and cohort-specific gains/losses. The full report is in the visible
benchmark checkout, `results/contact-penalties-20261018/REPORT.md`. This completes
the first approved resource/comparison; full chemical angular rules, line-of-sight,
broader PPI/atlas priors and independent validation remain open. The source
reconnaissance below records what was available before raw correspondence recovery.


Status: active goal; author approved the first implementation stage. Exact CLI
atom observations and the fresh native summary are complete and independently
checked; see `audit/2026-10-18/CONTACT_OBSERVATIONS.md`. Quantitative source
research is recorded separately. No learned table, penalty scorer or accuracy
improvement is claimed. Source request: manuscript `plan181026.txt`,
including the author's positive-penalty extension on 18 October.

## Goal and intended interpretation

Develop directed binder-to-target, literature-informed rules for contact plausibility and
test whether a frozen contact penalty improves energy-based ranking relative to
the original energy and matched ipTM. Assess artificial-contact contributions as
a mechanism, without assuming they explain every energy error. The primary new
candidate family adds positive energy for unnatural contacts; deletion/down-weighting
alone need not penalize a forced interface.

The native reference distributions and Bayesian posterior parameters are estimated
dependencies, frozen before scoring. This follows the permitted fixed empirical
potential route; it is not coefficient-free or training-free. No assay-fitted
classifier weights, cohort normalization, confidence mixing or retired fitted
readout enters final scoring. Existing scores and publication findings remain.

## What the existing evidence and implementation establish

- `contacts.geometry.all_atom_contacts(..., atom_pairs=True)` already enumerates
  heavy-atom pairs within the inclusive 5 Å cutoff using cKDTree. Its default is
  one closest atom pair per residue pair. Reuse it rather than write another
  neighbour search. `ContactMap.from_structure(..., atom_pairs=True)` and
  `interface(..., tcr_regions="cdr")` already provide annotated interface roles.
- `contacts.table.tidy_contacts` mirrors pairs internally; a physical observation
  must be counted once after selecting/canonicalizing the required interface.
  Residue counts and potentials retain binder-to-target orientation; no reverse
  counts or transpose averaging enter their estimation.
- Existing `contact_types` v2 exposes overlapping atom-type/distance heuristics,
  including residue-level aggregation over all atom pairs and ring stacking.
  These labels are starting annotations, not full angular or energetic validation.
  Reuse them and audit any chemical rules before extending them. The existing
  exact-export chemistry gap is recorded as library issue 59.
- `export --quantity contact-context` and `contact-summary` own the physical
  residue-edge context, original potential entries and canonical reductions.
  Preserve their identities, accounting and zero-contact models.
- Benchmark `CONTACT_DIFFERENCE.md` records 617 models and five matching native
  complexes. GLC's inversion is mainly in shared motifs; a universal deletion of
  unusual or polar contacts does not follow. Shared motifs can have different
  geometry, and native-compatible contacts can have positive potential entries.
  Motif equality is not atom-coordinate equality or physical-contact truth.
- Library `SOURCES` records 374 deposited Native2026 structures and 359 eligible
  structures in the current potential recipe. The benchmark's older source
  narrative still says 362; it is not the current derivation specification.
  The read-only local mirror inventory contains 374 structure files plus one
  metadata file in each of Native2026 and Canonical2026. Mirror HEAD observed:
  `6a56cbe2de5feb18f243ad121b72c3b03aff9ce5`, locally ahead of its tracking ref.
  No fetch, synchronization or upload was performed. Native and canonical files
  represent the same complexes, not independent observations.

## Scientific design

### Observations, context and direction

Start with the existing native inventory; derive the actual eligible set with
the current alpha/beta, standard-peptide and MHC-I/II rules. Report integer
structure, residue-edge and heavy-atom-pair counts separately; the source plan's
approximately 300 structures / 15000 contacts is a planning estimate. Do not
interpret approximately 50 residue edges as approximately 50 heavy-atom pairs.

Export each selected heavy-atom pair with structure/assembly identity, chain
role, author residue/insertion code, atom identity/altloc, residue type, distance
in Å, region and missing-geometry flags. Enforce a documented conformer policy
before estimating distributions; preserve the old contact definition separately.
Glycine contributes observed backbone atoms, not fabricated Cβ atoms. Missing
non-glycine side chains are missing observations, not glycine substitutions.

Individual atom-pair chemical identities may be exchangeable: canonicalize
whole amino-acid:atom tuples, never residue and atom names independently.
Residue observations, potentials and angular distributions retain binder→target
order. Theta/phi also retain the endpoint whose angle is measured. Biological
roles, context and physical endpoint identities must survive every reduction.
Keep TRA/TRB, CDR1/2/3, MHC class and groove-helix membership available as context;
do not declare every MHC-chain contact to be a helix contact. Internal CDR contacts
are a separate structural category, not interface energy added by accident.

### Priors and sparse-data updating

For each declared categorical observation space, use a Dirichlet prior with
alpha_k = tau * q_k, where q_k > 0 is a normalized verified external prior,
tau > 0 is its effective observation count, and k is a directed residue category or an explicitly exchangeable individual atom category.
For integer counts n_k from the selected native structures, the posterior
predictive category mass is (n_k + alpha_k)/(N + tau), N = sum_k n_k.
Diagonal categories count once. Raw contact counts are correlated within
structures; they are not N independent complexes. If existing redundancy weights
are used, expose their units and scale: fractional updates are a weighted
likelihood approximation, not an exact integer multinomial experiment.

Literature descriptions do not supply numerical q_k, tau, angular variances or
complete atom-pair tables by themselves. Retrieve exact published tables or
document a reproducible licensed reference-corpus recipe before implementing
an informative prior. Separate general packing, PPI and antibody evidence;
Tyr enrichment alone is not an atom-geometry prior. No invented prior strength,
empirical-Bayes calibration on test labels, or unrequested reference-corpus run.

Use a few marginal or low-order geometric distributions. Periodic sin/cos
embeddings preserve wraparound but are not automatically Gaussian: a two-coordinate
embedding lies on a circle. Choose a support-correct family or a justified
normalizing transform from the external/reference data, preserving its Jacobian
when comparing densities in different coordinates. Keep boundary angles,
degenerate vectors, equivalent atom names and missing hydrogens explicit.
Directed theta/phi descriptors must exchange with their endpoints. No full
high-dimensional amino-acid-specific angle covariance or GP is initially needed.
Posterior predictive checks requested by the source plan assess reference
distributions; they do not certify classifier performance or contact independence.

### Positive penalty candidate

Let R be the unique physical residue-pair contact set, E0 the original fixed
potential sum, and U_r >= 0 the plausibility penalty for residue pair r in R.
The primary candidate is Enew = E0 + sum_(r in R) U_r. Lower remains favourable.
Retain E0, each U_r and their total separately, with Å/radian inputs and explicit
dimensionless statistical-potential output. This is not binding free energy.

A possible U_r = -ln(s_r) requires a defined compatibility score 0 < s_r <= 1.
Posterior density may exceed one and has units; -ln(density) is not automatically
nonnegative or dimensionless. Native-only observations estimate native-contact
distributions, not P(natural | contact) without an alternative model/base rate.
Use a specified predictive compatibility construction and label it accordingly;
do not introduce an artificial-contact null corpus without author authorization.

Define atom-to-residue aggregation and the penalty scale before evaluation.
Adding a separate penalty per atom pair would favour or penalize large side
chains merely through multiplicity; it must not silently replace residue-level
counting. Likewise, combining marginal tail scores is not an exact joint model.
Choose one documented aggregation after the reference summary, then freeze it
across the finite penalty menu below; do not scan aggregations against labels.
The unit coefficient above fixes a
convention, not a physical calibration. Any different scale requires a reference
derivation and explicit provenance, not accuracy optimization.

Distinguish chemically impossible observations, severe clashes, unresolved
geometry and merely sparse observations. Rarity alone is not a hard rejection.
Line-of-sight needs a stated blocker/radius/covalent-neighbour definition; do not
invent a universal criterion. Preserve explicit missing/error states instead of
turning unavailable angles into zero penalty. Zero predictive support requires
a declared rejection or justified prior support, not an arbitrary numerical floor.

### Author-authorized finite penalty menu

The author subsequently requested brainstorming and comparison of fixed,
matrix-scaled and sign-flipped penalties. This explicitly expands the earlier
single-penalty candidate scope; it does not authorize a coefficient/threshold
sweep. Keep the same frozen unnaturalness assignment for every arm.

Let e_r be the original table entry for physical residue contact r, u_r in [0,1]
its frozen unnaturalness weight (binary flags are permitted), and
M = max_(a,b) |V_ab| over finite entries of the full selected frozen table V.
Define M independently of the scored cohort; keep absent entries explicit.
All expressions below give the resulting per-contact contribution e'_r.

| Candidate | e'_r | Increment relative to e_r |
|---|---|---|
| Fixed additive, one natural-log unit | e_r + u_r | u_r >= 0 |
| Matrix-scaled additive | e_r + u_r M | u_r M >= 0 |
| Literal sign flip, diagnostic | (1 - 2u_r)e_r | -2u_r e_r, potentially negative |
| Flip favourable entries only | e_r + 2u_r max(-e_r, 0) | 2u_r max(-e_r, 0) >= 0 |

Proposed fifth arm, for discussion: replace a fully flagged contribution with
+M, using e'_r = (1-u_r)e_r + u_r M. Its increment u_r(M-e_r) is nonnegative
for entries in V. This is a common positive replacement value, unlike adding M.
It is not yet part of the declared comparison. If M = 0, the matrix-scaled arm
is degenerate; report this without silently substituting another scale.

At u_r = 1 the literal flip produces -e_r and the favourable-only flip produces
|e_r|. A raw flip rewards unnatural contacts whose original entry was positive,
so it is a signed diagnostic, explicitly outside the nonnegative-penalty family.
All arms leave contacts unchanged at u_r = 0. Export e_r, u_r, the increment and
e'_r per physical residue pair, then sum once per model. No new nuisance weight,
cohort rescaling or label-chosen direction is needed. Keep the baseline and every
arm's matched results. Selecting among arms on these already inspected labels
is development evidence; independent predictive validation remains separate.

## Work sequence and verification

1. **Reference summary and exact CLI export.** Audit the eligible raw inventory,
   annotation batching, physical identities and interaction scopes. Extend the
   existing export with the minimum required atom/geometry fields. Verify reverse
   invariance, diagonal counts, glycine, altlocs, cutoff equality, missing atoms,
   empty/all-failed cases and unchanged legacy exports with small deterministic
   checks. Summarize observed atoms, distances, geometry and context without labels.
2. **External-prior dossier.** Retrieve the quantitative sources, support and
   strengths; document availability, license, overlap and reproducible recipes in
   SOURCES. Compare native shifts with the prior on matching definitions. If a
   required numerical prior is unavailable, report that exact gap rather than
   invent one. Keep every estimated dependency inspectable.
3. **Frozen plausibility resource and penalty.** Specify the small supported
   posterior families, predictive compatibility, residue aggregation, finite penalty menu and
   hard/missing cases using steps 1–2. Obtain review of those concrete choices,
   derive a versioned resource and expose application through the library CLI.
   Verify permutation/rigid-transform invariance, nonnegative increments for the
   penalty arms and the signed diagnostic's documented behaviour,
   conservation of Enew = E0 + total penalty, and immutable resource loading.
   Update API skill, docs and source provenance together. No scoring-time fitting.
4. **Matched ranking evaluation.** Pin an immutable noneditable wheel and resource
   digest in the benchmark. Rebuild all newly requested derived outputs from raw
   inputs into a new run directory. Compare every frozen penalty arm, original energy
   and native ipTM on exactly matched rows, with existing ROC/AP and established
   operating-point metrics. Retain labels, missingness, zero-contact models,
   separate CPL panels, separate YLQ/GLC, VDJdb cohorts/template strata and MHC
   classes. Report observation counts, all endpoint deltas and every regression.
   Do not extrapolate an untested class/stratum or choose directions/thresholds
   from the results. No additional controls, sweeps or significance tests.
5. **Interpretation and integration.** Inspect exact contact decompositions and
   selected geometry to distinguish changed contact plausibility from wrong pair
   potential signs or wrong poses. Update findings only with supported evidence;
   disclose any weakened standing claim before publication prose changes. Commit
   scoped library changes and deliver into the visible checkout; integrate
   benchmark/manuscript records through their owning repositories. Keep MD and
   release approval separate. Performance improvement is a measured outcome,
   never a promised completion criterion.

Use existing vectorized tables and one batched native search first. Profile the
new geometry stage before considering native extensions or process parallelism.
The authorized native summary is recorded in the stage-1 evidence above. Benchmark writes/compute
require the owning project workspace access; current writable roots cover the
library and manuscript. No computed-intermediate reuse or resume paths.

## Verified literature routes and remaining limits

Metadata checked through PubMed on 18 October 2026; retain identifiers rather
than reconstruct bibliography fields. These routes establish relevant methods,
not ready-to-use numerical Bayesian priors:

- Arpeggio: [PMID 27964945](https://pubmed.ncbi.nlm.nih.gov/27964945/),
  [DOI 10.1016/j.jmb.2016.12.004](https://doi.org/10.1016/j.jmb.2016.12.004).
  Full text confirms atom typing plus distance/angle criteria, and distinguishes
  proximity and angle-free polar labels from stronger interaction assignments.
- Protein Contacts Atlas: [PMID 29335563](https://pubmed.ncbi.nlm.nih.gov/29335563/),
  [DOI 10.1038/s41594-017-0019-z](https://doi.org/10.1038/s41594-017-0019-z).
  Metadata/abstract verified; numerical tables and full-text retrieval remain open.
- trRosetta orientations: [PMID 31896580](https://pubmed.ncbi.nlm.nih.gov/31896580/),
  [DOI 10.1073/pnas.1914677117](https://doi.org/10.1073/pnas.1914677117).
  Full-text geometry definitions distinguish atom-distance exchangeability from directed
  theta/phi. Its learned neural outputs are not adopted as a literature prior.

The [quantitative source dossier](audit/2026-10-18/CONTACT_PRIOR_SOURCES.md)
identifies the large antibody survey and its deposited contact table/producer,
plus PPI atomic-property and count supplements and pinned Arpeggio criteria.
The antibody deposit is now downloaded and checksum-verified: 856751 rows from
1833 PDB IDs, including 529 repeated physical atom-pair rows and two Native2026
PDB overlaps. The dossier records exact hashes, counts and required deduplication/
overlap handling. Its reference table contains antibody atoms only; complete
antigen geometry still needs raw structures. PPI supplement retrieval and
quantitative prior construction remain open. The atlas
book's numerical distributions still require retrieval. No numerical Bayesian
prior or penalty has been extracted from these sources.
Assay negatives remain functional labels rather than solved nonbinding poses;
native/test overlap and structural redundancy must be disclosed. None of these
limits removes the observed GLC inversion or authorizes changing existing results.
