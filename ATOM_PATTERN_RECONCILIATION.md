# Atom corrections and coherent joint scores

21 October 2026. Scientific design reconciliation; no new potential, structure
export, scoring run, assay fit, coefficient or package default. The two frozen
estimators remain defined by `ATOM_PATTERN_SPEC.md` and the API skill. Complete
reference/score attribution and source hashes are retained in the private audit.

## What edge balancing changes

For a directed atom cell x=(a,i,b,j), a,b are residue identities and i,j atom
labels. Let n_x be its integer physical-contact count and m_x its edge-balanced
mass, each reference residue contact supplying unit mass divided among its atom
pairs. Both references have identical physical observations, equivalences and
observed cell support. Query physical-edge averaging is unchanged. Consequently
this comparison changes the estimated table, not atom ownership, geometry,
unknown-cell support or query aggregation.

Let c^0_x=n_x and c^1_x=g m_x, where g=(total physical atom pairs)/(total reference
residue contacts). This fixed count-unit gauge matches their global totals.
It cancels from the final potential but fixes source bookkeeping. Let C^r_ab be
sum_ij c^r_abij, q^r_B(i|a) the corresponding directed binder atom marginal and
q^r_T(j|b) the target marginal, for r=0,1. Then exactly on observed support:

    D^r_abij = -ln(c^r_abij/C^r_ab) + ln q^r_B(i|a) + ln q^r_T(j|b)
    D^1-D^0 = -ln(c^1/c^0) + ln(C^1_ab/C^0_ab)
              + ln(q^1_B/q^0_B) + ln(q^1_T/q^0_T).

The four terms are cell numerator, residue-pair normalization, binder marginal
and target marginal. All terms are dimensionless. Individual cell/normalizer
pieces depend on the count gauge; their conditional-distribution sum and complete
D change do not. Equal edge mass is a modeling choice, not a correction of a demonstrated bug:
a one-atom edge receives the same reference mass as a dense many-atom edge.
The global native/antibody fraction alone cannot identify
source influence: edge balancing also changes each source's within-pair atom
distribution and the source mixture separately for each residue pair.

For a source-specific attribution, a positive group total t changes from t0 to
t1 along a declared linear interpolation of native and antibody counts. Define
L_inv=ln(t1/t0)/(t1-t0), or 1/t0 if t1=t0. A source contributes
(t1_source-t0_source)*L_inv to the log-total change. Source contributions sum
exactly, including opposing changes whose totals cancel. Apply this identity to
all six log-count terms in D. This is an integrated bookkeeping identity under
a fixed path/gauge, not a causal source ablation or source-only score. No such
alternative score is evaluated.

## Conditioning does not remove residue-dependent mean shifts

For P^r_ab(i,j)=c^r_abij/C^r_ab and Q^r_ab(i,j)=q^r_B(i|a)q^r_T(j|b),

    mu^r_ab = sum_ij P^r_ab(i,j) D^r_abij
            = -KL(P^r_ab || Q^r_ab) <= 0.

KL is Kullback–Leibler divergence, using natural logarithms on the finite atom
alphabet. P is zero outside observed cells; Q is defined by its endpoint
marginals. The expectation uses reference P, not query atom weights. Therefore
D is not generally mean zero conditional on a,b. It can shift residue-composition
preferences when added to a pre-existing residue score, even if a query follows
the reference distribution exactly. This is a property of the likelihood ratio,
not evidence of a coding error or a license to subtract an assay-selected offset.
The implemented neutral abstention on unseen cells is distinct from infinite
log odds when P=0; never claim a normalized full-support likelihood for that rule.

For a known query residue edge with K physical atom pairs and s supported pairs,
write its exact supported contribution as

    mean_query(D) = (s/K) mu_ab
                    + (1/K) sum_supported (D_abij-mu_ab).

Thus a change in a saved edge mean splits into a changed conditional reference
mean and changed query residual. Both are identifiable from stored edge scores,
support counts and the frozen table, even without every query atom identity.
They are algebraic parts, not independently evaluated component scorers.
Separate native/antibody contributions to mu are weighted expectations under
the corresponding pooled source probabilities, not isolated-source energies.

The full numerator/marginal source attribution requires the complete physical
atom-label histogram. Closest atoms and an edge average cannot recover it.
Where only a residue-edge ledger survives, retain the identifiable aggregate
pieces and explicitly leave the remaining source terms unidentified. This is an
evidence limitation, not a reason to stop other scientific work or rerun structures.

## Correction versus recalculating the whole score

These are not intrinsically different mathematical models. For one common
reference measure S, define the directed joint atom-mark probability
p_S(a,i,b,j), its residue distribution p_S(a,b), and corresponding binder and
target marginals p_SB and p_ST. On common positive support:

    E_res^S(a,b) = -ln[p_S(a,b)/(p_SB(a)p_ST(b))]
    D_S(a,i,b,j) = -ln[p_S(i,j|a,b)/(p_SB(i|a)p_ST(j|b))]
    E_joint^S(a,i,b,j)
        = -ln[p_S(a,i,b,j)/(p_SB(a,i)p_ST(b,j))]
        = E_res^S(a,b) + D_S(a,i,b,j).

This is an exact chain-rule identity, with no fitted coefficient. At the query
level an arithmetic atom-mark mean carries E_res once per fully supported residue
edge. Unknown-cell abstentions require their existing s/K accounting separately.

If E0 is an older residue table trained under another reference or observation
measure, the implemented hybrid is instead

    E0 + mean(D_S) = mean(E_joint^S) + [E0-E_res^S]

on fully supported edges. The bracket is the retained residue-prior difference.
It can be useful, but no exact single-reference likelihood interpretation follows
by calling both terms dimensionless. Dropping E0 and keeping only D_S would
discard residue-pair enrichment; it would not be complete atom-level rescoring.
Directly summing atom pairs would also introduce linear contact multiplicity,
whereas the present mean gives one evidence unit per residue edge.

A deliberate hybrid can nevertheless be probabilistically coherent. If the
original residue table has documented normalized P0_res and Q0_res with
E0=-ln(P0_res/Q0_res), define

    P_H(a,b,i,j) = P0_res(a,b) P_S(i,j|a,b),
    Q_H(a,b,i,j) = Q0_res(a,b) q_SB(i|a)q_ST(j|b).

Both are normalized distributions, and -ln(P_H/Q_H)=E0+D_S on common positive
support. Thus different reference sources do not by themselves make an additive
correction invalid or double-counted. The assumption is a deliberately chosen
hybrid conditional model; its usefulness for recognition still needs evidence.
For the actual empirical tables, their original normalization/support policies
must be checked before claiming this exact interpretation. Arithmetic averaging
over dependent atom marks remains a composite statistic, not a joint whole-edge
likelihood. This distinction favors preserving the useful residue prior rather
than requiring a common-source rebuild merely for aesthetic consistency.

**Recommended working route:** preserve E0 and formulate one explicit conditional
atom-pattern likelihood ratio, with the same context/support definition in its
numerator and denominator. In parallel, specify the coherent joint model using
that identical reference so its residue component is visible. Compare these
representations only after each is fully specified and frozen, never by choosing
an assay-fitted blending weight. Existing directed, edge-balanced and symmetric
chemical-core results remain separate experimental arms.

## A more informative energy specification

A complete contact pattern contains more than a sampled contacting atom mark.
For residue edge e, let A_e be its observed atom availability, K_e the count of
physical atom pairs, z_e their pattern, and h_e a declared physical context such
as receptor region and peptide site. Let P_N describe a compatible native
reference and Q a clearly specified opportunity/model measure on identical
support. A design, not an implemented or evaluated energy, is

    E = sum_e E_res(a_e,b_e|h_e)
        - sum_e ln[P_N(K_e,z_e | a_e,b_e,A_e,h_e)
                   / Q(K_e,z_e | a_e,b_e,A_e,h_e)]
        + Phi(interface configuration).

Every term is dimensionless. E_res is the chosen residue reference log ratio;
Phi is a separately defined dependency/environment correction, not a free fitted
coefficient. A sum over edges is a composite approximation when their atom
patterns share atoms or interact. An exact joint likelihood would require their
conditional dependencies. Splitting the conditional term into a count term and
a pattern-given-count term is another exact probability factorization, provided
both P and Q are defined. Do not turn it into independent Bernoulli atom pairs
without explicitly justifying that approximation.

One low-dimensional chemical construction is P_theta(z|a,b,A,h) proportional
to Q0(z|a,b,A,h) exp[-U_theta(z,a,b,A,h)], giving correction
U_theta+ln Z_theta. Q0 is the declared available-pattern measure, theta is the
reference-estimated chemistry and Z is the conditional normalizer. Symmetric
chemical interactions can coexist with directed residue/environment marginals.
The earlier symmetric chemical-core model already used an interaction-plus-
normalizer form and had substantial losses. Repeating that formula with the
same contact-endpoint reference is not a new scientific solution. A distinct
candidate needs actual new conditional information: feasible atom patterns,
geometry, environment or competition, with sufficient observed reference support.

Antibody data can support shared chemical interactions. They do not establish
TCR-specific paratope, peptide-site or MHC-pocket opportunity distributions.
Full atom inventories permit available versus unavailable endpoints; 5 Å pair
lists alone establish neither bond satisfaction nor buried unsatisfied polarity.
A physical Hamiltonian would also require consistent units, solvent/state and
protonation definitions. Adding kcal/mol terms numerically to contact log odds
is not justified by a label of 'energy'. No such term is fitted here.

## Why a coarse statistical potential can remain competitive

The following are explanations to discriminate, not results of new controls:

- A residue table aggregates many dependent atomic observations into fewer
  quantities and can tolerate uncertain predicted side-chain placement. More
  atom cells do not create more independent native complexes. An atom model can
  therefore add estimator and coordinate noise to a useful coarse signal.
- A native-versus-independent-endpoint log ratio is not the desired binder-versus-
  nonbinder likelihood ratio for generated complexes. Both labeled classes are
  generated by the same structural pipeline; common packing priors can make
  many apparently reasonable atom contacts nondiscriminating for recognition.
- Reference mixing changes within-pair preferences and directed marginals, not
  only a global antibody weight. An abundant antibody atom pattern may describe
  its own paratope context better than a TCR–peptide context.
- Conditional corrections retain residue-dependent reference KL offsets, while
  unit edge averaging discards linear atom multiplicity. Either may alter useful
  original rankings without measuring new binding information.
- The poor symmetric-core outcomes show that reducing cell count alone has not
  solved the reference/context problem. They do not prove atomic information is
  absent; the completed within-contact audit already observes such information.

TCRmodel2 uses modified AlphaFold with receptor and pMHC templates; that is a
structural modeling procedure, not a measurement of whether an arbitrary input
pair binds. Extending this fact to predicted binders and mismatches above is an
explicit inference about the benchmark task, not a newly tested claim.
[Primary paper](https://doi.org/10.1093/nar/gkad356),
[PubMed37140040](https://pubmed.ncbi.nlm.nih.gov/37140040/).
The reference-ratio viewpoint explains why a score's denominator and generating
measure determine its interpretation; it does not certify the new candidate.
[Verified arXiv1008.4006](https://arxiv.org/abs/1008.4006),
[DOI10.1371/journal.pone.0013714](https://doi.org/10.1371/journal.pone.0013714).

## Next scientific steps

1. Complete the frozen-output source/marginal/conditional-mean attribution and
   preserve every gain/loss and missing-identifiability boundary.
2. Keep E0 and all tested arms. Specify the conditional pattern state and reference
   denominator before considering more cells, distance regressions or weights.
3. Use the existing full atom ledgers where available to identify which pattern
   information is genuinely observed. Other candidates' missing atom ledgers are
   documented data issues, not permission to rebuild them during this audit.
4. Decide whether a coherent joint reference can preserve the desired residue
   preference, or whether the deliberate hybrid prior E0-E_res^S should remain.
   No replacement is promoted by algebra alone.
5. Any later estimator/data extension needs an explicit freeze before evaluation,
   exact cohort matching against original and both surfaces, and the existing
   all-task/legacy coverage record. Do not mask losses by reweighting assays,
   changing cohorts or publishing a new interpretation of an old score.
