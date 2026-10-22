# Inventory-aware whole contact patterns

## Frozen candidate, original residue score retained

This separate experimental arm adds a dimensionless conditional pattern log ratio
once per original scored CDR:peptide residue edge. Original TCRen and both bundled
potential tables remain unchanged. No assay labels, assay weights, source weights,
cohort normalization, distance regression or default promotion enter estimation.
The reference resource and immutable scoring package must be frozen before any
new benchmark score is inspected.

For a contacting directed residue pair with amino acids a,b, let A and B be its
selected observed heavy-atom inventories. Their Cartesian product contains M
physical atom-pair opportunities. Missing atoms are never built. For each physical
pair t=(i,j), z_t=1 when its distance is at most the existing 5 A contact cutoff,
and z_t=0 otherwise. Chemical equivalence pools parameter categories, retaining
every physical pair and its opportunity multiplicity. The observed state is the
whole binary vector z, conditional on A,B,a,b and K=sum(z)>0.

Native class-I CDR:peptide and curated antibody:peptide references supply integer
contact counts c_abij and available-pair counts n_abij across their contacting
residue edges. Actual counts are pooled, with separate source accounting. The
fixed Jeffreys Beta(1/2,1/2) posterior means are

    p_abij = (c_abij + 1/2) / (n_abij + 1)
    q_ab   = (sum_ij c_abij + 1/2) / (sum_ij n_abij + 1).

These are reference estimates, not coefficient-free quantities. P assumes
independent Bernoulli states with probabilities p_t. Q assumes exchangeable
Bernoulli states with the common q_ab. Both distributions are normalized on the
same physical inventory and conditioned on at least one contact. Define

    P_+(z) = product_t p_t^z_t (1-p_t)^(1-z_t) / (1-product_t(1-p_t))
    Q_+(z) = product_t q_ab^z_t (1-q_ab)^(1-z_t) / (1-(1-q_ab)^M).

The edge correction and retained-residue score are

    D_edge = -log(P_+(z)/Q_+(z))
           = -sum_t [z_t log(p_t/q_ab)
                     +(1-z_t) log((1-p_t)/(1-q_ab))]
             +log[(1-product_t(1-p_t))/(1-(1-q_ab)^M)]
    E_inventory = E_original + sum_scored_CDR_edges D_edge.

Every logarithm is natural. The correction has coefficient one; it is summed
once per whole residue pattern, without atom-pair averaging. No pattern is
silently restricted to contacting endpoints. Noncontacting atom pairs and the
nonempty normalization are essential terms. Counts, inventory size, contacting
pairs and whole-edge abstentions are reported separately.

An unseen amino-acid/category opportunity abstains for the entire edge with a
zero correction and an explicit unsupported-edge count. This is a support
policy, not a likelihood for an invented subset. Original framework and missing
matrix-entry policies remain in effect. Invalid or inconsistent inventories,
duplicate physical pairs and contacts outside their inventories are errors.

## Meaning and limits

This normalized conditional model is a deliberately sparse first approximation:
independent atom-pair Bernoulli variables cannot reproduce geometrical exclusion,
bonds or correlated packing. Empirical frequencies come from already contacting
residue edges; the model is not maximum-likelihood fitting of a truncated
Bernoulli family. Q describes exchangeable cells within these edges, not bulk
solvent, free atom exposure or a binding/nonbinding ensemble. Full observed
inventories do not establish protonation, hydrogen bonds or solvent exposure.
No physical free-energy units or assay-specific interpretation are assigned.

Antibody references augment actual available-pair chemistry; their greater
sample count remains visible. Transfer does not establish exchangeability of
antibody and TCR recognition contexts. No source-only classifier or source-weight
search is implied by the pooled resource.
