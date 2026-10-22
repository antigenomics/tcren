# Chemical core: frozen numerical specification

Version `chemical-core-count-1`, 20 October 2026. This is one experimental
candidate, specified before its numerical resource or assay evaluation. Keep
original TCRen, the directed atom-count candidate and both surface references.
No assay labels, fitted combination/source weights, radial/angular regression,
new controls, thresholds or default promotion enter this stage.

## Representation and inputs

Use the existing physical 5 Å native class-I CDR→peptide and curated antibody
peptide-target contact populations. Count each physical atom pair once before
chemical pooling. Retain actual integer native/external counts and each source's
hash, PDB/interface membership counts and estimator identity. Pool actual contact
counts with no new source weight; this exposes the antibody-dominated mixture.
Contacts within a structure/family are correlated, not independent examples.

Use the existing element/backbone/apolar/donor/acceptor/cation/anion/aromatic
signature with these versioned changes, without changing the general contact
classifier: split carbon by `atom == C` or `(AA,atom)` in `{D:CG,E:CD,N:CG,Q:CD}`;
separate `M:SD` thioether from `C:SG` unknown-state sulfur; replace `H:ND1,H:NE2`
by one unknown-protonation histidine ring-N type. These rules give 19 observed
types in the accepted reference. Freeze the explicit observed AA:atom→type map.
They do not establish charge, tautomer, hydrogen bonding or a force field. Retain
remaining coarse merges and terminal-state limitations explicitly. Unmapped
query AA:atom labels are unsupported, not assigned a guessed type.

## Symmetric estimator and directed context

Let `n_cd` count directed binder-type c→target-type d contacts, and N their total.
Let `b_c=sum_d n_cd/N` and `t_d=sum_c n_cd/N` be normalized global contact-endpoint
frequencies. For each unordered pair, pool both counts and reference measure:

    m_cd = n_cd+n_dc, Q_cd = b_c*t_d+b_d*t_c  (c<d),
    m_cc = n_cc,      Q_cc = b_c*t_c.

Both m and Q sum to N and 1 respectively over unordered cells. For observed
cells set `W_cd=m_cd/(N*Q_cd)` and `u_cd=-ln W_cd`; use the same values for c,d
and d,c. This is an empirical exchange-invariant interaction assumption, not an
exact fit to the directed contact distribution. Do not average log energies.

For zero-count cells set W=1 and u=0. This is an explicit neutral *core* completion,
not inferred compatibility or a pseudocount. No support floor or smoothing
parameter is selected. The resulting global sum QW need not equal 1; that scalar
normalization is a gauge and is unnecessary for the conditional correction.

For binder AA a and target AA b, estimate separately normalized conditional
endpoint frequencies `qB(c|a)` and `qT(d|b)` from the same directed counts. Define

    Z_ab = sum_c,d qB(c|a)*qT(d|b)*W_cd,
    p(c,d|a,b) = qB(c|a)*qT(d|b)*W_cd/Z_ab,
    Delta(c,d|a,b) = u_cd + ln Z_ab.

The correction is exactly `-ln[p/(qB*qT)]` on positive conditional reference
support. Core symmetry does not exchange qB and qT or the biological roles.
Every available conditional model sums to one. Multiplying every completed W,
including unseen cells, by one positive constant cancels in Delta; changing only
observed cells and resetting unseen W=1 is not a gauge transformation.

Important missingness distinction: an unseen chemical pair with positive qB/qT
gets Delta=ln Z, not zero. A query endpoint absent from the explicit map, a zero
conditional q, or an unavailable AA normalizer abstains: final increment zero,
flagged separately from unseen chemistry. This fallback makes no density-ratio
claim outside reference support. No logarithm of zero or infinite energy is used.

## Application and interpretation

On each physical residue edge with finite original potential support and binder
CDR1/2/3 membership, take the arithmetic mean of Delta over all observed physical
atom pairs. Unsupported pairs contribute zero and stay in this denominator.
Add one such edge mean with unit coefficient to the unchanged original residue
score. Do not collapse atoms sharing a type or replace a physical edge by its
closest atom. Missing all atom observations for an eligible scored CDR edge is
an error. Framework and missing-original-potential policies remain unchanged.
Export `chemical_core` separately from `atom_count`, plus edge increments,
physical-pair/unsupported/unseen counts and resource identity. Both original
residue potentials receive the same chemical increment; neither is refitted.

All values are dimensionless contact-composition log ratios. Unit addition and
edge averaging are frozen conventions, not assay calibration or a derivation of
joint likelihood/binding free energy. Endpoint frequencies are not exposed-atom
opportunity. The unchanged residue score was estimated under a different model;
conditional normalization does not prove its chemical signal is fully removed.
No burial, saturation, orientation or physical kcal/mol term is added here.

## Implementation, checks and evaluation

Reuse the existing physical-contact exporter, batched annotation, coarse worker
pool and original atom-edge aggregation. Freeze reference quantities through a
library CLI; the benchmark only reads CLI-owned energy outputs. A distinct
resource version/CLI flag and score arm preserve all existing definitions.
Vectorized count/group/join operations and small dense type arrays suffice;
no new dependency, native kernel, cache or nested pool is required.

Before evaluation verify hand-counted off-diagonal/diagonal conventions, directed
conditional normalization, whole-table gauge cancellation, unseen-core versus
unsupported fallback, physical multiplicity, original-score conservation, both
potentials' equal increments, input/resource validation, and serial/parallel
output equality. Validate source counts and all 19 mapping states. Bound resource
sizes and keep numerical source/hash validation explicit.

Derive and hash the immutable numerical resource only after this specification
and code checks pass. Record that identity before opening candidate assay results.
Fresh authorized evaluations rebuild from raw structures using all allocated
CPUs through one bounded layer. Preserve all matched cases and all four existing
metrics against original, directed atom counts and both surface scores. Keep CPL,
TCRvdb epitopes, VDJdb coverage/epitope strata and native ATLAS/MD endpoints
separate. Sequence-only legacy substitutions have no mutant atom geometry and
retain an explicit undefined atom arm. Never rerun accepted scores just for prose.
