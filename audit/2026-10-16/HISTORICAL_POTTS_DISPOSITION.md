# Historical Potts decomposition and CLI disposition — issue 32

Audited from library 69759c8, benchmark 8945dd7 and manuscript 86a617a.
This closes the historical disposition, not the unsupported feature implementation.
No fitted replacement, model refit or historical CPL cohort rerun was performed.

## Definitions and exact overlap

For partner position i, let S_i be its available receptor:peptide residue pairs;
a and b range over the twenty standard residues. Define C_i[f](a) as f_i(a)
minus its twenty-residue mean. For site s, eta_s(a)=r_s+h_par(a)+J_s(a), where
r_s is independent of the candidate residue, h_par is the fitted partner field
and J_s(a) is the fitted residue-pair term. Let p_s be the native-sequence
coupled marginal, held fixed, and sigma_s the observed contact indicator.

| Original arm | Centred quantity | CLI status |
|---|---|---|
| F_full | C_i[sum softplus(eta_s(a))] | scan dF, identical sites/model |
| F_nocomp | C_i[sum softplus(r_s+J_s(a))] | unsupported complete route |
| J_lin | C_i[sum J_s(a)] | unsupported complete route |
| F_full_p | C_i[sum p_s eta_s(a)] | scan --coupled dF, identical native marginals |
| J_lin_p | C_i[sum p_s J_s(a)] | unsupported complete route |
| J_lin_c | C_i[sum sigma_s J_s(a)] | unsupported complete route |
| F_full_c | C_i[sum sigma_s softplus(eta_s(a))] | unsupported complete route |
| J_lin_pc | C_i[sum sigma_s p_s J_s(a)] | unsupported complete route |

F_nocomp_p equals J_lin_p after centering: their difference before centering is
sum p_s r_s, a candidate-independent constant. It is not a ninth exported arm.
All original arm names, column names and numerical calculations are unchanged.

Positive log Z0=sum softplus(eta) is the exact independent-site log-partition,
in dimensionless natural-log units (nats). The statistical free energy of that
model is -log Z0, with the opposite sign; neither is a calibrated binding free
energy. The coupled derivative d(log Z)/d eta_s=p_s gives first-order response
sum p_s delta eta_s. Finite substitutions do not make it exact. In the existing
coupled scan, the column still named log_z0 contains the frozen-marginal linear
surrogate, including other positions' constant contribution; it is not coupled
log Z. Only centred dF overlaps the original corresponding decomposition arm.
CLI/API help now makes these distinctions without renaming or changing outputs.

The indicator and p-weighted predecessors use observed contact support rather
than all available sites. Missing marginals retain a unit weight, as originally
implemented. For n_i observed contacts at a position, the normalized weights are
n_i*p_ij/sum_j p_ij, with zero weights when the denominator is zero. The pseudo-
contact version is n_i*(p_ij+1/n_i)/(sum_j p_ij+1); zero or uniform marginals give
unit weights. The normalization fits no coefficients, but its input model is
fitted. Position-dependent scaling can change whole-peptide ranking; the old
unqualified invariance wording was incorrect. Within-clone means/scales used by
composed arms are estimated dependencies. These historical quantities do not
become fit-free because labels are absent at scoring time.

## Geometry, identity and examples

The benchmark already forwards model radius/cutoff into availability and codes
the marginal frame's own site order. Its physical-key join keeps distinct receptor
and partner chains and checks residue letters. Existing synthetic checks exercise
nondefault geometry and reversed input order without reannotating per item.

Fresh whole-batch CLI examples use public 1OGA (class I) and 1FYT (class II):

| Geometry | Class-I available sites | Class-II available sites | Scan rows per class |
|---|---:|---:|---|
| bundled 15 A / 5 A | 176 | 191 | I: 180; II: 260 |
| custom 12 A / 4.2 A | 82 | 83 | I: 180; II: 260 |

Default, custom and fixed-native-marginal scans each match the existing API
exactly on all 440 rows, for dF, log_z0, n_pairs and is_observed. Every available
pair maps to the exact exported chain, region, author residue/insertion identity
and amino acid. Both class-II groove chains remain in annotation. Classes are
reported separately; this is software evidence, not class-II ranking validation.
Tables and identity joins were manually inspected, including peptide positions
versus author numbering. No molecular image or final-paper figure changed.

The six CLI calls completed in 39.05362029199023 s, peak sampled concurrent tree
RSS 1,223,475,200 bytes (1.14 GiB), 13 live OS threads including bounded auxiliary
I/O. The selected budget was one CPU; actual MMseqs processes reached one thread
and --threads 1. Existing whole-set annotation APIs were used, with no outer
search pool or per-structure search calls. Sampling interval was 20 ms with eight
transient process-read misses; sampled peaks are not continuous proofs. The
300-second/8-GiB bounds were not reached. API equality checks are separate from
that six-command timing. Existing batch APIs and columnar joins sufficed.

The exact installed wheel is the previous immutable 4.0.1 artifact, SHA256
`d24f599861773b9d3ff69c174c6059e14aa245b46e64878e57954611197dceba`;
all 182 payload files match. Executable ASTs of scan and peptide_free_energy match
between that wheel and current source; only documentation changed. Eight focused
library scan/identity tests and two benchmark geometry/physical-copy tests pass.
Ruff and Sphinx -W pass. No release/install of a new version is claimed.

## Preservation and disposition

The private benchmark retains all 18 original decomposition/normalization tables
byte-for-byte, with their hashes, and all 430 recorded comparison rows: 398 finite
and 32 originally unavailable fields. Original matched baselines and every loss,
gain and tie remain. No population, score direction, threshold or baseline was
changed. The overloaded normalization-summary r column is explicitly identified
as Pearson R at cell level versus macro mean ROC-AUC at library level; its
r_lo/r_hi are per-clone extrema, not confidence intervals.

The owning quantity ledger records every exact/approximate equation and CLI gap.
Complete composition/support/normalized-CPL routes remain unsupported historical
producers. No final-paper dependency requires them; adding a fitted scoring arm
merely to close their feature gap would contradict the final specification.
Issue 32 is retired by that explicit disposition, not by a claim that scan is a
complete replacement. Reopening historical CLI support would require preserving
these exact quantities and separate provenance, not substituting another score.

The final manuscript's 112-file scientific graph, 40 selected historical macros
and all 15 molecular assets still pass their read-only gates. Fitted Potts remains
outside the final paper. All pre-existing dirty files are preserved; the long MD
campaign in issue 26, release approval and deterministic topology reductions
remain independent work.

## Limits retained once

- The diagnostic coupled example uses four chains, burn 16, 32 draws, thin 1,
  seed 0 and one worker. It verifies implementation equality, not convergence.
- The original rounded historical outputs retain incomplete original execution
  lineage. Fresh examples do not establish a full historical cohort reproduction.
- A radius-19 / cutoff-4.2 diagnostic was rejected by the existing fitted-block
  domain guard. Changing radius alone cannot supply missing categorical levels;
  the failed attempt remains separate. No silent clamp or substituted model.
- The first R summary attempt exposed factor-valued melt columns. Explicit string
  columns fixed the summary; the accepted output was written into a new directory.
