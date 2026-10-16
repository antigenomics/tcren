# Final-model and producer disposition — issues 33 and 35

Audited on 16 October 2026 from library master 98bb84d, benchmark ea784be and
manuscript c67a096. The owning private manuscript now records the exact scientific
include/drawing graph and primary producers: 112 project files and 126 edges,
including all 15 PNGs and 38 unique PDF inputs (37 outer assets plus a nested
vector plot). This excludes TeX runtime packages. Every active scientific
production family has its recorded CLI/R inputs and formatter named.

## Final formulas and estimated dependencies

Final readouts are individual interface sums Phi_I = sum_{(i,j) in C_I} V_I(a_i,a_j),
structural differences Delta_c Phi = sum_{I containing c}[Phi_I(x)-Phi_I(x_Ala)]
and individual structural descriptors. C_I is the observed heavy-atom contact set
at 5 Angstrom; V_I is the selected fixed empirical potential; x_Ala is the
coordinate-truncated reference with contacts rebuilt. Sums are dimensionless
contact preferences, not calibrated free energies. No universal fitted composition
is selected. Empirical counts, redundancy balancing, pseudocounts, potential
identities and fixed frame/chemistry conventions remain explicit dependencies.

Only 40 used definitions from 477 historical macros enter the paper, copied
verbatim: 17 fixed-template CPL values, three descriptive kinetics values and
20 potential values. Unused Gaussian/Q/T/GLM/Potts constants no longer load
through those mixed bundles. Complete original bundles remain unchanged.
Retired fitted DATs no longer enter active Make prerequisites via a wildcard.

The retained fixed-template reconstruction is mean_b Phi(x_i->b)-Phi(x_i->a),
where b ranges over the twenty standard residues on fixed partners. Its alanine
cell minus its native-residue cell gives Phi(native)-Phi(alanine) on that fixed
map. The former text conflated the alanine cell with a per-site scan; the
identity and distinction from rebuilt structural references are now explicit.
No reconstruction cell, observation or empirical correlation was changed.

The retained kinetic trend estimates an OLS slope/intercept and residual
variance on 23 class-I observations, with its recorded 95% mean-response band.
Those coefficients never enter a ranking arm or predicted rate. CPL reconstruction
regressions, bootstrap intervals and within-receptor centering are descriptive
evidence; they were not removed merely because they estimate a relationship.
Class-II and template-coverage evidence remain separate.

## Retired models are distinct historical implementations

The benchmark local contact model fits independent available-pair binomial
probabilities with an intercept, residue fields, asymmetric pair/geometry terms,
L2=1, gauge projection and sandwich uncertainty. The library Potts CLI fits
coupled contact configurations using neighbour-contact covariates in a penalized
pseudolikelihood. These are distinct statistical models, not equivalent frontends.

Local Gaussian diagnostics estimate power-transform parameters, location/scale,
class means/covariances/priors, shrinkage, projection directions and conditional
confidence moments. Equivalence to shipped assess is not established. Native
Q/T covariance weights, q_coupled Pearson/cohort moments, confidence GLMs, fitted
transition scales, self-consistent contact predictions and fitted multivariate
MD scores are historical. None is a final scoring producer. Their APIs and
artifacts remain inspectable; presence in the package is not final eligibility.

FINAL_MODEL_ARTIFACTS.tsv pins unchanged shipped empirical/reference inputs and
historical learned bundles. The private benchmark inventory pins 242 local
model source/output files and 17 accepted result tables. All 540 files of the
existing manuscript archive verify unchanged. All metric gains, losses, counts
and historical source discrepancies remain; no result is relabelled fit-free.

## Execution and acceptance boundary

The default library reproduce.sh already delegates to the existing fit-free
driver. The preserved benchmark recompute.sh now requires --historical before
execution, and the library historical caller forwards it. Future five-block
peptide/receptor CLI calls explicitly select structural references, matching
the later manuscript restoration. Earlier virtual-reference receipts retain
their source/commands/results. No scientific analysis was rerun for this audit.
Twelve focused synthetic command/population tests pass; shell syntax passes.
The manuscript graph/macro gates are read-only; a deliberate fitted-artifact
injection is rejected. Ordinary paper builds typeset recorded data without
fitting or scoring. Final build/visual checks live in the manuscript QA record.

Issues 33 and 35 close by explicit final/historical disposition and graph evidence,
not by implementing a replacement fitted scorer or claiming independent fits
are CLI reproductions. Issue 32's historical decomposition gap and issue 26's
long MD campaign remain independent. Release approval and deterministic topology
reductions remain on the roadmap. No paid review or release publication was requested.
