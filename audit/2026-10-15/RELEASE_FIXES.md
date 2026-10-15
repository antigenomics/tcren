# Fit-free release implementation

The author authorized follow-up fixes and fresh reproduction on 15 October 2026,
superseding the audit-only stopping point. Branch: `codex/fit-free-release-fixes`.

## Final-model boundary

Fixed empirical contact-potential tables are explicitly allowed, with their source
counts, pseudocounts, weighting, reference convention and table digest disclosed.
The existing contact sum and virtual reference difference remain the formulas:
Phi = sum over observed residue contacts of the selected potential entry;
DeltaPhi = Phi(observed identities) - Phi(reference identities), on the same map.
These are dimensionless statistical-potential sums, not calibrated binding energies.
Receptor reference replaces receptor identities; peptide reference replaces peptide
identities. They are distinct quantities. Lower Phi/DeltaPhi is the energy direction.

Learned Q/T covariance weights, cohort standardization, Gaussian/GLM readouts,
fitted Potts and transition links are excluded from final scoring. Their APIs,
parameters and recorded observations remain historical evidence. No parameter fit
or new scoring combination is authorized by this boundary.

## Implementation and verification plan

- [ ] Reconcile every open issue with its owning checkout and current source.
- [ ] Repair numerical/domain defects with focused synthetic regressions.
- [ ] Export exact receptor references and publication decompositions through CLI,
      with one batch annotation, complete input accounting and provenance.
- [ ] Reconcile final reproduction with the permitted formulas; preserve historical
      producers and results, and compare identical populations without changing floors.
- [ ] Verify mathematics, fast tests, lint, docs and independently installed wheel.
- [ ] Run all five authorized benchmark blocks fresh from raw inputs; record every
      metric gain/loss, counts and missingness. No resumed computed intermediates.
- [ ] Open one PR to master, request Copilot review, address findings and verify CI.
- [ ] Prepare the release at the human PyPI approval gate only after required checks.

The original dirty generated appendix and result files are not part of this change.
No complete reproduction, remote merge or release is established by this plan.
