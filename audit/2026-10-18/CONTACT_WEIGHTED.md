# Multiplicity-calibrated weighted contact scoring — 18 October 2026

Author scope: calibrate contact weights on native/external structures, freeze them,
compare weighted/original TCRen and ipTM on unchanged rows, and provide comprehensive
VDJdb/CPL interpretation with R plots. No default is promoted and no manuscript
standing finding is rewritten. The original plan's reference-informed geometry
remains the scientific context; this extension addresses the measured multiplicity
bias in the previous worst-tail reduction.

## Frozen definition

Keep the approved native/external transformed marginal resource unchanged:
SHA256 `63c141f0a8be2e090cda543193148dbd8b06f937dcaf5f4f428124c77d013f16`.
For each physical residue edge, let `s` be its minimum available predictive tail.
Group references by supported atom-distance count and distinct `(key,x)` omega,
theta, phi term counts. Thus repeated angular copies over atom rows count once.
Native CDR:peptide references define the target; only native-empty groups use
external antibody:antigen contacts. Native TCR:MHC and peptide:MHC edges are
excluded from the rank calibration. The marginal model itself is unchanged.

For a selected group of `N` reference edges, freeze the inclusive empirical rank
`w = (1 + count(s_reference <= s_query))/(N + 1)`. Ties retain their full inclusive
rank; the added one is the query's own rank contribution. There is no tuned
threshold, floor, interpolation, scale or assay-label fit. Native groups have
precedence at every positive count. Missing geometry or unseen groups retain
neutral weight one and explicit status. Export group edge/structure counts and
source, allowing coarse or transferred calibration to remain visible.

New CDR contribution is `w*e`, where `e` is the unchanged original table entry.
Framework remains unweighted; unknown entries retain legacy omission and counts.
Both positive and favourable magnitudes attenuate. Defining
`F_lost=Σ(1−w)max(−e,0)` and `B_lost=Σ(1−w)max(e,0)` gives the exact identity
`E_weighted = E_original + F_lost − B_lost`. This distinguishes gains driven by
removing positive entries from losses driven by removing favourable entries.

## Reference evidence and execution

The requested final native/external atom observations are calibration inputs;
no hidden cached feature rows or score-time updates are used. Native calibration:
63,849 atom pairs, 7,526 residue edges from 355 structures. External calibration:
854,947 atom pairs, 101,795 edges from 1,831 structures. Both have complete
minimum-tail support. There are 214 source-specific opportunity groups.

Frozen resource SHA256:
`c096b37a19f5937ce5d40afa4a0c84d4caf8d438a4dc32203e0774d05ebaf07e`.
The timestamped calibration record precedes benchmark scoring and binds exact
native/external inputs, marginal resource, installed noneditable wheel and producer.
Reference calibration ran in 3.532 seconds with 2,064,384,000-byte sampled peak RSS.
Scoring runs fresh raw models in one-thread batches of 32 with an 8-GiB/2400-second
process-tree guard. NumPy searchsorted within Polars opportunity groups sufficed;
no new native code, outer worker pools or dependencies were introduced.

## Interpretation and limitations

Weights are conditional reference-compatibility ranks, not probabilities of real
contacts, naturalness, occupancy, binding or free energy. Opportunity conditioning
addresses the declared count effect without claiming independent atoms or complete
exchangeability across chemistry, structural source or interface. Native calibration
includes available MHC-I/II CDR:peptide contacts; the benchmark comparisons remain
MHC-I only. Reference observations helped estimate marginals, so in-reference ranks
are not independent or conformal validation. Sparsity, external fallback and neutral
missing support remain explicit. Fixed matrix entries are numerically preserved,
but weighting is a new scoring operation rather than a conditional-potential refit.
Benchmark labels have already informed development. VDJdb mock mismatches are not
experimentally solved nonbinding structures; CPL library extremes are not individual
peptide affinity measurements. Predicted contacts alone do not establish AF's internal
causal choices. Chemistry heuristics and Cβ orientations do not provide a full angular
hydrogen-bond, solvation or steric validation model. Source-grounded interpretation:
[CONTACT_MECHANISM_LITERATURE.md](CONTACT_MECHANISM_LITERATURE.md).

## Delivery and checks

Completed: all 3808 models, 7616 model–potential rows and 88073 physical contacts
are retained. Original matrices/scores/graphs/raw identities and every comparison
row are exact. All 124 recorded arm/case/table AUROC differences pass pair-credit
accounting; all 744 metric rows have zero excluded models. Maximum signed/sum
identity error is 1.42e-14. Six R PDF/PNG figure pages were visually inspected.

Original 2022 TCRen AUROC: GLC .360013→.390772, YLQ .593062→.640932,
ILA1 .348356→.458946; losses include 4C6 .955124→.876630. TCRen2 has material
losses on SB27 (.861510→.726417), LLAGIGTVPI (.707500→.537500) and
SLYNTVATL (.538961→.386364). Both matrices and every case remain in the report.
Weighting attenuates statistical signs; it does not reverse them or uniformly
repair empirical-matrix inversions. No standing manuscript claim is changed.

Fresh scoring: 1142.642 s / 1505214464 bytes sampled tree RSS. Comparison:
2.137 s / 344539136 bytes; table audit: 3.422 s / 680591360 bytes. All guards
exit zero. Forty-two completed remote files match local hashes. Exact producer
sources from scoring wheel-a and diagnostic wheel-c are archived separately.
The latter additionally makes an entirely missing geometry observation explicit
neutral weight one; every assessed edge in this benchmark already had geometry.
Unused wheel-b is not a scientific producer. No completed stage was resumed.

Verification: 19 focused tests; Ruff src/tests; warning-fatal Sphinx; independent
mathematical/code review; transfer/source identity checks. An optional wider scripts
lint found eight existing issues in unrelated scripts; those files remain untouched.
Benchmark deliverable: `results/contact-weighted-20261018/REPORT.md`, all 31 linked
case records, long metrics, signed decomposition/support/motif tables and R plots.
