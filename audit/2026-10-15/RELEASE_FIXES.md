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

- [x] Reconcile every open issue with its owning checkout and current source.
- [x] Repair numerical/domain defects with focused synthetic regressions.
- [x] Export exact receptor references and publication decompositions through CLI,
      with one batch annotation, complete input accounting and provenance.
- [x] Reconcile final reproduction with the permitted formulas; preserve historical
      producers and results, and compare identical populations without changing floors.
- [x] Verify mathematics, fast tests, lint, docs and independently installed wheel.
- [x] Run all five authorized benchmark blocks fresh from raw inputs; record every
      metric gain/loss, counts and missingness. No resumed computed intermediates.
- [x] Open one PR to master, request Copilot review, address findings and verify CI.
- [x] Release 4.0.0 after the author approved the PyPI gate; verify publication and local installation.

The original dirty reduction figure and notebook result file are preserved outside
these commits. The generated version label was refreshed with its owning script.
Issue-by-issue evidence and explicit remaining boundaries are in
[ISSUE_DISPOSITION.md](ISSUE_DISPOSITION.md).

## Verification

- Fast local suite: 1,046 passed, one network skip. Broader non-slow suite:
  1,120 passed, five skips, one stale generated version expectation; the owning
  generator corrected it and 46 focused catalogue/interpolation checks passed.
- Ruff and warning-as-error Sphinx passed. PR CI on Python 3.10/3.12, lint and
  docs passed at `e83a1e1` after all three Copilot findings were addressed.
- Numerical review retains the exact contact/reference identities, distinguishes
  raw peptide-reference sums from reference-weighted historical totals, uses
  Cholesky Gaussian evaluation, and masks undefined reference arithmetic.
- Reviewed 4.0.0 wheel and source distribution built. All packaged Python files
  matched source; native payloads were unchanged by the Copilot review fixes.

Full fresh-run results and final integration/release state are recorded below
when verified. Earlier partial attempts are diagnostic evidence only.

## Fresh reproduction

The reviewed-wheel run completed all five blocks from raw inputs. Matched counts:
2,102 measured CPL rows (one unmeasured source row explicitly excluded), 618 TCRvdb,
1,089 VDJdb, 58 class-I ATLAS and 187 MD records (12 raw labels lack coordinates).
Compared keys, labels and score availability agree exactly. No fitted readout is
included; generator confidence is reported separately as a diagnostic.

All 1,552 paired metric comparisons are retained. Of 1,504 classification comparisons,
94 increased, 98 decreased and 1,312 were exactly unchanged. Of 48 correlation
comparisons, 16 numerically increased, 20 decreased and 12 were unchanged; numerical
correlation direction is not automatically an improvement. VDJdb cases remain
separate by epitope and native-template coverage, without pooled strata means.

The largest classification decrease is template-scarce PTDNYITTY receptor CDR1/2
peptide-reference AP: 0.4674881462 → 0.3916269933 (delta −0.0758611529; 54 rows,
27 positive). The largest increase is native-covered GILGFVFTL in the same score's
precision at recall 0.1: 0.9545454545 → 1.0 (delta +0.0454545455; 46 rows, 23 positive).
The loss remains reported. Its fresh scores round exactly to the recorded six-decimal
baseline; higher output precision resolves historical ties (33 → 44 distinct scores).
No rounding or alternative metric was substituted to conceal the decrease.

### Limitations

Historical potential-table identity is not recoverable for every baseline arm;
current hashes do not establish historical lineage or a causal explanation for every
delta. This is complete reproduction of the five permitted fit-free benchmark blocks,
not the retired fitted paper or the separate unfinished long MD campaign. Original
baseline snapshots, exclusions, commands, hashes, all per-task metrics and failed
attempt records are retained in the benchmark's local run deliverable.

Final producing source: `cf103f9`. The final metadata-corrected wheel
(SHA256 `c2f58510cabdbd3687b136c82e653426e973655a6f74e56f7c662c6371b76550`)
completed another full raw-input run. All five fresh dataset tables and all 3,672
metric rows exactly match the preceding successful reproduction; the 1,552 paired
comparisons above therefore describe this final package as well. Python 3.10/3.12,
lint and warning-as-error documentation CI passed on this source revision.

## Published 4.0.0 verification — 15 October 2026

PyPI serves 4.0.0 wheels and the source distribution. Release tag `v4.0.0` and
master resolve to `3cdbe1b`; tests, documentation and publishing workflows completed
successfully. The published CPython 3.12 macOS ARM64 wheel has SHA256
`effa605ec1cf2ac0096fb173600c8f83bdd31812e9c74715746fe8a2e6af3f88`.
It was reinstalled noneditably in the local repository environment: package metadata
and `tcren info` report 4.0.0, CLI help works, and all 178 installed package payload
files match the wheel. Focused installed-wheel checks: 114 passed, one slow test
deselected, 30.19 seconds. Both synthetic MD plot exports were visually inspected.

Another authorized fresh raw-input run retained all five benchmark populations and
availability masks. CPL, ATLAS and MD dataset tables exactly match the reviewed
pre-release run. Whole-receptor references also match; CDR/CDR1–2/CDR3 references do
not. Of 1,836 matched metric comparisons against that reviewed run, 1,444 are exactly
unchanged, 164 numerically increase and 228 decrease. All changes are retained;
correlation sign/direction is not automatically an improvement. The largest precision
at recall 0.1 decrease is 1.0 to 0.666667 on TTDPSFLGRY (62 structures, 31 positive),
and the largest increase is 0.625 to 1.0 on YVLDHLIVV (46 structures, 23 positive).
Both are template-scarce regional receptor-reference arms. Against the original
historical inputs, 1,552 paired metrics comprise 200 numerical increases, 280
decreases and 1,072 exact ties. These supersede the earlier comparison counts only
for the newly verified annotation environment, without deleting the earlier record.

Read-only index inspection explains the discrepancy: the earlier environment used
an ARDA cache whose searched MMseqs databases did not match the raw reference FASTAs.
Its human database has 15,690 records versus 15,069 FASTA records, with 10,084
sequence mismatches over 14,515 shared scaffold names. Mouse has 16,641 mismatches
over 18,048 shared names. The current environment's databases match their human and
mouse FASTAs exactly. ARDA code, version, FASTA/markup hashes and MMseqs version alone
were identical and therefore failed to capture the searched database identity.
ARDA's ancestor-based source-root detection and timestamp-only index validation
explain the different index roots and acceptance of stale content. Neither index
was modified during this verification. Issue #30 remains open for this dependency
and actual-index provenance gap; the earlier regional results cannot certify stability.

The requested endpoint extension rebuilds CPL activation from the seven raw response
plates and uses the same run's authoritative CLI outputs. Per-clone CPL activation
correlations reproduce the historical unweighted reference-sum arm to below 1e-12;
the endpoint is log2(1 + S_CPL), not log fold change. The reduced rupture-work
association with log10(koff) remains R = -0.5397988480 over 23 class-I WT crystals.
The historical three-interface MD correspondence remains R² = 0.4598920715 over
184 training structures, with three original held-out structures retained. This is
descriptive correspondence, not an ergodicity proof or a fitted final scorer.
The supplied delta_sum/ddG correlation remains R = 0.8006247560 over all 181 rows;
conflicting source descriptions still prevent assigning experimental/free-energy
provenance. Statistical reference differences retain their own units and convention.

Closed with verified release linkage: #14–20, #22–25, #27–29, #31 and #36 (16 issues).
Remaining open: #21, #26, #30 and #32–35 (seven). Final publication migration,
historical CLI non-equivalence and the long MD reconstruction are not completed by
this release check. No manuscript claims, final formula, fitted parameters or
confidence-fusion rule were replaced.

## Searched-index repair and receptor reproduction — 15 October 2026

The authorized follow-up repairs ARDA's root selection and stale-index acceptance.
ARDA PR [148](https://github.com/antigenomics/arda/pull/148), producing revision
`97159efbecc80625d3eb4d80f7c82339cfb19aea`, merged as
`2c9962674b4a2c44e0c826ffdb55c856cc4ca82c` after the author's explicit admin-bypass
authorization and successful CI verification. It recognizes a source checkout only when
its `src/arda/paths.py` resolves to the loaded module. Explicit `ARDA_HOME` remains
supported. Reference indexes now use content-validated immutable generations keyed
by source FASTA, database type and MMseqs version, with sequence, header, index and
split-payload hashes checked before reuse. Legacy indexes remain preserved.

tcren PR [38](https://github.com/antigenomics/tcren/pull/38) merged as
`b593c421c4dee280382ec592b3819988e9160dd7`. Annotation captures actual ARDA and MHC
searched-index identities before and after a search, rejects changes during search
and mixed generations across chunks, and carries the identities through scoring,
features, recognition and exports. The geometry scoring path records both annotation
passes; recognition from features inherits the input stamp. MHC hashes identify the
searched bytes; they do not assert correspondence to an MHC source FASTA. Missing
backend identity is explicitly unavailable rather than certified.

The installed candidate versions are tcren 4.0.1 and ARDA 2.36.2; these patch versions
have not been published to PyPI. Their locked wheel SHA256 values are respectively
`32e0559325e997aba4c1db379a6d9a12a4317b46164c94883cec654cdf61e4a7`
and `3e6d0700b72489050325e22b306f93108100a41e78b42b31ce818c6879ede5cf`.
All 178 tcren and 63 ARDA installed payload files match their producing wheels.
Installed-wheel tcren checks: 1,132 passed, six skipped, 141 slow tests deselected;
ARDA checks: 1,350 passed, 16 skipped. Both Ruff and warning-as-error Sphinx gates
passed, as did tcren Python 3.10/3.12 and ARDA pull-request CI.

Direct parsing of repaired MMseqs sequence/header/index payloads confirms identical
FASTA and searched name sets with zero sequence mismatches: human 15,069 records,
mouse 18,800 records. The independent fresh receptor runs retain this semantic
validation, locked wheels, commands, raw-input identities and searched-index hashes
in their private benchmark deliverable. Geometry, interface and topology are direct
descriptor arms; fitted readouts and both banned pitch fields remain excluded.

The affected receptor reproduction retains 618 TCRvdb models (309 positive) and
1,089 VDJdb models (523 positive), with both TCRvdb epitopes and all 22 VDJdb
epitope/template-coverage cases separate. It evaluates 124 eligible direct
placement/interface/topology descriptors, six energetic arms, 12 exact receptor
reference arms and two native confidences. There are no missing/unexpected CLI rows,
annotation errors, changed labels or baseline key-membership differences. Per-arm
undefined descriptors retain their actual finite support.

All 3,000 matched metric rows against the recorded published-4.0.0 route are exact
ties, with exactly equal raw values for every shared same-definition score. Against
the earlier stale-index attempt, 3,000 comparable rows contain 2,216 ties, 328
numerical increases and 456 decreases. Against the historical feature snapshots,
26,000 same-name metric rows contain 25,728 ties, 120 increases, 144 decreases and
eight undefined correlations for constant GLCTLVAML contact-type quantities. These
are repeated metric rows across endpoints and separately matched ipTM/pLDDT
comparisons, not independent experiments. No arm becomes newly finite or unavailable.
Every signed delta, finite count and unavailable metric is retained in the benchmark's
`results/receptor-index-20261015-a3/metrics.csv` and `availability.csv`.

All combined-interface `g_*` contact-graph and contacted-pMHC `fp_*` footprint values
match the historical feature snapshots exactly; `m_*` spectral differences are at
most 1.5543122344752192e-15. The largest other feature difference is 2.3874e-12;
floating-point changes can resolve old ranking ties. Historical feature snapshots
lack verified catalogue/annotation lineage, so their same-name comparisons remain
cross-version diagnostics rather than causal estimates of this repair.

Topology is retained explicitly alongside placement and chemistry. The combined
CDR-to-pMHC graph, pMHC footprint complex and CDR-to-peptide/MHC-helix soft proximity
spectra are distinct objects. As exploratory examples among 124 descriptors,
`m_erank_tm` ROC-AUC is 0.861677 versus ipTM 0.776993 on the same 423 YLQPRTFLL
models, and 0.695062 versus 0.813580 on the same 195 GLCTLVAML models. The existing
CDR3-length-dependency warning remains attached to this descriptor. RAQAPPPSW
`g_loop_overlap` has ROC-AUC 0.754630 versus matched ipTM 0.087963 on 33 supported
models (nine positive, 24 negative); the 11 unsupported models are all positive and
remain explicitly unavailable. No full-cohort graph-performance claim follows.
These observations support case-specific inspection, with chemistry, placement,
topology and native confidence reported separately; no sign optimization, threshold
or new scalar confidence fusion is introduced.

Both authorized raw-input reconstructions completed with identical population
accounting, labels and finite support. The independently rebuilt ARDA indexes have
identical source and searched payload hashes across all six cohort/CLI stages
(18 stage/index identities including the common MHC input). MHC was a shared,
hash-recorded reference input, not an independent source-validation experiment.
Both cohorts' energy and exact receptor-reference tables agree exactly. All
highlighted `g_*`, `fp_*`, `m_*` and `D2_pep24` values and metrics also agree exactly.

Whole-output bitwise reproduction does **not** pass: four other topology columns,
`co_pep`, `co_mhc`, `partcoef_tcr` and `partcoef_pmhc`, differ by at most
1.1102230246251565e-16. Source inspection identifies unordered Polars group
reductions followed by floating means in `topology/literature.py`; summation is
not associative. Across 114,400 repeated metric rows, 448 current values differ,
112 historical delta fields differ, and 40 decrease flags differ. The largest
metric delta is template-scarce ATDALMTGF `co_pep` precision at recall 0.1:
0.631578947368421 to 0.6111111111111112 (delta -0.020467836257309857; 42 models,
21 positive). No score was rounded to hide tie breaking. All 864 availability rows
retain their finite counts; only three historical maximum-difference summaries
change, by at most 2.7756e-17. Complete value/metric/availability differences and
the false strict-exact-reconstruction flag remain in the benchmark deliverable.
This separate numerical limit does not invalidate the verified index repair or
certify all topology rankings as bitwise stable. It was recorded without another
library change or scientific rerun.

Issue [#30 was closed with the complete verification and numerical-limit record](https://github.com/antigenomics/tcren/issues/30#issuecomment-6193696693).
Both merged master CI/documentation workflows passed. Open issues are now #21,
#26 and #32–35; this closure does not complete manuscript/figure migration, restore
missing historical lineage or establish bitwise stability of the four literature
topology columns. The private benchmark report supplies task-specific descriptor
recommendations and links every raw-value, metric and availability delta.
