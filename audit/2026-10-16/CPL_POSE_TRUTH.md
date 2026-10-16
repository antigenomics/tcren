# CPL structural-reference validation

The follow-up audit recovers historical DockQ interface definitions and verifies
new candidate/reference identities before interpreting spatial diagnostics.
Existing `export --quantity pmhc-cells` supplies authoritative cells and residue
markup; the external DockQ 2.1.3 CLI supplies structural-reference measurements.
No scoring, descriptor, contact-potential or classification definition changed.

Historical source tables cover 1,108 models / 6,648 interfaces. Every recovered
Total_DockQ matches the deposit and the six-interface arithmetic mean to reported
precision. Its interfaces include TCR internal and peptide:MHC contacts; the four
recognition interfaces must be reported separately. A mismatched peptide or
substituted receptor reference measures reference deviation, not coordinate truth
for the candidate's own complex. Missing original coordinate bytes/executable
identity and unrecovered interfaces remain explicit provenance gaps.

The fresh subset comprises 17 raw structures, ten model/reference comparisons and
60 interface records. Exactly three candidates pass the observed-sequence gate,
one per clone. Peptide identity, physical chain roles and observed receptor/MHC
substitutions/internal gaps are verified; terminal expression residues and
constant-domain extensions remain recorded. The other comparisons retain their
identity exclusions without removing assay labels or ranking results. Dataset
prose is not authoritative membership: the later-added cognate is already in the
accepted statistics, and remains identified as such.

All 500 paired descriptor/ipTM association records use exact finite rows within
clone, class and template stratum; Pearson/Spearman estimates, missingness,
weaker associations and sign changes remain. The independent oracle checks raw
identities, interface arithmetic, matched membership, fresh spatial equality and
unchanged original classification exports. No fitted score, threshold, permutation
or new control is introduced. Three singleton checks cannot estimate detector
sensitivity/specificity; multiple independently verified coordinate outcomes per
intended stratum remain necessary. No VDJdb/TCRvdb/class-II truth is inferred.

The isolated external DockQ environment respects its declared NumPy<2 dependency;
the immutable TCRen wheel/environment was not changed. One annotation batch uses
at most four MMseqs CPUs; subsequent kernels, DockQ and R are serial. The complete
subset plus R summaries took 15.926 s, with sampled peak process-tree RSS 1.383 GB
(decimal), 0.1-s sampling and five transient unreadable samples. Peak 19 live
threads includes bounded auxiliary I/O pools; no outer annotation worker pool.
Installation, later independent QA and publication compilation are outside that
measurement. The 1,200-s / 24-GB target passed. All attempts remain recorded.

Benchmark `CPL_POSE_TRUTH.md` owns private source paths, exact mappings and full
precision results. Manuscript `POSE_VALIDATION_AUDIT.md` and `POSE_FINDINGS.txt`
own the publication record and remaining gates. Ten tables are generated from
CLI/DockQ exports through R and a formatting-only Python producer. Standing
fit-free findings and regressions remain. Stratos MD was untouched and uninspected.

Identifiability follows from the existing definition: fixed peptide/MHC coordinates
give identical pMHC covariances under any rigid TCR movement. Those covariances
can supply context but cannot by themselves detect the receptor displacement;
existing receptor-relative coverage/gap/clash/placement diagnostics are needed.
No new control was run. Fresh/deposited total disagreements are retained as
reference/configuration differences, never silently promoted to score regressions.
