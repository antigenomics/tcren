# Contact-penalty mechanism investigation — 18 October 2026

Author-requested mechanism investigation; complete on the unchanged benchmark rows. The unchanged 3808-model comparison and its
frozen resource remain authoritative inputs. Report chemical heuristic support
versus proximity-only contacts separately from continuous geometry equivalents;
neither defines physical truth. Keep every epitope, MHC class and template stratum.

Steps → verification:

1. Audit exact formulas and matrix/filter callers → reproduce signed contribution
   identities without fitting, score-direction changes or new candidate metrics.
2. Export only fresh geometry/chemistry diagnostics from the same raw coordinates,
   using bounded batches and the unchanged frozen resource → exact physical-edge,
   source-hash and weight agreement with the accepted contact decompositions.
3. CLI-owned table audit → retain every model (including zero contacts), attribute
   means/motifs and exact gained/lost binder–nonbinder rank pairs per cohort/table.
   Audit all VDJdb cases and CPL ipTM/original extremes; no pooled VDJdb accuracy.
4. Consult verified primary structural biology/chemistry literature → distinguish
   actual chemistry, statistical contact representations, AF confidence and inference.
5. Deliver full case records and mechanistic conclusions into the visible benchmark
   checkout; preserve every earlier gain/loss and all unrelated work.

Existing owners: `contact_observations` and v2 classifiers own chemistry;
`marginal_samples`/predictive tails own frozen geometry; accepted decomposition
owns every original/modified score; `associations._metrics` owns the accepted
metric definitions. New diagnostics expose lineage and arithmetic, not a new score.

## Source and mathematical audit

The implemented formulas match the declared four-arm menu. Let E₀ include the
unchanged framework terms and original omitted-entry convention. Over known-entry
assessed CDR edges, define U=Σu, Fᵤ=Σu max(−e,0), Bᵤ=Σu max(e,0). Then exactly:

- Efixed=E₀+U; Ematrix=E₀+MU.
- Efavourable=E₀+2Fᵤ; Eliteral=E₀+2Fᵤ−2Bᵤ.
- Eliteral−Efavourable=−2Bᵤ. Any literal-versus-favourable difference therefore
  comes exclusively from rewarding original positive entries.
- For positive/negative models, D=Enegative−Epositive is their favourable
  separation. Adding P changes D by Pnegative−Ppositive. A nonnegative penalty
  can systematically hurt real positives if their contact burden is larger.

Here e is a dimensionless frozen matrix entry and M its maximum finite absolute
value. These identities explain arithmetic, not physical causality. Positive
statistical log-odds entries are not automatically repulsive chemical energies.
The TCRen derivation uses observed joint counts relative to endpoint marginals;
unchanged entries after filtering preserve the original table, not necessarily
its interpretation as the selected population's log-odds. The literature note
provides the exact conditional-probability identity and verified source evidence.

The worst-supported reduction is u=1−min(sⱼ). Identical repeated orientations
are idempotent under minimum; distinct atom distances add opportunities for
extremes. Even an ideal calibrated single continuous predictive tail is uniform
rather than concentrated at one. Under an illustrative independent K-tail model,
E[u]=K/(K+1); this is mathematical context, not a calibration claim for correlated
contacts. As d→5 Å from below, log[d/(5−d)] diverges and its tail tends to zero;
exact d=5 transforms are missing. Atanh(cos angle) similarly has boundary effects
and discards torsion sign. The recorded geometry audit will establish actual
winner families/distances; none of these observations alone labels a contact
physically impossible.

Existing `scoring --drop-untyped` uses closest-atom v2 classification, whereas
`residue_pair_types` unions chemistry over every atom pair. The new diagnostics
quantify this difference without changing the old filter or baseline. v2 polar,
hydrogen-bond and cation–π flags retain their documented heavy-atom heuristics;
CB-frame orientations are not hydrogen-bond donor/acceptor angles.

Verified literature and access scope:
[CONTACT_MECHANISM_LITERATURE.md](CONTACT_MECHANISM_LITERATURE.md).

## Live execution and focused verification

The first six-model probe correctly failed on an unset installed-wheel data root;
no score was produced. The existing pinned MHC FASTA was present; the driver now
passes its explicit TCREN_HOME and records its hash rather than rebuilding it.
The corrected six-model probe completed in 6.66 s with peak sampled process-tree
RSS 1335181312 bytes. Its 111 CDR edges match 222 accepted potential/edge rows
with exact zero weight difference. The full fresh geometry/chemistry-only export
completed with one compute thread, 32-model batches, an 8 GiB sampled process-tree
guard and 1800 s cap. Accepted scores, resource and labels remain unchanged.

Independent review repaired two boundary issues before full launch: canonical
endpoint chemistry now swaps with physical IDs, and single-class strata preserve
unavailable AUC deltas. Reversal/zero-pair regressions and the zero-contact full
CLI audit pass. The reviewer rechecked both repairs. Completed focused gate:
56 passed, three slow tests deselected; new mechanism/plausibility gate eight
passed. Ruff and warning-fatal Sphinx passed after one heading-underline repair.

The source transfer's initial directory comparison produced a read-only deletion
preview of destination-only files; no deletion ran. Automatic review rejected the
bulk push, so only individual files were transferred with backups and tested from
a new immutable wheel. A default Python 3.13 build was not installed; the actual
wheel is explicitly Python 3.12. Earlier unused builds and failed probes are not
accepted scientific outputs.

## CPL and model-generation context

The benchmark source ledger identifies VDJdb negatives as de novo mock mismatched
TCR/pMHC pairings, not experimentally solved nonbinding structures. CPL's binary
endpoint is archive best/worst candidate membership derived from positional-scan
activation-library responses; it is not an individual binding-affinity assay for
every candidate. The accepted counts, labels and metrics remain unchanged.

Recorded original-2022-TCRen/ipTM AUROC extremes (0–1 scale): 1E6 0.972712/0.731101
(325 candidates), 4C6 0.955124/0.675427 (321), MEL5 0.973981/0.737350 (325),
868 0.537189/0.975388 (321), ILA1 0.348356/0.976594 (325), MEL8
0.689850/0.999111 (260). TCRen2 results remain separate in the accepted ledger.

This audit reread the twelve existing top views for the highest-ipTM negative
and lowest-ipTM positive examples of those six clones from
`results/contact-context-20261017-a/renders/`. Both labels visibly make groove/
peptide contacts; compact-looking networks are not confined to positive labels.
4C6's selected negative has a visibly denser local network than its selected
positive; ILA1's selected positive has a denser central network than its negative.
MEL5 and MEL8 show different loop paths despite both models engaging peptide.
These image observations concern twelve preselected discordant models, not the
whole cohort or an energetic/steric validation. The prior manual ledger owns
all image/scene hashes and its larger forty-four-view review; no new render or
coordinate score was substituted.

The matching native-motif investigation already records Tyr→Cys at GLC peptide
position 3 in native 3O4L. Its unchanged table values are +0.654505003725954
(TCRen) and +0.8270901793661077 (TCRen2). A positive empirical entry can therefore
occur in a native-compatible motif. The new weighted-positive decomposition will
quantify its contribution, without equating a positive log-odds with repulsion.
The prior fixed-template CPL reproduction and deposited-model ranking remain
distinct routes; recorded strong ILA1 threading performance is not substituted
for its current deposited-model score.

## Completed evidence and deliverable

All 3808 original raw sources succeeded. The fresh diagnostic stage completed
in 1016.07 s at 1446461440 bytes peak sampled process-tree RSS (1.347 GiB).
The table-only audit completed in 3.40 s at 774848512 bytes. It reused the requested
accepted final decompositions and did not recalculate benchmark scores.

Exactly 85669 assessed CDR edges match the accepted frozen weights with maximum
absolute difference zero; 2404 framework edges remain unassessed. Independent
verification checks all 88073 physical pairs, 7616 model/table rows, unchanged
matrix entries, complete raw/source hashes, labels, all zero-contact models and
372 exact accepted AUROC differences. Signed-sum errors are ≤2.85e-14.

Visible report and 31 case records:
`private-evaluation-record`.
`CASE_INTERPRETATIONS.md` gives individual notes for both assay epitopes, all
22 VDJdb cases and all seven CPL clones. Detailed contact/model examples remain
local; aggregate records do not publish restricted per-structure assay labels.

The mechanism is quantitative: original-2022-TCRen scored U means are
20.242582/16.578086 for GLC positives/negatives and 22.903728/17.091133 for YLQ.
Their similarly high per-edge weights cause nonnegative penalties to punish larger
positive networks. GLC literal flip rewards its larger Bᵤ (3.780629/2.740786),
with CDR1α Tyr→Cys3, CDR3α Ser→Thr4 and CDR1α Tyr→Thr4 among leading contrasts.
This is a signed statistical-score mechanism, not proof of attractive chemical
energy. Near-cutoff atom tails and angular extremes drive large weights in both
labels; chemical support and geometry equivalents are demonstrably different.

The complete VDJdb report distinguishes contact-burden reversals (e.g. RAQAPPPSW)
from composition changes at similar counts (e.g. PTDNYITTY/TPRVTGGGAM). CPL notes
trace positive/negative residue-pair contributions in the six ipTM/TCRen extremes
and preserve the separate fixed-template and deposited-model routes. Every gain
and loss is retained. No new weight, threshold, classifier, default or manuscript
result was promoted; a future reliability weight requires independent residue-edge
calibration and explicit conditional-potential interpretation.

Final review also closed raw-ledger/null-hash coverage and absent-class motif
contrast gaps, with focused regressions. Ruff, warning-fatal Sphinx, 56 focused
tests (three slow deselected) and the independent full-output checks passed.
