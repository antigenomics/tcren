# Roadmap and example audit — 16 October 2026

ISALGO coding review and manuscript figure guidance were applied to the current
software, actual callers, documented commands, live GitHub work items and recorded
publication evidence. This is a bounded example/continuity audit, not a claim that
all modules, optional renderers or scientific findings were revalidated.

## Identity and execution

Source master: `c04cb943c60c76f4307ffcbe6003e39f9e9dd381`.
Installed package: noneditable tcren 4.0.1, ARDA 2.36.2. All **180** installed
package payload files match the immutable wheel, and every packaged Python file
matches source. Wheel SHA256:
`8367bee127271f61dffb0f76d64d6479f9a30a5a7d39d2f446571358cc302961`.
No environment replacement, global dependency installation or parameter fitting.

Public Native2026 reference inputs: 1AO7, class I, LLFGYPVYV; 1FYT, class II,
PKYVKQNTLKLAT. A two-file manifest supplied both structures in one call for each
export and feature command. The annotation reference was already present. The
first sandboxed attempt could not acquire ARDA's external reference-index lock;
it is not a scientific failure. The subsequent authorized run used the ordinary
environment without user-added thread-limit or data-root overrides.

Full local outputs, exact argv, logs, annotation/potential sidecars and both
rendered previews: `/tmp/tcren-example-audit-20261016/`. No completed benchmark
was rerun. Fresh calculations were limited to these requested examples and
explicit invalid-input probes; outputs are audit deliverables, not a resume cache.

| Example | Result / observations | Wall time, s |
|---|---|---:|
| `scoring --delta --reference-mode structural --no-superimpose` | Exit 0; two rows; mode explicitly structural | 3.686 |
| `export --quantity pmhc` | Exit 0; two rows; class I and II retained separately | 7.078 |
| `export --quantity receptor-reference --reference-mode structural` | Exit 0; two rows | 7.117 |
| `export --quantity position` | Exit 0; 16 contacted-position rows | 6.824 |
| `export --quantity residue` | Exit 0; 49 contacted-residue rows | 6.946 |
| `export --quantity contacts` | Exit 0; 53 TCR:peptide residue-pair rows | 6.789 |
| `export --quantity placement` | Exit 0; two rows | 6.751 |
| `score`, three nine-residue candidates, 1AO7 | Exit 0; exactly three candidate rows | 1.868 |
| `cpl`, 1AO7 | Exit 0; 9 positions × 20 amino acids = 180 rows | 2.958 |
| `surface --svg --cells`, 1AO7 | Exit 0; one map and 1,820 occupied-cell rows (64 × 32 grid) | 2.937 |
| `features -i placement,interface,topology --no-metadata -t 1` | Exit 0; two rows, 128 columns including identifier | 6.049 |
| `surface --svg`, 1FYT | Exit 0; one class-II map | 2.915 |
| `annotate`, 1FYT | Exit 0; 597 region-markup rows | 2.920 |
| `export --quantity pmhc`, nonexistent input | Expected exit 1; accounting/diagnostics retained | 0.595 |
| `surface --grid 0x0`, valid 1AO7 input | **Wrong exit 0**, one error row, zero mapped structures | 2.914 |

These are whole CLI process timings on this host, including Python startup and
annotation. They do not isolate contact-table scoring or establish multicore
scaling. Export adds physical-residue/provenance work, and pMHC export includes
SASA. Do not infer a kernel performance regression from these unmatched timings.
The vectorized table/native APIs already exist; no new parallel implementation
or computed cache was needed for the audit.

## Independent checks and manual inspection

- Read both scoring rows and both pMHC feature rows, candidate scores, CPL cells,
  region markup, contacts, position and residue exports. Every export accounting
  table reports two successful inputs and no failed input. Physical author residue
  IDs and insertion-code columns are present. Candidate identity is unchanged.
- 1AO7 raw TP Phi = −1.020363318891668; position contributions sum to
  −1.0203633188916683. 1FYT raw TP Phi = −0.7241452366597441; position contributions
  sum to −0.724145236659744. These are dimensionless selected-potential sums,
  one structure per comparison. All 53 contact rows independently join the exact
  exported potential table and recover the same scores.
- Endpoint residue contributions sum to −2.040726637783336 and
  −1.448290473319488, respectively: twice the corresponding interface sums because
  each contact contributes to both endpoints. Eight peptide positions per structure
  are contacted; the position export is not a list of every resolved peptide residue.
- Receptor-stripped peptide SASA is 419.8253519617044 Å² (1AO7, nine residues) and
  744.5810711111652 Å² (1FYT, 13 residues). Their exposed-area-weighted hydropathy
  is −0.12126494972867548 and −1.926183414886936 in Kyte–Doolittle scale units.
  These are two distinct example measurements, not a class comparison.
- Class-II annotation retains TRA and TRB FR1–4/CDR1–3, both groove-chain roles,
  helix/floor markup and the peptide. The scoring allele representative/candidates
  remain ambiguous sequence matches, not experimental genotyping. In particular,
  a representative DRA allele alone is not a complete class-II restriction.
- All 180 CPL cells satisfy `phi = phi_tcr + phi_mhc`; the nine WT cells have
  `effect_wild_type = 0`. These are fixed-template statistical-potential effects,
  not measured activation or a reproduction of the later seven-clone experiment.
- Opened both SVG-derived PNG previews individually at native 520 × 760 pixels.
  Peptide identity, orientation text and maps are visible without clipping. Found
  the presentation limitations below; these native example surface rasters are
  not manuscript figures. No manuscript heatmap was added or restored.
- Focused installed-wheel/API suite: **33 passed**, 3.44 s:
  `tests/unit/test_exports.py`, `test_associations.py`, `test_cli_entrypoints.py`,
  `test_prediction.py`, `-m 'not slow' -q --tb=short -p no:cacheprovider`.
  This includes paired finite rows, class/template separation, missingness,
  association CLI/provenance/overwrite handling and early interactive-parallel
  errors. Synthetic cases are software checks, not new scientific controls.

## Confirmed remaining software and presentation work

| Priority | Evidence / consequence | Smallest owning follow-up |
|---|---|---|
| P1 | `src/tcren/cli.py:surface` catches an invalid-grid `ZeroDivisionError`, writes an error-only CSV and returns success. The initial external-lock failure likewise returned 0 with zero maps. A pipeline can accept an all-failed run. | Validate positive grid dimensions before annotation; preserve failures and exit nonzero when no map succeeds. Cover both failed and ordinary batches. |
| P1 | `surface` calls `classify_chains` inside its structure loop, whereas `annotation.batch` and `exports.export_batch` already provide batched annotation. It has no CPU-budget flag and MHC search explicitly requests all cores. | Reuse the owning batch API; one query set per species/MHC, explicit budget and complete input accounting. Do not wrap per-structure searches in another pool. |
| P1 | Shared `_typed_batch`/`annotation.batch.iter_typed` and descriptor defaults still select all annotation cores. `recognition_table` can also spawn a `ProcessPoolExecutor`; runtime policy now excludes executor APIs. Passing descriptor `threads=1` alone is not evidence of a one-core process tree. | Trace propagation from CLI through annotation/native pools before a scoped fix; verify actual process-tree threads and memory. This audit did not measure scaling or assert CPU compliance. |
| P2 | Both native surface previews lack a visible colour-channel/unit label and source-outline legend. The height colour range is fit independently per map. Signed-channel `_normaliser` ignores explicit `vmin`/`vmax` despite the renderer's shared-range documentation. | Label quantity/units and outlines; honor explicit shared limits with a focused signed-channel check. Use the established publication backend and flat tables for the manuscript. |
| P3 | `export --help` still begins “Exact virtual receptor reference” although structural mode is supported. | Describe both modes; retain historical virtual default and explicit manuscript flags. |

The all-failed, batching and budget findings are tracked together in
[issue #43](https://github.com/antigenomics/tcren/issues/43). These are audited
findings, not implemented runtime fixes. The shared callers and
existing APIs above define the scope for a focused follow-up. Neither numerical
scores nor regression flags were changed to produce this report.

## Reconciled publication / issue state

Initial live checkpoint on 16 October: six open issues (#21, #26, #32–35), no open PR; latest observed
GitHub release 4.0.0. Master tests/docs workflows passed at `c04cb94`. No new paid
review, release approval or public-install claim is implied.

`ROADMAP.md` now separates completed reproduction/restoration from real remaining
work: release of later fixes; exact CLI-to-render lineage (#34); final producer
and archive disposition (#33/#35); historical Potts CLI non-equivalence (#32);
notebook-safe multicore backend (#21); long MD reconstruction (#26). The last two
are distinct from the supported serial/guarded workflow and existing short MD
exports. Four literature-topology columns' recorded tie-sensitive nondeterminism
remains visible and is not fixed by changing printed precision.

The recorded benchmark/manuscript heads are `8657b89` and `26b912c`; their output
and source hashes are authoritative for historical claims. The original procedure
on later CPL data and the migration exceptions are already reproduced, not pending.
Earlier audit discovery notes are explicitly superseded by their follow-up report.
README task entrypoints now lead with fit-free CLI routes and label fitted APIs
historical. No established result, archived model or regression was removed.

Unverified here: optional PyMOL/DOPE engines, every historical notebook, every CLI
subcommand, first-use reference bootstrap on a clean machine, multicore scaling,
long trajectories, publication data deposits and a new public package release.
The two pre-existing dirty research outputs were preserved unchanged.

## Documentation verification

Ruff over `src tests` and warning-as-error Sphinx HTML passed after the documentation
reconciliation. Source/model files were not changed; no broader scientific rerun
was required. `example_commands.tsv` records all 15 exact CLI argument lists with
portable `$CLI`, `$RAW` and `$OUT` placeholders, exit codes and measured times.
`example_inputs.tsv` records both raw input hashes. These fixtures are public PDB
references; machine-specific annotation/provenance paths stay in the local outputs.

To replay the receipt, set `CLI` to the intended installed executable, `RAW` to the
public Native2026 mirror, and `OUT` to a new output directory. Write `structures.txt`
with the absolute paths of `1ao7.pdb.gz` and `1fyt.pdb.gz`, one per line. Write
`candidates.txt` with `LLFGYPVYV`, `ALFGYPVYV` and `AAAAAAAAA`, one per line. Run
the receipt commands in order; the failure probes must be evaluated against the
expected statuses above. Convert each native SVG using `rsvg-convert` for inspection.
