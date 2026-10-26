---
name: tcren
description: Use TCRen to annotate TCR–pMHC structures, rank candidate receptors or peptides, export contact evidence, and reproduce structural analyses through its CLI.
---

# TCRen

Use the installed CLI or public API for calculations; keep dataset selection,
private identity maps and publication formatting in their owning projects.
Start with `tcren --help` and the chosen command's `--help`. Do not reconstruct
scoring formulas in notebooks or manuscript scripts.

## Choose the task

| Input and question | Command | Guide |
|---|---|---|
| One pMHC, candidate receptors | `screen` | [Screen contract](../../docs/screening.rst) |
| One receptor, separately modelled peptides | `screen` | [Screen contract](../../docs/screening.rst) |
| Candidate peptides threaded on one fixed template | `score` | [Getting started](../../docs/getting-started.rst) |
| Existing structural epitope-energy/ipTM table | `rank-epitopes` | [Epitope ranking](../../docs/ranking.rst) |
| Contacts, atom inventories or score attribution | `contacts`, `export` | [Features](../../docs/features.rst) |
| Canonical placement and pMHC surface | `superimpose`, `surface` | [Inspection](../../docs/kit.rst) |
| Recorded observations and matched comparisons | `evaluate`, `associate` | [Evaluation](../../docs/evaluation.rst) |
| Dynamics and mechanical descriptors | `md`, `mechanics` | [MD](../../docs/md-validation.rst) |

## Install and batch inputs

Use the repository `setup.sh` for development, with its `.venv`. Installed wheels
contain native kernels. `tcren build-mhc-ref` prepares the MHC reference; ARDA owns
its TCR reference/search bootstrap. See [installation](../../docs/getting-started.rst).

Structure inputs accept PDB/mmCIF, gzip, directories, tar.gz, globs and manifests.
Pass the whole query set in one call: annotation is batched internally. Use the
command's CPU budget (`-t 0` for all allocated CPUs where supported, `1` for serial,
positive N for a bound). Never launch one search process per structure or nest
an outer pool around the search engine. See [performance](../../docs/performance.rst).

## Default TCR ranking

```sh
tcren screen -s models.tar.gz --confidence confidence.tsv -o scores.tsv
```

`confidence.tsv` has unique `complex.id`, `iptm` on [0,1], and `plddt` on [0,100].
Identifiers must match the input structure stems. Metadata can supply confidence;
TCRen does not calculate ipTM. Do not substitute experimental B factors for pLDDT.

The frozen rule is `S=(1-.41q-.22(1-q))*ipTM+.41q*B+.22(1-q)*M`:

- B: recognition graph beta1/E over CDR:peptide, peptide:MHC and CDR:MHC contacts.
- M: all-six-CDR/peptide local contact-cycle evidence calibrated against 2,056
  whole mock graphs. Frozen histograms and compressed calibration ship in the package.
- q: inclusive CDF of unlabelled cohort mean pLDDT/100 against five frozen means.

No query labels or score-time fit are used. Ambiguous outcome rows remain in
cohort context. Missing active values use flagged frozen imputation; changing
cohort composition can change weights. Coefficients were selected during labelled
development, so the score is not fit-free or a calibrated binding probability.

Read `screen_score`, `score_unavailable`, `missing_*`, component values and
`gate_q`/component weights together. Higher score ranks first. Zero-contact models
retain rows/reasons but no score; they do not count toward the default 30-model
floor. Flags report by default. `--weight` preserves the Betti:motif ratio while
changing the confidence weight. Epitope and legacy references have distinct rules.
Template status is context, not a new default coefficient selection.

APIs: `screen.screen_scores`, `contacts.motif.bundled_enrichment`,
`contacts.motif.contact_enrichment`. Reference equations and output fields live
in [screening](../../docs/screening.rst), not in duplicated notebook code.

## Confidence context for individual complexes

The measured cohort mean remains default. The explicit `motif-confidence-proxy`
command offers frozen MHC allele-context lookup with support counts and a flagged
global fallback. See [confidence options](../../docs/motif-confidence.rst) for
single-structure and small-set alternatives, their input requirements and limits.
For raw singleton structures use `screen --task tcr --confidence-context nearest-plddt
--min-cohort 1 --confidence confidence.tsv`; `nearest-mean-plddt` uses a small
same-pMHC set's observed mean. These modes select a frozen representative, report
observed support/distance and withhold when pLDDT is missing. They do not predict
an independently validated true cohort mean. Do not describe a representative allele match as experimental MHC identity.

## Structural and chemical conventions

CDR3 excludes conserved junction anchors. Preserve paired receptor, peptide and
MHC identities when joining results; an annotated sequence is not a binding label.
MHC role/groove annotation is needed by the graph; an exact allele name is not a
term in the default score. Sequence-equivalent allele names may be unresolved. Read [MHC identity](../../docs/mhc-annotation.rst)
for supplied metadata, sequence support and ambiguity fields.

Use `scoring --delta --reference-mode structural` for peptide alanine truncation
and `export --quantity receptor-reference --reference-mode structural` for
contacted receptor residues. These rebuild the contact map. Fixed-map virtual
substitution and the historical feature-table dPhi columns are different quantities.
Do not silently switch them during reproduction.

Chemical v2 flags are multi-label and use all atomic witnesses; the nearest
atom pair can hide other interactions. Hydrogen-bond flags are heavy-atom-compatible
contacts, not explicit-hydrogen/angle verification. Preserve backbone/side-chain
scope and count atom pairs, residue pairs and structures separately.

## Reproduce and inspect

[Reproduction](../../docs/reproduction.rst) owns CLI input/output provenance.
Record version, source/reference hashes, exact memberships and masks. Compare
scores and metrics on identical populations, retaining both gains and losses.
Raw-coordinate reproduction and frozen-table replay are distinct checks.

Historical `recognize`, `assess`, Potts and fitted contact commands remain available
for their existing contracts; they are not the default screen. Read
[reference dependencies](../../docs/fit-free.rst), [experimental ranking APIs](../../docs/experimental-ranking.rst),
and the relevant module docs only when reproducing those methods. `SOURCES`
records bundled parameter identities. Preserve the original classic potential
separately from updated classic tables and rederived/full/LOO/LCO variants.

## Verify changes

Run focused `.venv/bin/python -m pytest`, then affected broader tests. Use
`.venv/bin/ruff check src tests` and, with the venv active,
`make -C docs html SPHINXOPTS=-W`. Rebuild after native-source changes and verify an
installed wheel outside the checkout. Scientific reproduction is separate from
unit tests; never infer a publication or dataset upload from local build success.
