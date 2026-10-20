# Native contact observations — 18 October 2026

Stage 1 of `CONTACT_PLAUSIBILITY_PLAN.md` is implemented and its new native
summary is complete. This establishes the observation/export layer, not a
Bayesian prior, unnaturalness rule, penalty or classifier improvement.

## Implementation and verification

`tcren export --quantity contact-atoms --interface all` reuses the existing
batched TCR/MHC annotation and cKDTree neighbour search. The new opt-in quantity
retains physical atom/author residue/insertion identities, coordinate/altloc
context, whole unordered amino-acid/atom tuples, observed C-beta distance and
omega/theta/phi geometry. It restricts MHC endpoints to annotated helices and TCR
endpoints to CDR1/2/3. Missing/degenerate orientation geometry remains null.
The documented deterministic single-conformer convention applies only to this
quantity. Compressed Parquet avoids a large plain-text contact export.

Summaries separate atom-pair from residue-pair frequencies and retain interface,
MHC class and receptor region. Residue chemistry uses unions over all atom pairs;
closest-atom distance is summarized separately. Existing chemical heuristics and
ring stacking are reused. Numerical arrays and table operations sufficed; no new
dependency, native kernel or parallel worker pool was added. Legacy contacts and
scoring defaults were not changed.

Fresh read-only review found a parsing/eligibility defect: the default parser can
drop HETATM amino acids before standard-peptide validation. A real-PDB regression
test demonstrated incorrect acceptance; the new quantity now opts into the
existing `keep_modified=True` path. The test then passed, and the reviewer checked
the fix and latest scoped code with no remaining material findings.

- Focused geometry/contact/export scope: **49 passed**, 1.63 s pytest time,
  2.20 s including process startup, 120 s timeout. Exact selected files are the
  nine contact/geometry/export files in the implementation plan's verification.
- Existing chemical typing suite: **41 passed**, 13.14 s, 120 s timeout.
- Ruff `check src tests` and `git diff --check`: passed.
- Sphinx HTML with `SPHINXOPTS=-W`: passed after correcting a new heading length
  and marking the internal implementation plan orphaned.
- The broad unit command `pytest tests/unit -m 'not slow' -k 'not refitting and
  not reproduc' -q -x --tb=short -p no:cacheprovider` **did not complete**: its
  process reached the 300 s timeout after progress beyond 56%. Collection placed
  that section in historical Potts/ensemble tests. No completed full-suite pass
  is claimed. A process check found no leftover pytest/annotation subprocess.
  The subsequent targeted 49-test scope and 41-test chemistry scope are the
  completed evidence for this change, not substitutes labelled as the full gate.

The repository venv initially imported a previously installed wheel. It was
restored to the documented editable development install using the offline uv
cache before source verification. This summary is source-development evidence;
future benchmark comparisons still require an immutable noneditable wheel.

## Native population and counts

Raw input: the existing `paths.native_dir()` Native2026 inventory, 374 compressed
structure files totaling 36557448 bytes. The four entries in the bundled withdrawn
table were excluded before annotation, leaving **370 inputs**. The batched CLI
accounted for every input: **357 successful complexes**, **13 eligibility errors**.
Successful complexes comprise **263 class I** and **94 class II** complexes.

| Class | Interface | Successful complexes in class | Complexes with contacts | Atom pairs | Residue pairs |
|---|---|---:|---:|---:|---:|
| I | CDR:peptide | 263 | 263 | 48827 | 5641 |
| I | CDR:MHC helix | 263 | 263 | 74321 | 9565 |
| I | peptide:MHC helix | 263 | 263 | 123452 | 11681 |
| II | CDR:peptide | 94 | 92 | 15022 | 1885 |
| II | CDR:MHC helix | 94 | 94 | 27624 | 3532 |
| II | peptide:MHC helix | 94 | 94 | 50675 | 5531 |

Total: **339921 unique physical atom pairs** and **37835 unique physical residue
pairs** across the three interfaces in 357 complexes. These contact counts are
not independent sample sizes. Two successful class-II complexes have zero
CDR:peptide contacts and remain in the class denominator. The machine-readable
six-row table is `NATIVE_CONTACT_COUNTS.tsv`; missing-omega counts there are atom
rows, not independent residues. Class-II native observations establish no
modelled class-II classifier performance.

The 13 errors comprise ten missing required-chain cases, one non-alpha-beta case,
and two strict standard-peptide exclusions. In 5JHD and 5JZI the preserved peptide
includes `AMN` at author residue 0, represented as `X`. This explains the difference
from the older 359-structure recipe under historical parsing. No historical
potential was regenerated, no established score was changed, and no manuscript
result was rewritten to conceal this population difference.

## Raw-coordinate verification and provenance

The verifier reread raw coordinates without repeating annotation or scoring.
Every endpoint of every exported pair matched physical chain/author residue/
insertion/atom/altloc identity. All **339921** raw-coordinate distances matched
the export with maximum absolute error **0 Å**. Pair identities were unique,
all distances were within the inclusive 5 Å cutoff, unordered chemical keys were
ordered, successful accounting summed to the row count, and separate atom/residue
summary counts matched the physical observations exactly.

Requested deliverables live under `data/contact-observations-20261018-a/`:
`contacts.parquet`, its summary/accounting/residue/annotation/provenance sidecars,
raw-input manifest, `SOURCES.json`, `RUN.json` and `VERIFICATION.json`.
Total directory size measured: **35482557 bytes**. The observation Parquet is
**13323022 bytes**, SHA256
`8947bdd7e2947f7fdafe9c4fa1b8b97b0a537b8bf99c42dd1a7517f602130539`.
`SOURCES.json` pins every raw structure, the withdrawn table and producing source
files by hash. Provenance separately pins the actual searched annotation indexes.
These are requested final observations and source records, not a resume cache.

Exact production command (from the repository root):

```sh
.venv/bin/tcren export -s data/contact-observations-20261018-a/native.txt --quantity contact-atoms --interface all --cutoff 5 -t 1 -o data/contact-observations-20261018-a/contacts.parquet
```

The caller used a new output path, a 300 s subprocess timeout and measured elapsed
time **62.8611 s** and peak child RSS **1875066880 bytes** (about **1.75 GiB**).
No new scientific output was reused to build the cohort; every selected raw
structure was parsed/annotated for this authorized run. Source inputs and old
results were preserved. No Stratos, MD or publication/release action occurred.

## Remaining scope

Retrieve and inspect the external quantitative prior files documented in
`CONTACT_PRIOR_SOURCES.md`; construct and check low-order native/reference
distributions, including proper support and endpoint symmetry. Full hydrogen-bond
angles, line-of-sight and chemically justified unnaturalness assignments remain
to be implemented. Native rare-contact counts are not physical-invalidity labels.
Freeze the assignment and declared penalty menu before matched benchmark
evaluation. All original energy/ipTM baselines, strata and regressions remain
requirements of that later stage.
