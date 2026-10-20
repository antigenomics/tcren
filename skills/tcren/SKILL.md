---
name: tcren
description: tcren — TCR-pMHC contact potential (TCRen) pipeline; conventions and public API
---

# tcren Skills Guide

## Classic reference identities

The `karnaukhov2022` CLI alias / `potential.tcren()` retains the bundled May-2026
updated classic table. The original December-2022 root table, matching saved
Figure 3 and the March-2026 root, is a distinct frozen resource; use an explicit
CSV when comparing that baseline. Native populations and redundancy/debiasing
recipes must be named separately, along with full/LOO/LCO status. See `SOURCES`;
no coefficient/default changes follow from this provenance clarification.

## Symmetric chemical core, directed composition

`contact-chemical-core-fit --native contacts.parquet --external curated.parquet
-o chemical-core.json -t 0` freezes `chemical-core-count-1`: a potential-specific
19-type map, unordered physical counts and symmetric enrichment W against pooled
role-specific endpoint products. Native class-I CDR:peptide and curated antibody
peptide contacts use actual integer counts at 5 Å; no assay/source weights.
Conditional role marginals define `Z_ab=sum qB(c|a)qT(d|b)W_cd` and
`Delta=-ln W_cd+ln Z_ab`. General contact classification stays unchanged.
`contact-penalties ... --chemical-core-resource chemical-core.json` adds one
physical-atom-pair mean per known CDR edge, unit coefficient, as `chemical_core`.
Unseen W=1 yields ln Z, not final zero. Unknown mappings/zero conditional support
abstain with final zero, distinct counts and unchanged mean denominator.
Resource/input hashes and original scores remain explicit; `contact-compare`
retains the separate arm. No new default, physical free energy or accuracy claim.
See `docs/contact-penalties.rst` for equations, provenance and support semantics.

## Within-residue atom-pattern layer

`contact-atom-pattern-fit --native contacts.parquet --external antibody.parquet
--equivalence-resource prior.json -t 0 -o pattern.json` freezes edge-balanced
atom marks: each physical residue edge contributes one unit divided among its
atom pairs, before saved chemical-equivalence aggregation. Source raw counts
and fractional masses remain explicit. D=-ln[P(i,j|a,b)/(qB(i|a)qT(j|b))]
uses the resulting directed endpoint marginals. Native class-I CDR:peptide and
antibody peptide references only; no assay/source weights or distance regression.
`contact-penalties --atom-pattern-resource pattern.json` adds one mean D per
known CDR edge to original TCRen, fixed coefficient one. Unsupported pairs abstain
with zero and remain in the denominator. Separate export from atom-count resource;
`contact-compare` retains `atom_pattern` and all original arms. This categorical
within-contact measure is not full atom occupancy or a joint graph likelihood.
See `ATOM_PATTERN_SPEC.md` and contact-penalties docs for support/provenance.

## Directed residue+atom count layer

`contact-atom-count-fit --native contacts.parquet --external curated_atom_contacts.parquet
--equivalence-resource prior.json -o atom-count.json` freezes directed conditional
atom enrichment from existing native class-I CDR:peptide and curated antibody
peptide contacts, pooling actual integer counts. Saved chemical equivalences
apply to both references and queries without collapsing physical observations.
`D=-ln[(n_xy/N_ab)/(qB(i|a)qT(j|b))]` uses contact-endpoint marginals, not exposure.
`contact-penalties ... --atom-count-resource atom-count.json` adds one mean D per
known CDR residue edge to each original score; unit coefficient, no assay fit.
Unknown cells contribute neutral zero with explicit counts and remain in the
average denominator. Framework/missing original entries retain their policy;
missing CDR observations fail. `contact-compare` reports `atom_count` alongside
all original arms. Class I only; no angular/distance regression, surface-area
measure, physical-validity probability or new default. See contact-penalties docs.

## Experimental frozen contact penalties (18 October 2026)

`contact-prior-geometry --source DEPOSIT --exclude native_ids.txt -o NEW_DIR`
recovers external raw frames only through unique exact float32 chemical/coordinate
correspondence; numbering is not assumed equal. It handles one raw mmCIF at a time.
`contact-prior-fit --external external_geometry.tsv.gz --native native.parquet
-o RESOURCE.json` freezes approximate normal marginals of logit(d/5 Å) and
atanh(cos(angle)), external/native sufficient statistics and declared chemistry
backoff. The external actual counts supply prior strength; there is no label input.
Angular torsion sign is discarded; missing/boundary geometry remains unavailable.

`contact-penalties -s manifest.txt --resource RESOURCE.json --chunk 512 -t 0 -o NEW_DIR`
exports original TCRen2 and bundled classic TCRen energies plus fixed, matrix-scaled, literal
flip and favourable-only flip contributions. It keeps the legacy all-TCR peptide
baseline, adding one increment per physical CDR:peptide residue edge. Framework
contacts remain explicitly unassessed. Worst supported atom/marginal predictive
tail defines the shared weight; no score-time fitting, probability-of-naturalness
claim or independent-marginal product. Missing potential entries retain legacy
omission with explicit counts/null contributions. Batched annotation is bounded
by `--chunk` (default512); parent-only batched annotation precedes coarse
contiguous spawn-worker scoring. `-t 0` (default) uses all allocated CPUs; `1`
is explicit serial and positive N is a bounded allocation. Native kernels stay
serial inside workers; pool failure propagates. No computed-part cache.

`contact-compare --labels labels.tsv --scores NEW_DIR/scores.tsv -o COMPARISON`
uses exactly the same finite IDs for all arms/original/ipTM within every potential,
cohort, MHC class and template stratum. Existing ROC/AP/precision-at-recall metrics,
membership and all gains/losses are exported. See `docs/contact-penalties.rst`.
The resource is an estimated experimental dependency, not a new final-paper model.

``contact-geometry-audit`` exports fresh bounded all-atom chemistry plus exact
worst-tail lineage with the frozen resource; it does not recalculate energies.
``contact-mechanism --run ACCEPTED_RUN --geometry NEW_DIAGNOSTICS -o NEW_DIR``
verifies unchanged row membership, weights, signed sums and accepted AUROC deltas
via pair-concordance accounting. Outputs separate chemical/proximity counts from
continuous Σ(1−u)/Σu geometry equivalents; motif/matrix selection and changed-rank
examples retain cohort/class/template strata. No binary geometry cutoff or new
score is selected. See the mechanism section in `docs/contact-penalties.rst`.

## Multiplicity-calibrated weighted contact scores (18 October 2026)

`contact-weight-calibrate --external external_geometry.tsv.gz --native native.parquet
--resource marginal-resource.json -o calibration.json` freezes conditional empirical
ranks of the existing minimum predictive tail. The opportunity group counts supported
atom-distance rows and distinct `(key,x)` omega/theta/phi terms. Native CDR:peptide
strata take precedence; external antibody:antigen is used only for native-empty
strata. Weight `w=(1+#reference_tail<=query_tail)/(N+1)` retains ties inclusively.
Export native/external source, edge/structure counts and missing support. Unseen
strata/missing geometry retain neutral weight one explicitly. No labels or tuned
thresholds/coefficient enter calibration. The marginal resource hash is bound.

`contact-penalties ... --calibration calibration.json` additionally exports weighted
entries `w*e`, with unchanged framework/missing-entry policies and original matrices.
`contact-compare` then compares original/weighted/ipTM on identical finite IDs.
`contact-weight-summary --run NEW --accepted OLD --mechanism PREVIOUS -o AUDIT`
verifies exact original scores/graphs/raw membership and exports signed losses,
CDR-only equivalents, support, motifs and exact pair-concordance accounting.
Reference ranks do not establish contact-validity probabilities, independent
calibration or binding free energies. See `docs/contact-penalties.rst`.

## Directed reference-partition correction

contact-normalizer-fit --native PARQUET --resource V2 --mhc-class I --out JSON
freezes logZ=log(mean_native(exp(-u))) in exact native CDR→peptide opportunity
strata (n_distance,n_omega,n_theta,n_phi). contact-penalties --normalizer JSON
adds g=u+logZ once per edge; contact-compare retains all original arms.
Empirical mean exp(-g)=1 and constant penalties cancel. Source version/hash/class
are checked. Unseen groups have zero increment with a flag; missing geometry
stays unavailable; framework and missing-table policies are unchanged.
Partitions are estimated reference dependencies with no assay-fitted weights.
The relative factor may have either sign and does not supply binding free
energy or SASA/pocket normalization.

## Opportunity-aware contact-class audit

`export --quantity atom-inventory --interface all -o atoms.parquet` exports all
selected observed heavy atoms in full TRA/TRB/peptide/MHC/B2M components, including
zero-contact atoms, with the existing contact conformer/source-index policy.
No potential, SASA, bonds/protonation inference or energy is used. Inventory
membership includes context/framework/Gly-extra atoms; contact_role/eligible
retain the separate historical scope. Source index is parsed-heavy-atom ordinal,
not PDB serial; peptide.pos is observed one-based ordinal.

`contact-class-audit --manifest labels.tsv -o NEW_DIR -t 0 --chunk 128` performs
fresh batched annotation/coarse-worker structural export, then opportunity-aware
statistics in bounded feature blocks, retaining all models/clusters per feature.
Final Parquet paths can be partitioned datasets read by `polars.read_parquet`.
Required manifest: id,label,receptor_cluster,dataset,cohort,epitope,
mhc,mhc_class,template_coverage,source,sha256. All selected TCR regions are separate
binder contexts; targets are peptide and annotated MHC helices. Physical atom and
residue contact counts use fixed5A. Sparse opportunities multiply full matching
endpoint availabilities, retaining C0/Opositive and explicit O0 missingness.
Site views stay case/allele-specific; overall chemistry contrasts weight cases
equally and retain MHC class/coverage and separate TCRvdb epitopes. Exact support,
prevalence, density, unsmoothed ratios and receptor-clustered CR1/t uncertainty
are exported; unreplicated/degenerate inference stays null. All failures retain
accounting and block statistics. No score, fitted weights, significance filtering,
FDR/permutations, cached feature inputs or default promotion. See
`docs/contact-penalties.rst` for full definitions and observational limits.

## Contact observations (18 October 2026)

`tcren export -s native.txt --quantity contact-atoms --interface all -t 1 -o atoms.parquet`
uses the existing batched annotation/search and one neighbour search per structure.
No potential required. `all` selects CDR1/2/3:peptide, CDR1/2/3:MHC helices and
peptide:MHC helices; the existing three interface names select individual arms.
MHC groove-floor/B2M contacts are excluded. Missing MHC helix annotations fail
explicitly. Existing quantities and legacy conformer handling are unchanged.

The new quantity selects the lexicographically first nonblank altloc per residue,
plus blank atoms, with selected-label atoms winning name ties; this is not an
occupancy selection or proof of inter-residue conformer compatibility. It retains
altloc/source atom indices and omitted-atom counts. Glycine uses observed backbone
atoms only. `chem.a`/`chem.b` canonicalize whole residue+atom tuples; physical sides
retain chain/author residue/insertion, coordinates and biological region labels.
Observed N/CA/CB define trRosetta omega/theta/phi in radians; glycine/missing CB and
degenerate geometry remain null. Current v2 chemistry flags remain heuristics;
stacking uses the existing ring geometry. No angular hydrogen-bond assignment,
plausibility probability, penalty or energy is calculated.

The new quantity writes compressed Parquet for a `.parquet` output suffix,
otherwise TSV. `atoms.parquet.summary.tsv` separates exchangeable atom chemistry and directed binder-to-target residue frequencies by
class, interface and receptor region. Residue distances are closest retained atom
distances and type flags are unions over atom pairs, not closest-pair labels.
Overlapping type counts do not sum to totals. The observation API is
`contacts.observations.contact_observations`; table-only reduction is
`summarize_contact_observations`. Accounting includes zero-contact structures;
summary rows describe observed groups, not missing groups or independent samples.
Summary hash, annotation/index lineage and conformer policy remain in provenance.

## MD diagnostics (issue 26)

`tcren md inspect --trajectory protein-continuous.pdb --reference first.pdb
--potential potential.xvg --components components.json --mhc-class I
--pbc-whole --docking -t 1 -o NEW_DIR` exports the exact existing
`validation.md.inspect_md_trajectory` readouts: partner-relative backbone RMSD
(Å/ps), contact counts/occupancy, optional docking, first/final/minimum-potential
saved frames (whole-system kJ/mol) and source/output provenance. Component JSON
maps disjoint names to chain-ID lists; `all` is reserved. Class II uses a separate
export with both groove chains in MHC; class is an explicit caller declaration.
Exact correspondence includes author residue/insertion and altloc identities.
PBC preparation must be done beforehand; the flag acknowledges it, not proof.
`scripts/inspect_gromacs_trajectory.py` retains the GROMACS whole/compact-anchor/
continuous-nojump export workflow. No active Stratos job is restarted by inspection.

`tcren md quality -s model.pdb --mhc-class II --cutoff 2 -o NEW_DIR` exports
`coordinate_diagnostics` and `close_atom_pairs` with physical atom/altloc identities,
possible-disulfide flags and hashes. Covalent classification remains unknown;
adjacency exclusions are heuristic, not a MolProbity clashscore.

`tcren md compare-forces --reference-energy REF.xvg --reference-forces REF_F.xvg
--energy E.xvg --forces F.xvg --atoms N --topology COMMON.tpr
--coordinates SAME.gro --same-system --mhc-class I -o NEW_DIR` exposes
`compare_same_frame_energy_forces` (kJ/mol and kJ/mol/nm). One identical timestamp,
System-group order and shared input provenance are required; the acknowledgement
does not independently prove XVG origin. Preserve the rerun commands/MDPs.

The trajectory API/export enforces `threads=1|0|N` within allocation, serial native
kernels and a separate batched reference-annotation phase. Quality/force diagnostics
are serial. No fitting, inherited predictor confidence or binding-energy estimate.
Long campaign endpoints and missing original MD inputs remain separate gates.
See `docs/md-validation.rst` for tables and limitations.

## Structural alanine reference (manuscript convention)

Use `scoring --delta --reference-mode structural` for peptide references and
`export --quantity receptor-reference --reference-mode structural` for receptor
references. These reuse the existing coordinate substitution, preserve backbone/Cβ,
truncate other side-chain atoms, construct Cβ for glycine, and rebuild contacts at
the selected cutoff. The default `virtual` preserves the older fixed-contact-map
identity substitution for historical reproduction; it is not the manuscript
calculation. Rows and provenance identify the mode. Compare both on matched rows
and retain regressions; never repair a scientific mismatch by rewriting the intended
Methods to describe an accidental CLI path.

Performance: build the structural peptide reference map once per structure and
reuse it across interfaces. Profile contact construction separately from table
scoring and coordinate/SASA descriptors. CPU limits are owned by the package, not user shell exports. The repaired surface/scoring/feature/export and annotation paths
use `threads=1` (serial), `0` (all allocated) or a positive count within allocation.
Annotation and coarse spawn-worker phases run sequentially, each worker using
one numerical thread. Import tcren before initializing Polars; an existing
multithreaded Polars pool cannot be resized and is rejected with an actionable
fresh-process error. Runtime BLAS/OpenMP limits and environment settings are restored.
`surface -s MANIFEST -t N` uses one TCR query batch per species and one MHC batch,
keeps per-input failures, validates positive grid dimensions before annotation,
and exits nonzero for empty/all-failed output. Mixed runs retain explicit errors.
Polars async I/O is limited to one thread and blocking I/O to four auxiliary
workers per process; a blocking pool of one can deadlock and is not used.

## pMHC feature–performance associations

`tcren export -s structures.txt --quantity pmhc -o pmhc.tsv` calls
`topology.pmhc.pmhc_features`: surface relief/protrusion/chemistry plus peptide
and side-chain SASA (Å²) and exposed-area-weighted Kyte–Doolittle hydropathy.
The receptor is removed before defining the frame; supplied coordinates remain
unchanged. No potential required. Class II retains both MHC chains; no guessed
core/register. Accounting and searched-index lineage sidecars apply.

`tcren associate --results recorded.csv --features pmhc_by_id.tsv --scores scores.json
--feature-columns relief,frac_above_ridge,sasa_hydropathy -o associations` compares
explicit fixed-direction scores against ipTM on pair-specific finite rows.
Tables require unique `id`; results additionally require dataset/task/cohort/
mhc_class (I/II)/template_coverage, binary y and iptm. Unknown coverage is explicit.
`associations.associate` returns metrics, feature summaries, descriptive Spearman
associations and row membership. Classes/template strata stay separate. No fitted
score, label-selected direction, p-value or threshold sweep. Missing features do
not drop score rows. New output directories only; all-unavailable metrics exit
nonzero. See `docs/associations.rst` for equations, schemas and interpretation.

Joint scalar context: repeat ``--feature-product relief:sasa_hydropathy`` and
``--feature-product frac_above_ridge:sasa_hydropathy`` to append raw signed
products (``LEFT__x__RIGHT``), automatically included in associations. The API
``associations.add_feature_products(features, pairs)`` uses vectorized Polars;
no labels, centering, fitted weights or thresholds. CLI ``joint_features.tsv``
and provenance retain exact per-model products/source pairs. Nonfinite/overflow
products remain missing without changing score/ipTM membership. Summaries use
the median of products, not the product of medians. These scalar descriptors do
not establish spatial co-localization, independent interaction or a new score's
performance gain. Compare against component associations on the same cohorts.

## Final-model boundary and exact exports (2026-10-15)

Fixed empirical contact potentials are permitted with explicit derivation and table
identity. Q/T/S reference moments/correlation weights, Gaussian/GLM, fitted Potts and
transition readouts are historical APIs, excluded from final models. Older fit-free
labels below mean no fitting at scoring time and do not establish final eligibility.
Benchmark-local available-pair contact GLMs and Gaussian diagnostic fits are not
equivalent to library Potts/assess commands by name. Preserve their historical
results separately; final producers read fixed-potential and individual-descriptor
CLI tables. `audit/2026-10-16/FINAL_PRODUCER_DISPOSITION.md` records this disposition.

`tcren export -s INPUT --potential tcren2 --quantity receptor-reference -o OUT.tsv`
exports exact virtual receptor poly-Ala differences on unchanged contacts for all,
CDR, CDR1/2 and CDR3, TP/TM and their totals. `--quantity position|residue|contacts|placement`
uses the exact existing decomposition/placement APIs, not leave-one-out explain.
Inputs include files, directories, manifests/globs and supported archives. Full stems
and physical residue identities are retained; duplicates fail. Accounting, residue,
annotation and potential/provenance sidecars retain failures and estimated dependencies.
All-failed/empty commands exit nonzero. Benchmark consumers must check every input key.
Molecular panels must read these exported rows, never recompute the quantity in PyMOL
or substitute `explain`. Join physical chain/author residue/insertion identities,
preserve exact region membership, and retain potential, raw-coordinate, scene, camera
and transformed-geometry hashes. Verify actual atom colours against the export.
Keep class II, epitopes and template strata separate; inspect every final PNG and
its composed page. The issue-34 acceptance record is
`audit/2026-10-16/MOLECULAR_CLI_LINEAGE.md`.

Annotation requires ARDA >=2.36.2. `Structure.annotation_provenance` captures the
resolved searched database prefix and sequence/header/index SHA256 values around
each batch. Scoring, feature, recognition and export sidecars retain these under
`searched_indexes`; mixed generations or mid-search changes raise. An empty record
means no searched identity was supplied, never version-only proof of reproducibility.

`tcren inspect-models manifest.tsv output/` joins native predictor confidence with
batch annotation, contacts and descriptors. Historical fitted assessment is disabled
by default; `--historical-assessment` explicitly opts in. AF3, TCRmodel2, Boltz2,
tFold-TCR and AF2-multimer have format-specific support; Protenix/RF3 scalar summaries
retain native names, with unverified token mapping explicitly unavailable. No ranking
score is silently relabelled ipTM. See `docs/prediction.rst` and `docs/fit-free.rst`.

`score_table` has per-readout availability flags and JSON `missing_descriptors`.
`recognition_table(..., annotate=False)` reuses supplied annotated Structures.
Parallel descriptor calls from stdin/notebooks use an importable subprocess entrypoint
for the existing bounded spawn pool; no notebook replay, fork or serial fallback.
Whole-set annotation remains in the caller, before workers start. Ordinary scripts
still need a `if __name__ == "__main__":` guard. Import tcren before initializing
Polars, including in notebooks. Temporary transport is private and deleted after each call.

`tcren` reproduces and extends the TCRen contact-energy potential (Nat Comput Sci 2022)
on a pure-Python pipeline (structure parsing → contacts → TCR/MHC annotation → potential
derivation → epitope-ranking benchmarks). Annotation uses the `arda` package
(mmseqs2-backed), a runtime dependency published to PyPI as `arda-mapper` (imports as `arda`,
`>=2.36.2`) — no separate checkout and **no `ARDA_HOME`** (arda auto-fetches its reference and a
static mmseqs binary on first use). Repo-local `.venv` via uv, no conda (`bash setup.sh`).

## The flat middle became four named sub-packages (2.29.0) — old imports still work

`docking` (renamed from `orient`), `topology`, `energetics` and `mechanics`. **Every old top-level
name is a transparent re-export and keeps importing**, so nothing breaks; prefer the canonical home
when writing new code, since that is where the docstrings and the layering audit
(`scripts/audit_architecture.py --check`) live.

| deprecated location | canonical home |
|---|---|
| `tcren.orient` | `tcren.docking` (`orient.docking` → `docking.angles`) |
| `tcren.footprint` | `tcren.topology.footprint` |
| `tcren.surface` | `tcren.topology.surface` |
| `tcren.pose` | `tcren.topology.pose` |
| `tcren.interface_graph` | `tcren.topology.graph` |
| `tcren.ddg` | `tcren.energetics.mutation` |
| `tcren.rotamers` | `tcren.energetics.rotamers` |
| `tcren.scoring` | `tcren.energetics.scoring` |
| `tcren.dynamics` | `tcren.mechanics.dynamics` |
| `tcren.stability` | `tcren.mechanics.stability` |

`energetics.ddg` is `energetics.mutation` because as a submodule it shadowed the `ddg` *function*
of the same name on the package. Section headings below still use several of the deprecated names;
they resolve, and this table is the map.

## Current audit boundary (2026-10-14)

The geometry prior, coordinate-chain and cross-product experiment, its bundled geometry
resources and CLI commands `potts geom`, `chain`, `cross`, `predict` are archived under
`attic/2026-10-audit`. The retained contact model still uses supplied coordinates;
`potts contacts`, `agreement`, fitting and self-consistent regeneration remain supported.
Historical release notes are in `CHANGELOG.md`; publication gates are in `ROADMAP.md`.

`features --chunk N` limits structures held in memory and always recomputes from raw
inputs. There is no resume mode or computed part cache. Empty/all-failed feature runs
exit nonzero and retain their diagnostic table.

`structure.Atom.altloc` retains source alternate-location labels without changing atom
order or contact inclusion. `validation.structures.close_atom_pairs` exports both labels,
`alternate_location_involved`, `possible_disulfide`, and `bond_classification="unknown"`.
Adjacent-residue exclusions are heuristic; alternate labels do not establish cross-residue
incompatibility. PDB/mmCIF writers keep the existing one-conformer selection and its label.

Native array boundaries reject invalid dimensions, indices and scalar domains before access.
Frozen model loaders reject stale catalogues. A model without sufficient valid confidence
training returns an unavailable (NaN) confidence residual. Feature/manifest joins validate
unique non-null keys and binary labels; fit metadata reports unmatched and incomplete counts.

The 2026-10-14 definition migration preserves frozen coefficient arrays and numerical fit
metadata; `audit/2026-10-14/frozen-metadata-migration.json` records their hashes and both
catalogue digests. This is not a refit or proof that corrected class-II/shear features reproduce
historical outputs. Mechanics defaults have reduced stiffness Å^-2, force Å^-1 and dimensionless
work; no SI energy calibration is supplied. `fp_chi_*` is β₀−β₁, equal to full Euler characteristic
only when higher Betti numbers vanish. `exp_lost` is a linearized nearest-pair loss sum.

Potts `available_pairs` exports `chain.id.rec`, `chain.id.par` and heavy-distance `probe`.
Kernel and marginal joins retain physical chains. Legacy peptide-only tables can infer one partner;
legacy MHC tables must be regenerated. `threaded_ensemble` requires peptide-only sites and exact
candidate `par.len`, and carries every Gibbs chain's final state into the next candidate.
`n_parameters` counts stored coefficients, including orientation, not identifiable degrees of
freedom. The imputation fitter combines Gibbs imputation with conditional pseudolikelihood;
it has no likelihood-EM monotonicity theorem. Coupled-marginal `p_any` remains an explicitly
labelled independence approximation. Sampling trajectories after these repairs require fresh
reproduction before scientific promotion.

`scripts/check_fresh_holdout.py --features <fresh harmonization parquet> --out <report.json>`
requires every shipped manifest ID and compares a temporary full refit against the frozen arrays
and numerical metadata under the existing exact full-scale contract. It never installs the refit.
The driver owns raw-input freshness; invoking this checker on an old table is not reproduction.

## Batch annotation — never loop (mmseqs2 is the parallel layer)

**All structure annotation (TCR chain typing AND MHC allele mapping) must gather every
sequence first, make ONE batched mmseqs2 call, then map the output back for downstream
per-structure analysis.** mmseqs2 parallelises internally across threads — that is the
parallel layer; Python orchestration is a single call.

- Each per-structure annotate call pays a fixed ~825ms mmseqs2 process+index-load cost;
  a batch of 300 sequences costs the same ~930ms total.
- A `ProcessPoolExecutor(fork)` over structures **deadlocks** (fork after mmseqs2/BLAS
  spawn threads). A `ThreadPoolExecutor` runs but still pays the fixed cost N times.
- `annotation.batch.annotate_batch` (also exposed by paper helpers) does TCR
  annotation for a whole dataset in 2 arda calls (human + mouse). MHC annotation uses the same pattern: `mhc.annotate_mhc_batch(structures)`
  — ONE mmseqs search over every candidate MHC chain, sliced back per structure.

Reference: `arda.annotate_sequences([(id, seq), ...])` — one call, threads internally.

## Parallel phases — batch searches, then structural work

- **Annotation (TCR + MHC) is never Python-threaded and never per-structure.** It is one
  batched call per reference with the selected CPU budget passed explicitly to MMseqs.
  No Python pool wraps the search. The descriptor phase then reuses the annotated
  structures in bounded spawn workers, with one numerical compute thread each.
  Notebook/stdin calls use the importable entrypoint described above.
- **Outer parallel work is restricted to mmseqs-free stages:** structural
  alignment (Kabsch/SVD superposition), peptide mutation, relaxation, and rendering — i.e.
  pymol / Rosetta / FlexPepDock and figure generation. `orient.run_folder(threads=…)` threads
  the parse and the align+write stages (default `os.cpu_count()`); annotation between them is
  the single batched pass. `tcren orient -t N`. **`superimpose` is the same**: `run_superimpose`
  batch-annotates all inputs, then threads the mmseqs-free ensemble alignment + write
  (`superimpose(..., annotate=False)` on the pre-annotated structures). `tcren superimpose -t N`.
  `-s` accepts file / dir / .tar.gz / glob; `-o` is a directory, or a single structure file
  (one input) whose extension must match `--mmCIF`/`--compress` (validated by `_output_target`).

## Two orientation commands — `superimpose` vs `orient`

- **`tcren superimpose` / `tcren.orient.superimpose(s, db_dir=…)`** — bring a NEW structure into
  the canonical frame against a canonical *database* (default `data/Canonical2026`). It detects
  the input's MHC class + species, selects every DB member of that class+species (from the DB's
  `orient_metadata.json` — bundled in the package for the shipped `Canonical2026`, which
  `fetch-data` populates with structures only), superposes the query groove Cα onto each, and
  **averages** the rigid transforms (chordal/SVD mean rotation + mean translation) into one
  consensus placement. The matching DB subset is batch-annotated once and cached per process.
- **`tcren orient` / `tcren.orient.run_folder(...)`** — BUILD a canonical DB from native
  complexes using the per-class derived frame (how `Canonical2026` is produced). Not for orienting
  a single new structure — use `superimpose` for that. Writes `<out>/orient_metadata.json` (the
  format `superimpose` reads); `--metadata foo.csv` writes CSV instead.
- **HF upload is NOT a user command.** `--push-to-hub` was removed; maintainers run
  `scripts/push_canonical_to_hub.py` instead.

## Structure scoring — `tcren.run_pipeline` / `tcren scoring`

- **It is called `tcren scoring`, not `pipeline`** (renamed 2026-07; the old name errors with a
  pointer). It only scores: canonicalisation, region mapping and the Cα/contact/atom-distance
  matrices are `tcren annotate` / `superimpose` / `contacts`.
- `run_pipeline(structure, superimpose=True, db_dir=…)` → `PipelineResult`: import → annotate
  (alleles + chains + MHC groove) → superimpose onto the canonical DB (canonical Cα) → resmarkup
  + 5 Å contacts → per-interface energies. Scores: **TCRen** for TCR↔peptide, **MJ** for TCR↔MHC
  and peptide↔MHC, plus `total` (sum of the residue-pair potential over each interface's contacts).
- **Φ and ΔΦ with the TP / TM / PM breakdown come from this one command.** Columns `Phi_tcr_pep`,
  `Phi_tcr_mhc`, `Phi_pep_mhc`, `Phi_total`; `run_pipeline(…, reference_aa="A")` /
  `tcren scoring --delta` adds `dPhi_*` and `dPhi_total` (each interface's
  `tcren.ddg.reference_delta`; `dPhi_tcr_mhc` ≡ 0, the peptide is not in that interface). Use ΔΦ,
  not Φ, when each candidate carries its **own** generated pose. **Column names are shared with
  `tcren recognize`** — one vocabulary — but the **keys differ**: `scoring` emits `pdb.id`,
  `recognize`/`features` emit `complex.id`, so rename one before joining (`scoring --geometry` does
  that rename internally).
- **`F_*` → `Phi_*`, `dF_*` → `dPhi_*`, `F_TERMS` → `PHI_TERMS`, `f_score` → `phi_score` at
  2.26.0** — in the code and the emitted columns, not only in prose, and with no alias. **`d` in
  `dPhi` is the reference difference** ΔΦ = Φ(sequence) − Φ(reference), never a derivative, and the
  package carries **no `dd` quantity** except `ddG`, the change in binding free energy on mutation.
  (Manuscript notation: Φ for the energy, φ for a potential matrix entry; `cohort.phi_score` is the
  binder-direction channel `−Φ_TCR:pep`.)
- `tcren scoring --geometry` appends the interface descriptors + `Q` by calling
  `recognition_table` + `cohort.q_score` — it does **not** reimplement them.
- `-s` accepts a file, directory, `.tar.gz`, quoted glob, `.txt`/`.list` manifest (one path per
  line, `#` comments), comma-separated list, or a repeated flag. `structure_paths` handles globs
  and manifests; `resolve_sources` does the comma/repeat splitting.

## Compiled extensions — `_align`, `_refine`, `_relax`, `_fold`, `_geom` (pybind11 / scikit-build-core)

- FIVE C++ exts, all in `CMakeLists.txt` (`pybind11_add_module` each, `install(TARGETS _align
  _refine _fold _geom _relax …)`): `src/_align/align.cpp` (MHC pseudoseq fitting alignment),
  `src/_refine/refine.cpp` (potential-guided peptide refinement), `src/_relax/relax.cpp` (DOPE
  interface energy for `tcren energy`/ΔΔG), `src/_fold/fold.cpp` (CCD loop closure) and
  `src/_geom/geom.cpp` (the clash and contact-fragility kernels: `interface_clashes` behind
  `tcren.clashes.interface_clashes`, `contact_stability` behind
  `tcren.mechanics.stability.contact_stability`). Adding a sixth = same pattern; a new *function*
  in an existing module (as clashes/stability were added to `_geom`) needs no `CMakeLists.txt` change.
- **Three of `_geom`'s five entry points have no Python caller left.** `shrake_rupley`,
  `interface_hbonds` and `contact_descriptors` built the interface-geometry table for the
  `tcren binder` command, which 2.26.0 removed with its frozen model; they are exercised only by
  `tests/unit/test_geom.py`. Do not reach for them expecting a supported path.
- The MHC-pseudosequence fitting-alignment hot path is a C++ ext (`src/_align/align.cpp`,
  `CMakeLists.txt`). Build backend is `scikit-build-core` (not hatchling); `uv pip install -e .`
  builds it once at install. Funcs: `fitting_score`, `best_hit` (GIL released over
  candidates), `align` (traceback). Scoring matches Bio.Align's fitting config EXACTLY (BLOSUM62,
  placed-gap open -11/extend -1, free target + end gaps), so `tcren.mhc.pseudo` falls back to
  Biopython transparently when the ext is absent. ~40 ms vs Bio 59 ms vs pure-Python 15 s for 4k
  candidates (a modest 1.5x — Bio's aligner is already C). `editable.rebuild = false` in
  pyproject: the ext is built once at `pip install -e .` (do NOT rebuild on import — that needs
  cmake on PATH at import time and breaks pytest/CI). CI (uv, `astral-sh/setup-uv`) installs
  `uv pip install --system -e . --no-deps` + explicit runtime deps (so arda-backed tests skip) and
  `uv pip install --system cmake ninja` to build the ext.

## Peptide substitution + refinement — `tcren.refine` (`tcren refine`)

- `substitute_peptide(structure, new_peptide)` — backbone-preserving identity swap on the peptide
  chain (keep N/Cα/C/O+Cβ, drop side-chain beyond Cβ); pure data-model, no atoms moved. Region
  markup is carried over onto the new residues — without it the contact map has null `pos.from/to`
  and every scorer fails. `score_peptides`
  is *virtual* (matrix lookup over the fixed contact map) — substitution is only needed to then refine.
- `refine_peptide(structure, restraint_w=0.5, …)` → `(structure, energy)`: knowledge-based rigid-body
  **Metropolis MC** of the peptide via the `_refine` C++ kernel. Energy = **DOPE** atom-level
  distance-dependent statistical potential (Shen & Sali 2006) over all peptide↔partner heavy-atom
  pairs (chain-agnostic, atom-class × distance lookup with linear interp) + **harmonic restraint to
  the input pose**. DOPE is used ONLY for refinement — deliberately decoupled from the TCRen/MJ
  potentials tcren *scores* with (no circularity: don't refine and score against the same quantity).
  DOPE's short-range bins are repulsive, so it supplies its own clash term (no separate clash). The
  restraint keeps it local (without it the search can drift to another favourable pocket). Partners =
  all non-peptide chains within a 12 Å shell. NOT physics MD; Rosetta FlexPepDock is the deferred path.
- DOPE data: `src/tcren/data/dope_potential.npz` (158 atom classes × 29 bins 0.75–14.75 Å), built by
  `scripts/build_dope.py` from the pymod/altmod MODELLER libraries; `tcren.refine._dope()` loads it.
- `tcren refine -s … -o … [--substitute PEP] [--steps N] [--restraint W]`. Native pose ≈ stays
  (RMSD ~0.2 Å); a buried/clashed peptide relaxes locally. Deterministic given `seed`.

## Clash detection + contact stability + register fix — `tcren.clashes`, `tcren.stability`, `tcren.refine.register`

QC for **generated** (AlphaFold/TCRmodel) complexes: their peptide-swap poses are routinely non-physical
(forced poses), which corrupts the contact energy the score reads.
- `interface_clashes(structure, tolerance=0.4, severe=0.6) -> ClashReport` — **native `_geom` kernel**
  (`_geom.interface_clashes`; the numpy path `_clash_pairs_numpy` is kept as the reference/fallback).
  Heavy-atom vdW overlaps (Bondi radii) between the peptide chain and its partners, broken down by
  partner `chain_type` (`by_partner`), with `n_clashes`/`n_severe`/`max_overlap`/`clash_score` + the
  worst residue pairs. `has_clashes(structure)` is the bool convenience. Self-check `python -m tcren.clashes`.
- `contact_stability(structure, cutoff=5.0, delta=1.0) -> StabilityReport`
  (`tcren.mechanics.stability`; `tcren.stability` still imports) — TCR:peptide
  contact fragility read straight off the contact map (native `_geom.contact_stability`, numpy reference
  behind it). Per contact `margin = cutoff − dmin`; report has `n_contacts`, `mean_margin`, `frac_robust`,
  `frac_marg_lt1`, `exp_lost` (expected contacts lost under a `delta`-Å shift). A coordinate-only interface
  positional-confidence readout. Self-check `python -m tcren.stability`.
- `tcren recognize` appends both as output columns: `n_clashes`, `clash_score`, `exp_lost`, `mean_margin`,
  `frac_robust`. **Four of the five are catalogued descriptors now** — `n_clashes` and `clash_score`
  in the `interface` family, `exp_lost` and `mean_margin` in `kinetics` — so they reach the score
  set like any other coordinate rather than riding along beside it.
- `check_register(model, reference=None) -> RegisterReport` — always reports the clash burden; with a
  correctly-registered `reference` (crystal / trusted pose) it adds the **anchor-Cα RMSD** in the
  MHC-groove frame (`peptide_rmsd`) → `wrong_register` (True/False; `None` without a reference).
  **A heavy clash burden alone is NOT a register call** — AF swap models are routinely clashy; register
  needs the reference. (The ila1 CPL lesson: forced pose → raw TCRen ROC 0.35, recovered to ≈0.77 on the crystal register.)
- `fix_register(model, template, engine="ccd") -> ModelResult` — re-threads the model's peptide sequence
  onto `template`'s correctly-registered backbone and re-refines via `model_peptide` (the
  FlexPepDock-functional path; needs the `_refine`/`_fold` kernel). Equal peptide length required.
  Template-free re-docking (CCD to canonical anchor targets) is the future extension.

## 3D alanine scan, both sides — `tcren.ddg.alanine_scan` / `tcr_alanine_scan` (2.25.0)

- **Both sides now scan in 3D.** One residue at a time is truncated to alanine through
  `tcren.refine.substitute.substitute_residues`, the contact map is recomputed and the interface
  rescored. Alanine needs no rotamer: its heavy atoms are exactly backbone + Cβ, so truncating at
  Cβ *is* the alanine.
- `alanine_scan(cm, native, pot, structure=s)` — the **peptide**, one row per position
  (`pos`, `wt_aa`, `ddG`). Without `structure=` it takes the virtual path (no atoms move).
- `tcr_alanine_scan(cm, s, pot, tcr_regions="cdr")` — the **receptor**, one row per *contacted*
  CDR residue (`chain.id`, `chain.type`, `region.type`, `residue.index`, `pos`, `wt_aa`, `ddG`).
  There is no virtual path: truncating a receptor side chain without moving atoms would leave every
  contact it made in place, which is the failure this exists to fix. Scored over `tcr_peptide`
  alone, since a receptor substitution cannot change the peptide:MHC energy.
- `tcr_alanine_reference(scan) -> dict` — `dPhi_ala_cdr12`, `dPhi_ala_cdr3a`, `dPhi_ala_cdr3b`,
  `dPhi_ala_tcr`. Each is the **sum of per-residue** ΔΔGs of that loop, so the three partition the
  total exactly.
- CLI: `tcren ddg -s <in> --native <PEPTIDE> --alanine-scan --side peptide|tcr|both`, plus
  `--virtual` for the fast peptide path.
- `ddG = E(native) − E(Ala@residue)` throughout: **positive = the residue was earning its place**.
- ⚠ **The pre-2.25.0 structural path was wrong.** It threaded the whole peptide through
  `substitute_peptide`, which truncates *every* residue to backbone + Cβ, so each position was read
  against a poly-stub baseline. On 1ao7 the native sequence threaded back through it keeps **14 of
  29** TCR:peptide contacts, and that offset appeared at every position — including positions with
  no contacts, which must read exactly 0. `substitute_peptide` stays correct for the poly-alanine
  *reference*, where every residue genuinely is mutated.
- `substitute_residues(structure, {(chain_id, seq_index): aa})` is the general primitive: any chain,
  any residues, everything else byte-identical. An empty mapping is the identity, which is the
  baseline every ΔΔG rests on.

## Poly-alanine reference score — `tcren.ddg.reference_delta` (geometry-normalized TCRen)

- `reference_delta(cm, peptide, pot, interface="tcr_peptide", reference_aa="A") -> float` = ΔΦ = Φ(peptide)
  − Φ(all-`reference_aa` peptide), subtracting the pose's identity-independent geometry baseline
  Φ(polyAla). It equals the sum of `alanine_scan().ddG` **on the virtual path only**: with a fixed
  contact map the energy is additive over contacts, so the parts sum to the whole. Once atoms move
  they do not — truncating every side chain at once loses contacts each residue alone retains — so
  the 3D scan and the 3D reference are two different measurements. Both are correct; say which.
- **CLI: `tcren score --delta [--reference structural|virtual]`**, one `dphi` column per candidate.
  `structural` (the default) passes `structure=` and rebuilds the poly-Ala contact map;
  `virtual` re-indexes poly-Ala on the template's map. **They are different quantities** — on
  2b4/3QIB under `karnaukhov2022` the references read **-1.3916** and **-0.0166** natural-log units
  against a scored energy of -2.3104, and poly-Ala keeps 62-85 % of the crystal's TCR:peptide
  contacts. Φ(reference) is candidate-independent, so it is computed once per template.
- **CLI: `-c` takes a per-template manifest**, a TSV/CSV headed `complex.id` + `peptide`; each
  structure scores only its own rows. 204 units with their own pools = ONE call and ONE batched
  mmseqs pass, not 204 of each. A plain one-per-line file is unchanged. `tcren rank` reads it too.
- **`interface="complex"` sums both peptide-bearing interfaces** (TCR:peptide with `potential`,
  peptide:MHC with `mhc_potential=`, default MJ), matching `cpl.response_matrix`'s cell convention.
  `ddg`, `neoantigen_ddg` and `reference_delta` all take it; CLI is `tcren ddg --interface complex`.
  **Use it to rank whole library peptides.** `tcr_peptide` alone is blind to presentation — on 1ao7
  the C-terminal anchor swap `LLFGYPVYV -> LLFGYPVYA` reads ΔΔG **0.0000** over `tcr_peptide` and
  **-0.9740** over the complex. The two effects are NOT separable in a library varying every
  position; report the per-interface terms beside the complex.
- **Use for GENERATED poses only.** On a *fixed* contact map ΔΦ = Φ − const → ranking unchanged
  (no-op); it differs only across candidates with their *own* structure (AF swap models). There it
  normalizes out the per-pose interface geometry: **rescues forced/wrong-register poses** whose geometry
  corrupts raw Φ (CPL ila1 TCR-ranking ROC **0.35 → 0.83**), but costs a little where the AF geometry is
  itself informative (mel5/mel8) — so it is a mode, not the default. Sweep: `scripts/ala_reference_sweep.py`
  (manuscript repo). It is **NOT** an affinity ΔΔG (dimensionless contact-preference, no free energy).
- **Ranking, not affinity (ATLAS).** Both raw Φ and ΔΦ are *within-receptor* peptide rankings: they order
  a peptide panel against one fixed TCR:pMHC (Garcia B*27:05 EC50, Spearman ρ≈0.75–0.9) but predict
  equilibrium binding *across* complexes only weakly — on ATLAS SPR raw Φ and ΔΦ track dG/Kd/koff/kon at
  |ρ|≤0.3 (cache-free recompute, `scripts/atlas_tcren.py` + `atlas_within_series.py`, manuscript repo).
  Turning ATLAS Garcia-like (thread a within-TCR panel on one pose) partly recovers the signal — proof
  it's the task setup (within-receptor ranking vs cross-complex affinity), not the potential.

## The background in 3D — `tcren.energetics.background` / `tcren.refine.sidechain` (3.10.0)

- **Why alanine and not "any other residue": `substitute_residues` truncates at Cβ, and alanine's
  heavy atoms ARE N/CA/C/O/CB.** For an alanine target the truncation *is* the mutant — no rotamer,
  no relaxation, no choice. Every other target is a Cβ stub with the right label and the wrong
  reach, which is why the shipped scores reference poly-alanine and why a 20-letter reference had to
  wait for a side-chain builder. What the reference removes is measured: over 370 crystals and
  **8,002** TCR:peptide contacts at 5 Å, only **37.77 %** survive truncation to backbone + Cβ and
  **45.03 %** are made by side-chain atoms on both sides (`bench/eda/sidechain_reach.py`).
- **`smoothed_reference` is a VIRTUAL scan and `background.scan_structure` is the structural one.**
  The first holds σ fixed and re-indexes the potential, so a contact a long arginine reaches across
  is still charged to the alanine replacing it. The second **builds** every candidate at every
  rotamer of a χ grid. Report which; never swap them silently.
- **Three outcomes, not one.** `p = 0` is *no contact*; `p > 0` is *contact*; and every rotamer
  clashing is **collision**, which makes the candidate INADMISSIBLE — it leaves the background
  average with the prior renormalised over survivors, rather than entering at its least-bad rotamer.
  That last is the whole point: `w / w.sum()` reports an impossible substitution as the least-bad
  one and lets it set `F_i`. `reference_free_energy` takes `F_i` over the admissible set `A_i`, and
  `n_admissible` = |A_i| is a read-out — on 1ao7 the A\*02:01 P2 pocket admits **3 of 20** and the
  exposed P8 admits 19.
- **Measured over 60 crystals, 596 positions, 585,720 cells: 14.48 % collision, 81.48 % no contact,
  4.04 % contact**; |A_i| mean **17.17** of 20 (p05 12, p95 20).
- **`sidechain.TOPOLOGY` is chemistry and carries no number; the geometry is fitted** from 374
  Canonical2026 crystals over 176,475 residues into `data/sidechain_geometry.json`
  (`scripts/fit_sidechain_geometry.py`). Round-trip through a residue's own χ: **0.0784 Å** mean
  heavy-atom RMSD over 54,954 residues. `nerf` and `dihedral` invert to 2.4e-14 — **a sign error in
  `dihedral` leaves every bond and angle right and every atom 180° out, reading as a plausible 4 Å
  RMSD rather than an exception**, which is why `_selfcheck` exists.
- **ONE batched call per structure** into `_relax.substitution_scan`; it releases the GIL and threads
  over positions itself. Never per residue, never a Python pool. Agrees with `scan_numpy` to
  **2.2e-16**, 100× faster, **27.5 ms** a structure. Threads run over the `20 P` (position,
  candidate) cells from an atomic cursor, not over positions — the peptide side has only nine —
  and a clashed rotamer is abandoned at the first overlap. Peptide 505 ms serial to 23.1 ms on
  16 cores; receptor 6,372 to 234 ms.
- **The grid is calibrated against one invariant — a deposited residue must be admissible at its own
  position.** `step=60°`, `clash_tol=0.8 Å` reads **0.9725** over 400 residues in 40 crystals;
  30° does **not** beat it (0.9700) at 15× the cost, so the residual is strain, not coarseness.
  `scripts/calibrate_background_scan.py`.

## Potential decomposition + the MJ reference tables — `tcren.potential`

- **A contact matrix is not purely an interaction.** `Potential.decompose()` splits it exactly by
  double-centring: `e(a,b) = mean + H(a) + H(b) + J(a,b)`, `J` with zero marginals. `mean` and both
  `H` depend on one residue each, so an additive per-position model already absorbs them — **`J` is
  the only part a sum over positions cannot express.** Report `J`, not the raw sum, when the claim is
  "this is pairwise chemistry".
- `Potential.hydrophobicity_fit()` → `C0 + C1(q_a+q_b) + C2 q_a q_b` (Li–Tang–Wingreen, PRL 79:765,
  1997: the MJ matrix is nearly rank one). R² = **0.98** on `mj1996`, **0.85** on the bundled `mj`.
  Consequence to state when using such a potential: the interaction term is only `C2·q_a·q_b`, so it
  **cannot prefer one pair of side chains over another of equal hydrophobicity**.
- **Both refuse a directed potential.** TCRen is TCR→peptide; decomposing it is meaningless.
- `Potential.components()` returns that same split as **three scorable potentials**, because an
  interface score is a sum over contacts and the split carries through to it: `size` (the grand mean
  everywhere → `mean × contact count`, an *area* term), `comp` (`H(a)+H(b)`, degree-weighted
  composition) and `pair` (`J`, the interaction). Score a structure with each in turn to say which of
  the three a result is actually reading. **A matrix with no positive entries has a large negative
  mean, so its interface sum is dominated by the contact count** — that is how an interface-area
  effect ends up wearing a chemical name.
- `mj1996()` = MJ 1996 Table 3 raw contact energies (AAindex MIYS960101, transcribed). **Five of its
  210 unique pairs disagree with the AAindex record** by 0.04–0.28 (M–V, D–M, E–M, H–R, A–P; r =
  0.99978). Left untouched, pinned by `tests/unit/test_aaindex.py`; `aaindex("MIYS960101")` is the
  curated alternative. `mj_partition_energy()` = MJ 1985 one-body scale (MIYS850101, larger = more
  hydrophobic; **opposite sign convention** to a contact energy). Cross-check that validates both:
  r = +0.98 between the partition scale and the `q` axis recovered from the 1996 matrix.
- **The bundled `mj()` is Miyazawa–Jernigan 1999, identified 2026-08-29** — AAindex3 `MIYS990106`,
  400 of 400 cells exactly, runner-up off by 0.65. It is **not** 1985 and **not** 1996 Table 3, which
  is what the old "upstream table unrecorded" warning here said. `keskin()` is `KESO980101`, the
  solvent-mediated interfacial form, likewise 400/400 with the runner-up off by 2.77. Both files are
  left byte-for-byte untouched; what changed is that they can be cited. See `SOURCES`.
- **Reference state decides what a comparison measures.** `mj` (−0.079) and `betancourt` (−0.057) are
  *pair-contact* matrices with the one-body term removed; `keskin` (−3.547) and `mj1996` (−3.166) are
  *raw contact energies* that keep it. The like-for-like pairs are `mj` ↔ `betancourt` and `keskin` ↔
  `mj1996`; comparing across the groups compares reference states, not derivations.

### Every published contact matrix — `tcren.potential.aaindex`

`src/tcren/data/aaindex3.txt` is the whole AAindex3 flat file, bundled verbatim (47 records, 80 kB),
so adding a matrix to a comparison costs a string rather than a transcription.

- `catalogue()` → one row per entry with `kind`, `symmetric`, `n_missing`, `mean`, `min`, `max` and
  the citation fields. **Read `mean` for reference state** (above).
- `aaindex(acc)` → one entry as a `Potential`. It **refuses** the 2 contact-*count* tables
  (`TANS760102`, `MIYS960103`) and the 3 side-chain-*distance* tables (`BONM030104`–`BONM030106`),
  because scoring a contact map with a count table is a silent category error; `entry(acc)` still
  returns them deliberately. 42 of the 47 are usable energies.
- `identify(pot)` → accessions ordered by max |Δ|. **An identification needs an exact match AND a
  distant runner-up**; that is how `mj()` and `keskin()` were pinned down.
- Three `ZHAC*` entries are asymmetric by construction (row secondary structure vs column) and
  `decompose()`/`components()` refuse them.
- **AAindex's PMID field sometimes cites the paper that tabulated a matrix, not the one that derived
  it** (`MIYS850102` carries Bastolla 2001). Check the entry's own author/title/journal before
  citing.

## Historical contact-map Potts model — `tcren.potts` / `tcren potts`

Fitted Potts readouts are historical, excluded from final fit-free scoring.
`potts scan` emits positive uncoupled log Z0 (nats) or fixed-native-marginal
first-order response. Its coupled `log_z0` column is a linear surrogate, not the
coupled normalizer. Only full-field centred arms overlap the benchmark; no complete
CLI equivalence exists for composition/support/normalized-contact combinations.
Preserve those original quantities and losses; see the issue-32 disposition audit.

Every other scoring path reads the contact map a structure **has**. This one models the map itself:
the random variable is the whole configuration `sigma`, one binary per *available* residue pair.

```
E(sigma) = - sum_a eta_a sigma_a - 1/2 sum_ab A_ab sigma_a sigma_b ,   P(sigma) = exp(-E)/Z
```

A **site** `a = (i,j)` is a receptor residue and a partner residue whose **Calpha** atoms lie within
`radius` (15 A); `sigma_a = 1` iff a heavy-atom contact formed within `cutoff` (5 A, the unchanged
TCRen definition). The availability mask is **Calpha-only on purpose** — a side-chain- or
type-aware radius would make the single-body field circular. That reference state is the point:
a TCRen potential is a Boltzmann inversion *conditioned on a contact existing*, so a residue that
could have reached the peptide and declined contributes nothing to it; here the non-event is the
observable.

```python
from tcren.potts import PottsModel, available_pairs, fit_potts, score_sites, contact_probabilities

pairs = available_pairs(structure)                 # or partner="mhc" (needs annotate_mhc first)
model = PottsModel.bundled()                       # or fit_potts(pairs, weights=...)
score_sites(pairs, model)                          # energy, log Z, log-likelihood, psi
contact_probabilities(pairs, model)                # p_independent / p_model / p_conditional
```

```
tcren potts fit      -s structures/ -o potts.json --balance both
tcren potts score    -s structures/ -o scores.tsv          # bundled model by default
tcren potts contacts -s complex.pdb -o contacts.tsv
tcren potts contacts -s complex.pdb -m potts_tcr_peptide_sc -o contacts.tsv   # a bundled name
```

`-m/--model` takes a **bundled name** — `potts_tcr_peptide` (default), `potts_tcr_peptide_sc`,
`potts_tcr_mhc`, `potts_pep_mhc` — or a JSON path, on `score` / `contacts` / `map` / `ensemble` /
`scan`.

**Fitting needs no partition function.** `P(sigma_a | rest)` is logistic in
`eta_a + sum_k K_k n_k(a)` with `n_k(a)` the count of contacting neighbours in coupling class `k`,
so the coupled fit is a weighted-binomial GLM with a few extra integer covariates and stays convex
(Besag pseudolikelihood; consistency arXiv:1010.0311, same recipe as plmDCA arXiv:1211.1281).
**Scoring does**, and its reference is exact: `ais_log_z` anneals only the coupling term from
`beta = 0`, where the model IS the uncoupled one and `log Z_0 = sum_a log(1+exp(eta_a))` in closed
form. **Always read `ais_ess`** — effective sample size out of `--particles`; small means the
schedule was too short.

**Gauge: penalise then project.** The design is over-parametrised and identified by an L2 ridge,
*then* projected to zero-sum. An L2 penalty silently picks its own gauge (plmDCA flags this), so
the projection is a separate step, never a constraint. After it, `J` is directly comparable with a
double-centred `Potential`.

**Four coupling families**, all on sequence offsets so no extra coordinates are needed:
`K(di,dj)` within one loop (`|di|,|dj| <= 2`, 12 classes), `L(|dj|, same chain?)` across two
hypervariable loops (6), `M` for the same receptor residue against both partners (1, joint models
only), and `X(same receptor residue, |dj| > 2)` (1, `far=True` only, new in 3.8.0). Everything else
is asserted uncoupled.

**`X` is the exclusion class and it is opt-in.** `K(+0,+1)` and `K(+0,+2)` say a receptor residue
grips a contiguous run of peptide positions; nothing then says the run stops. `fit_potts(...,
far=True)` fits the class that does, at **K = -0.457, s.e. 0.049, z = -9.34** on 367 crystals,
+81.8 nats for one parameter. Every shipped model has `far=False`, so it changes nothing already
written. It bites under a *predicted* Ca geometry, not on the crystal's own.

**The result to know.** On the 362 alpha-beta Native2026 crystals every **axial** class is positive
and every **off-axis** class negative — `K(+1,0) = +0.79`, `K(0,+1) = +0.66` against
`K(+1,+1) = -0.82`, `K(+1,-1) = -0.81`. A made contact recruits its own sequence neighbours onto
the *same* partner residue and suppresses the diagonal one. **This is invisible in raw data**: the
diagonal offsets look positive unconditionally (odds ratio +1.43) and only turn negative once the
axial terms and the Calpha distance profile are held fixed. The couplings buy +505.7 nats of
pseudo-log-likelihood for 18 parameters.

**`--coupling-matrix` is the fair way to compare potentials.** It fixes `J` to one scale on a
bundled matrix, so competitors carry identical parameter counts and identical designs and their
pseudo-log-likelihoods compare directly; and because the matrix is double-centred it contributes
nothing to the one-body marginals, so the 40 field parameters absorb all composition content in
every fit and the comparison is about **pair structure alone**. Measured on Native2026, the ranking
**inverts** between interfaces: TCRen2 beats MJ by 103.3 nats on TCR:peptide, MJ beats TCRen2 by
35.5 nats on the TCR:MHC groove, and TCRen2's fitted scale falls 5.4-fold across that move
(+1.131 -> +0.209) while MJ's barely moves (+0.803 -> +0.974). That is the measurement behind
scoring `Phi_tcr_mhc` with MJ and reserving TCRen for TCR:peptide — **do not reuse TCRen2 on the
groove.**

**`pin_centred=False` is the gauge that reproduces a referenced score.** Double-centring is right
for *comparing* matrices and wrong for *reproducing* one. `reference_delta` is a difference of
one-body sums, so a centred pin re-injects `n_i * c(a)` — the position's contact count times the
potential's partner-residue column mean (s.d. 0.0668 on TCRen2, `n_i` from 1 to 54) — and the
identity fails. `fit_potts(..., coupling_matrix="tcren2", pin_centred=False)` pins the raw matrix,
so the coupling **is** the potential and any linear read-out of the field reduces to the potential's
own score up to the fitted scale.

**The free energy has three limits, and the interface picks one.** `dPhi/d eta_a = p_a`, so
`Phi = sum softplus(eta)` is an interaction sum weighted by contact probability, and a fixed contact
map is the `p in {0,1}` case. Hard contact reproduces `reference_delta`; smoothed (`p` free) is for
a *plastic* interface; saturated (`p` driven to 0/1 by the interface itself) is where a fixed map is
already exact. The groove is measurably in the saturated limit — over four 100 ns trajectories,
36 of the 38 engaged peptide positions (94.7%) reach a maximum groove-pair contact frequency above
0.98, against 14 of 38 (36.8%) on the receptor side. See `docs/potts.rst`.

**Two bundled models**, `PottsModel.bundled("potts_tcr_peptide")` (the default, 64,622 sites /
7,865 contacts) and `"potts_tcr_mhc"` (239,093 / 15,451 — the TCR makes twice as many contacts with
the MHC as with the peptide). Both carry the alpha-beta HARD RULE, as `derive-potential` does.

**Bound versus unbound** (`bound_unbound`, `count_profile`, 2.14.0). `E(empty) = 0`, so the scores
`score_sites` already emits decompose exactly:

```
-E(sigma_obs)  =  log Z  +  L(sigma_obs)
[binding log-odds]  [capacity]  [typicality]
```

`bound_unbound` gives three readings of the whole-interface two-state contrast from ONE Gibbs pass:
`df_empty` = `log(Z - 1)`, exact; `df_threshold` = `log[P(N>=x)/P(N<x)]`, in which `Z` cancels so no
AIS is needed; and `mu_star`, the chemical potential at which the model's mean contact count matches
the observed one (`nan` outside the sampled support — that is not a failure, it is refusing to
extrapolate). They are not competing estimates: no sampler reaches `N = 0` for a docked pose, so
`df_empty` must come from `log Z`. `count_profile` gives the pooled `F(N) = -log p(N)` landscape, so
a threshold is read off it rather than assumed.

**Constraining a statistic of the whole configuration.** `gibbs(..., observer=fn)` calls `fn` on
each kept draw with the `(chains, n)` configuration matrix — that is the hook for accumulating the
kind of statistic a Lagrange multiplier couples to (Jaynes; Tkacik et al. 2014's K-pairwise model is
this with `O` = the total activity). Because such a statistic depends on sigma only through a
low-dimensional summary, one sampling pass serves every step of a moment-matching fit: reweight with
`tilt_mean` rather than resampling. A LINEAR tilt in `N` is exactly a constant added to every field,
`E - mu N = -(eta + mu).sigma`, which is what makes the reweighting checkable against direct
simulation.

**Gotchas.** `partner="mhc"` needs the groove regions, which chain typing alone does not assign —
run `annotate_mhc` / `annotate_mhc_batch` first or every MHC residue is silently dropped for want
of a within-region coordinate. Distance-bin edges are **global** (base 0), never derived from a
frame's own minimum, or a model fitted on one set indexes the wrong coefficients on another.
`psi` (log-likelihood per available pair), not `log_lik`, is the column to compare across
interfaces of different size.

## Does the residue turn toward the peptide? — `g_orient` (3.9.0)

The Ca geometry says a pair is close enough; whether the residue **faces** it is a separate
question the site table could not ask, carrying `d_heavy` and `d_ca` and no Cβ.

```python
model = fit_potts(pairs, coupling_matrix="tcren2", orient=True)   # needs `cos_cb` on `pairs`
```

- `geometry.orientation_cosine` — cosine between the Ca→Cβ direction and the direction to the
  partner Ca. `+1` aimed at it, `-1` away. Held at a fixed Ca distance it separates contact from
  non-contact at **ROC-AUC 0.69–0.84** over 4–12 Å.
- Fitted as `g_orient` (8 bins, 369 crystals, `J` pinned to `tcren2`): **+787.0 nats for 8
  parameters**, against the 3.8.0 exclusion class's +81.8 for one. Monotone across 3.05 log-odds,
  −1.29 at cos ∈ [−1, −0.75) to +1.76 at [+0.75, +1].
- **The Cβ is reconstructed, not read.** `geometry.virtual_cb` builds it from three consecutive Ca
  (0.547 Å mean residual), and `available_pairs` uses the same reconstruction on a crystal — so
  fitting and scoring describe **one** quantity rather than two. Pearson 0.81 against the real side
  chain. Glycine is not a special case: this is where a Cβ *would* sit, so nothing is missing.
- **It pays on a crystal, not on a drawn geometry.** Held out over 24 crystals, on the complex's
  own Ca: js **0.7052 → 0.7290**, ROC-AUC **0.9599 → 0.9691** (the 3.8.0 exclusion class moved js
  0.7133 → 0.7162 and left AUC flat). On a `potts.chain` draw it is neutral (0.4200 → 0.4203) —
  `cos_cb` is a deterministic function of the Ca trace, so it is only as receptor-specific as the
  trace. **Use it when you hold coordinates.**
- **Opt-in on both sides.** A model with `g_orient == []`, or a frame with no `cos_cb`, carries six
  blocks and scores exactly as before. All four bundled models have `g_orient == []`, and
  `SITE_COLUMNS` is unchanged — `cos_cb` rides beside it.

## EM over the borderline contact calls — `tcren.potts.em` / `tcren potts fit --em`

`fit_potts` treats `sigma` as data. A contact is called by thresholding the closest heavy-atom
distance at 5 A, so a pair within half an angstrom of that line is **a coin toss recorded as a
fact**. Those sites are latent, and Monte-Carlo EM is what a latent indicator asks for.

- **E-step** — sample `sigma_latent | sigma_observed` by block Gibbs with every confident site
  CLAMPED. The couplings are why this is worth doing: a borderline pair beside three made contacts
  is imputed differently from one standing alone.
- **M-step** — the existing penalised pseudolikelihood on the completed data, `draws` completions
  stacked at weight `1/draws` so the effective sample size, and therefore what the ridge means, is
  unchanged. Each completion gets a `#m` suffix on its `pdb.id`, because every coupling edge in
  `tcren.potts.kernel` joins on that column and unsuffixed copies would link sites across draws.

```python
from tcren.potts import available_pairs, fit_potts, fit_potts_em, latent_mask, consistency

pairs = available_pairs(structure, probe=5.5)   # NOTE the probe -- see below
ref   = fit_potts(pairs)
model, trace = fit_potts_em(pairs, init=ref, draws=5, steps=10, seed=0)
consistency(ref, model)                          # per block: how far it moved, in what units
```

```
tcren potts fit -s <structures> -o potts_em.json --em --em-draws 5 --em-trace-out em_trace.tsv
```

**THE LATENT BAND MUST REACH BOTH WAYS, and by default it does not.** `available_pairs` measures
`d_heavy` only out to its own `cutoff`, so every pair that MISSED carries a null distance. A band
built on that table is one-sided: on the 17 bundled test structures it selects **44 sites, every
one of them a contact**, so EM could remove contacts and never add them — a systematic bias, not a
small one. Pass `probe=cutoff + margin` (which populates `d_heavy` past the cutoff and does NOT
move `sigma`) and the same band selects **115 sites, 44 contacts and 71 near misses**.
`latent_mask` raises rather than proceeding one-sided.

**This adds a fit; it does not replace one.** The shipped pseudolikelihood models are unchanged and
are the point of comparison. `consistency(reference, other)` scales each block's shift by the
reference fit's own cluster-robust s.e. where it has one (`kernel`) and by the block s.d. where it
does not; a block with no spread (`g_class` on a single-partner model) reads `degenerate` rather
than being divided by ~0.

**Gotchas.** MCEM with a fixed `draws` does NOT converge to a point — it random-walks in a ball
whose radius is the Monte-Carlo noise, so `max_shift` has a floor `tol` may sit below. Read
`max_shift_vs_init` beside it: a run that is only jittering has a flat one. Everything is seeded,
so a run is exactly reproducible regardless. `margin=0` makes every site observed and reduces the
whole thing to `fit_potts` exactly, which is the test that it is an extension and not a rewrite. A
redundancy `weights` dict is keyed on the UNSUFFIXED `pdb.id` and is multiplied by the `1/draws`.

## Self-consistent TCRen — `tcren.potts.selfcon` / `tcren potts selfcon`

**The matrix it produces is shipped: `tcren.potential.tcren2_em()`, `data/TCRen2_sc_potential.csv`,
manifest key `tcren2_em`, and the model pinned to it as `potts_tcr_peptide_sc`. Both are reachable
by name from the CLI** — `-p tcren2_em`, `--tcr-peptide-potential tcren2_em`,
`--coupling-matrix tcren2_em`, `-m potts_tcr_peptide_sc`. **`tcren2()` is still the default
potential** — TCRen2-EM ships beside it,
Pearson r = 0.944 against TCRen2 over the 380 cells against r = 0.867 between TCRen2 and the 2022
matrix, so it is a smaller step than the re-derivation was and the two are still not
interchangeable.

`derive_tcren` tallies the amino-acid pairs of the contacts a crystal **shows**; `fit_potts` reads
the finished matrix back as its pair term. The loop was never closed. `fit_selfconsistent` closes
it: each step recomputes the contact posterior under the current model, re-derives the potential
from those contacts **weighted by their probability**, damps it against the previous matrix, and
refits with the new matrix pinned as the pair term.

```python
from tcren.potts import available_pairs, fit_selfconsistent
from tcren.potential import balanced_weights

pairs = pl.concat([available_pairs(s, probe=6.5) for s in structures])   # NOTE the probe
pot, model, trace = fit_selfconsistent(pairs, margin=0.5, damping=0.5, steps=15,
                                       weights=balanced_weights(markup))
pot.to_csv("TCRen_sc.csv")          # a drop-in Potential: -p TCRen_sc.csv anywhere
```

```
tcren potts selfcon -s <structures> --balance both --margin 0.5 \
    -o TCRen_sc.csv --model-out potts_sc.json --trace-out sc_trace.tsv
```

**`--tol` sits on `max_abs_de`, which is a maximum over 380 cells and therefore has a high
Monte-Carlo floor. Do not lower it without raising `draws`, and do not expect `draws` to rescue an
arbitrarily small value.** Measured on Native2026: the plateau falls 0.00391 (`draws=100`) →
0.00217 (`draws=400`), a ratio of 0.553 against the 0.500 that quartering the variance predicts, so
1e-3 needs roughly 1,800 draws and would measure nothing — across those same steps `r_vs_initial`
moves by 6.2e-5 and `beta_matrix` by 2.1e-4. Read **`rms_de`** instead when the question is whether
the matrix is still moving; averaging over the cells rather than maximising puts its floor about
`sqrt(380)` lower. The loop breaks the moment `max_abs_de <= tol`, so a tol near the floor also
decides how many steps you get.

**Weighting each configuration by its probability costs nothing.** The tally over whole contact
configurations, `sum_delta P(delta | obs) * n_ab(delta)`, **equals** `sum_a p_a` by linearity of
expectation — the same identity behind `energy_ensemble`. So configurations never have to be
enumerated to get the counts; `map_ensemble` lists them as a read-out, not as a dependency.

**The estimator is untouched.** `derive_tcren(weight_col=...)` already took a per-contact fractional
weight, so the log-odds, its contact-fraction reference state and its pseudocount are the shipped
ones. Only the counts change.

**The observation is one additive field on `eta`**, `logit = (5 - d_heavy)/margin`, clipped at +-40.
So the E-step is the existing block-Gibbs sampler with a modified vector — no clamping machinery,
no second code path. One knob spans the family: `margin -> 0` pins the crystal map, `0.5 A` releases
the threshold band, a large margin lets the crystal enter only through the fit.

**The pair term is PINNED, not free** — `J` carries one coefficient (`beta_matrix`) instead of 400.
That is what closes the loop: a free `J` absorbs the pair structure and leaves the re-derived matrix
no route back into `P(sigma)`. `fit_potts(coupling_matrix=)` now accepts a `Potential` object, not
only one of the five bundled names, and a model pinned that way stores the matrix so it round-trips
(`PottsModel.coupling_pinned`).

**THE SAME PROBE RULE AS EM.** Build the table with `probe = cutoff + 3*margin` or `d_heavy` is null
for every non-contact and the field can only ever REMOVE contacts. `tcren potts pairs --probe` and
`tcren potts selfcon` both handle it.

**Measured, Native2026, `--balance both`, 359 structures after the withdrawn drop.** At `margin = 0`
the loop reaches its fixed point at step 2 with `max_abs_de = 0.0` **exactly**, contact mass
7,799.0, and reads **R = 0.9977 against the shipped TCRen2** (mean |delta| 0.0075 over 380 cells,
largest 0.53 on K:K, the rarest pair). The residual is the availability set, not the loop: the 15 A
Ca radius holds **7,799 of the deposited derivation's 7,816 TCR:peptide contacts (99.78 %)** — the
17 it misses sit one or two apiece across 17 structures, side chains touching while their Ca sit
further apart than the radius.

**Gotchas.** This is EM-SHAPED BUT IS NOT EM: the potential update is a re-derivation, not the
maximiser of one objective, so `pseudo_loglik` need not rise and on the bundled crystals it does not
(-337.94, -335.05, then drifting to -335.74 over six steps). Convergence is declared on
`max_abs_de`, which contracts at the damping rate. Read `contact_mass` against the crystal's own
contact total — a loop inventing contacts shows up there — and `mass_from_uncertain` for how much of
the tally came from sites the crystal did not settle. Never mutate a `Potential` in place:
`as_matrix()` caches with no invalidation and the bundled loaders are `lru_cache`d, so build a fresh
one each step (`fit_selfconsistent` does).

## Contact ensembles — `tcren.potts.ensemble` / `tcren potts ensemble`

The marginals settle the **mean** interaction energy and nothing else. `Phi(sigma) = sum_a sigma_a
phi_a` is linear in `sigma`, so by linearity of expectation

```
<Phi> = sum_a <sigma_a> phi_a          for ANY joint over sigma, coupled or not
```

Dependence between sites does not break that. It fixes the two conditions, and **both are easy to
violate**:

1. `<sigma_a>` must be `p_model`, the marginal of the COUPLED model by block Gibbs — not
   `p_independent`, the factorised `sigmoid(eta_a)`.
2. It must be read **at the sequence being scored**. `ddg(weights=...)` and
   `cpl.predict_response(tcr_weights=...)` take ONE weight vector computed on the template and
   apply it to every threaded candidate. A candidate residue enters `eta_a` twice, through
   `h_par` and through `J`, so the ensemble itself moves under substitution:
   `sum_a p_a(template) phi_a(x)` is not `sum_a p_a(x) phi_a(x)`.

```python
from tcren.potts import PottsModel, available_pairs, energy_ensemble, map_ensemble, threaded_ensemble

pairs = available_pairs(structure)              # the sites
m     = PottsModel.bundled()
energy_ensemble(pairs, m)                       # phi_mean, phi_var, phi_free, map_entropy, hamming_obs
map_ensemble(pairs, m, by="position", top_k=20) # the configurations, ranked, with p_hat
threaded_ensemble(pairs, m, candidates)         # per candidate; marginals RECOMPUTED
threaded_ensemble(pairs, m, candidates, coupled=False)   # the shipped frozen-marginal approximation
```

```
tcren potts ensemble -s <structures> -o ensemble.tsv --maps maps.tsv --by position -k 20
```

**What the joint buys over the marginals.** `phi_var` needs `<sigma_a sigma_b>` and is invisible to
any set of marginals; `phi_free = -1/beta log <exp(-beta Phi)>` is a soft minimum over
configurations (`beta -> 0` recovers the mean, `beta -> inf` the best map the geometry admits);
`map_entropy = <E> + log Z` is how many maps the interface is choosing between. On the 14 alphabeta
complexes in `tests/assets/pdb`, **Var Phi runs 0.341 to 1.286 squared TCRen units** — a 3.8-fold
range, so plasticity varies between complexes as much as the mean energy does.

**A configuration needs coarse-graining before it has a probability, and the granularity is
measured.** On 1ao7 (174 available sites, 6,400 draws) the raw site-level map recurs NEVER — 6,400
distinct configurations out of 6,400 draws, each at `p_hat = 1.6e-4`, because the space is `2^174`.
Coarse-grained to which **peptide positions** are engaged, the same draws give **41 distinct
patterns, the most likely at `p_hat` = 0.3219 — and that modal pattern is the crystal's own.**
`by="loop"` (39 groups, 3,308 patterns, top 0.0216) is at the edge of what 6,400 draws resolve;
`by="pair"` is only meaningful where `exact_log_z` can enumerate instead. Across the 14 complexes
the modal pattern is the observed one in **6 of 14**.

**The identity is the cheapest available correctness check.** `phi_mean` must equal
`sum_a p_model_a * phi_a` from an independently seeded `contact_probabilities` run to within
Monte-Carlo error; disagreement means one path is wrong or a chain has not mixed. `energy_ensemble`
carries `phi_mean_from_marginals` so a caller can check it without a second Gibbs run.

**Gotchas.** `anchor="observed"` starts the chains at the structure's own map — that is the
"configurations that do not differ much from the observed one" ensemble, and it is NOT a
Hamming-ball constraint: an unconstrained chain started there still samples `P(sigma)`, and
`hamming_obs` reports how far it travels. A site the potential does not define (TCRen2 has no
`from = C` row, so it is 380 cells and not 400) contributes **zero** and is counted in
`n_undefined` per structure, never swallowed. `threaded_ensemble` is serial over structures and
uses its current candidate-initialization policy (it does not retain every candidate’s final
state as a warm start); `coupled=False` is there so the shipped approximation and
the ensemble mean can be measured against each other on identical rows.

## Ring stacking — `tcren.stacking.ring_stacking` (geometry, NOT an energy)

- A contact potential scores a pair by identity alone, so two rings face-to-face at 3.5 Å score the
  same as the same two residues brushing past edge-on. `ring_stacking(source, cutoff=7.5,
  min_seq_sep=1)` measures what identity cannot carry: `centroid_distance`, `interplanar_angle`
  (0 = face-to-face, 90 = edge-to-face), `vertical` (gap between planes) and `lateral` (sideways
  slide). Parallel-displaced stack = small vertical + a few Å lateral.
- **Pro is in `RING_ATOMS`** although not aromatic — its pyrrolidine ring packs face-on against
  aromatics via CH–π, which is the interaction the module exists to measure. Trp is represented by
  its six-membered ring (the face that stacks).
- Returns **no energy** by design. It says the rings are or are not arranged the way a stack is;
  valuing that is left to whoever has a potential that can.

## Intra-peptide term — `tcren.intra_peptide_energy` (the contacts a chain makes with ITSELF)

- Every Φ in the package sums over contacts between two **different** chains, so a peptide held in its
  bound conformation by its own side chains scores the same as one that is not. This term is that
  omission, and it is **off everywhere by default** — `scope="inter"`, `peptide_internal=False`,
  `intra_weight=0.0` all reproduce the previous output byte-for-byte.
- `all_atom_contacts(..., scope="inter"|"intra"|"all")`; `peptide_internal_contacts(structure,
  cutoff=5.0, min_seq_sep=3)` is the intra case plus the separation filter. Same 5 Å as everywhere
  else, so an internal and an interface contact mean the same thing; `|i−j| ≥ 3` is what does the
  filtering, because sequence neighbours touch by covalent geometry — over `tests/assets/pdb` the
  totals are 18 contacts at `|i−j| ≥ 3` vs 134 at `|i−j| ≥ 2`, the jump being `i`/`i+2` pairs.
- `intra_peptide_energy(cm, pot, peptide=None, contact_weight="residue") -> float`; needs
  `ContactMap.from_structure(..., peptide_internal=True)` (which stores the pairs *beside* `contacts`,
  never in it). `peptide=` threads a candidate onto the structure's peptide positions.
- **The potential is symmetrised**, `(F + Fᵀ)/2`. An intra-chain pair has no from/to orientation, and
  which residue lands on which side is the contact table's canonical `(chain.id, residue.index)`
  ordering, not chemistry. Matters for TCRen (directed), no-op for MJ. **Default MJ**: TCRen is derived
  from TCR↔peptide contacts and says nothing about a chain's contacts with itself.
- Wired in at three layers: `score_peptides(..., intra_weight=w, intra_potential=)` /
  `tcren score --intra-weight` (score = Φ + w·E_intra); `run_pipeline(..., intra_weight=w)` /
  `tcren scoring --intra-weight` (reports `Phi_pep_int` raw, folds w·it into `Phi_total`; potential
  overridable via `potentials={"peptide_internal": …}`); `tcren recognize --full` (`Phi_pep_int`,
  `n_pep_int`, catalogued `involves_tcr=False` — the peptide's own conformation is cohort identity).
- **Expect a sparse term.** A canonical extended class-I 9-mer makes 0–2 internal contacts, so it only
  separates candidates where the peptide is genuinely bulged or self-packed. Untested as a ranking
  signal — it is exposed so the assumption can be measured rather than inherited.

## Interface mechanics — `tcren.mechanics` (koff/kinetics, NOT ΔG)

- The TCR↔pMHC contact map as a network of breakable Cα-anchored Hookean springs (per-contact
  stiffness from heavy-atom-pair multiplicity; default `weight="invdist2"` = multiplicity/dist²).
  Pure-numpy, single structure, no MD. Public API:
  - `interface_springs(structure, cutoff=8.0, weight='invdist2') -> InterfaceSprings(a, b, k, rest, axis)`
    — the TCR-side/pMHC-side Cα anchors, spring stiffnesses, rest lengths, and the unit docking axis
    (TCR→pMHC). Raises if no peptide chain.
  - `stiffness_tensor(structure, cutoff, weight) -> dict` — linear-response `K = Σ kᵢ ûᵢ⊗ûᵢ`:
    `S_tot` (trace), `K_tens` (along docking axis), `K_shear` (`S_tot − K_tens`), `aniso`
    (`K_shear/K_tens`), `lam_max`/`lam_min`, `n_spring`. All `nan` if < 3 springs.
  - `rupture(structure, direction='tensile'|'shear'|'auto', cutoff, weight, break_strain=0.5, steps=80) -> dict`
    — steered-unbinding cartoon (rigidly pull pMHC off, break springs past strain): `rupture_force`
    (peak resisting force), `rupture_work` (∫ force·displacement), `n_spring`, `break_strain`.
    `"auto"` = min-force of tensile/shear.
  - `coupling_residues(structure, cutoff=5.0) -> dict` — residues in both an intra-body scaffold
    contact and the interface: `couple_pep`/`couple_mhc`/`couple_tcr`, `couple_total`, `n_interface`.
  - `interface_mechanics(structure, ...) -> dict` — **the one definition of "the mechanics row"**:
    the union of the three above under their shipped defaults. `tcren mechanics` and
    `tcren recognize --mechanics` both go through it, so the two agree by construction.
- **For a cohort, prefer `tcren recognize --mechanics -t 0`** over running `recognize` and
  `mechanics` as two commands. Both need the same annotated structure, so the second command repeats
  the parse and both mmseqs searches and returns a second table (CSV, keyed `pdb.id` rather than
  `complex.id`) that then has to be joined; the flag reuses the annotation and costs only the
  mechanics arithmetic. The structure must be MHC-annotated either way — an unannotated MHC chain
  silently empties the TCR:MHC half of the spring network (this is what the 2026-07-28 sharding bug
  did: `couple_mhc` = 0 in 523/523 rows, 67 springs against 279).
- **Caveat: these track the dissociation off-rate koff / kinetic stability (Bell–Evans rupture
  resistance ~ r0.5 on ATLAS), NOT the equilibrium ΔG/Kd** — rupture reflects the dissociation
  barrier, not the well depth (physically apt for the TCR mechanosensor / catch bonds). Use them
  **between structures** (one value per complex) to rank/compare; do not pool many per-residue or
  per-spring rows from one structure as independent samples — that is pseudo-replication.
- Self-check (no PDB): `python -m tcren.mechanics`.

## Docking geometry — `tcren.docking.angles` + `tcren.docking.tcrdock_geometry`

- **Two interpretable angles** (existing): `docking_angles(structure) -> DockingAngles(crossing_angle,
  crossing_angle_signed, incident_angle, ...)`. `crossing_angle` = the groove-plane "scanning" angle,
  `incident_angle` = the tilt. Computed from the Vα→Vβ axis in the groove frame; no reference DB.
- **Full rigid-body pose** (`tcren.docking.tcrdock_geometry`): `docking_geometry(structure) -> DockingGeometry(d,
  torsion, tcr_unit_y, tcr_unit_z, mhc_unit_y, mhc_unit_z)` — native reimplementation of **TCRdock**
  (phbradley/TCRdock, MIT, commit `c5a7af4`; see `THIRD_PARTY_NOTICES.md`). MHC + TCR symmetry stubs (β-sheet
  floor / Vα-Vβ two-fold), MHC-I core by BLOSUM-align to TCRdock's template, TCR core by conserved IMGT
  framework positions from arda region markup. Needs a chain-typed + MHC-annotated structure; **class-I only**
  (class-II raises).
- **What we use / provenance finding (2026-07-05, validated on 618 TCRvdb models):** the upstream AF/TCRmodel2
  annotation table's `scanning_angle` **is** reproducible (= `crossing_angle`, r≈0.88), but its
  **`pitch_angle` is NOT** any clean geometric angle (best correlate `d`, r≈0.42) and out-discriminates every
  clean docking feature on TCRvdb (macro-PR≈0.72 vs `d`≈0.64 / `torsion`≈0.62 / tilt≈0.58) → its extra signal
  is **AlphaFold-confidence contamination, not geometry**. Prefer the documented `d`/`torsion`/`crossing`/
  `incident` over the opaque upstream `pitch_angle`.

## Publication figures — `tcren.viz.pymol` (headless PyMOL)

- `render(scene, png, size=(1200,1200), dpi=300, corner='bottom-left', gizmo=True)` ray-traces a
  PyMOL script body. Scene presets: `overlay_scene(ids, canon_dir, limit=8)` (ensemble, side-on),
  `groove_scene(pid, canon_dir, surface=False)` (peptide in the cleft, top-down; `surface=True` is
  the histo.fyi look), `interface_scene(pid, canon_dir, cdr_resi)` (peptide + CDR loops).
- **`CANONICAL_AXES` is the naming of `tcren.docking.frame`, and the gizmo prints it**: x=`width`
  (groove width, α1↔α2, PC3), y=`N→C` (groove axis toward the peptide C-terminus, PC2), z=`TCR`
  (docking normal, MHC floor→TCR, PC1). Same three directions as SwiftTCR / TCR3d; the PC *ranking*
  differs because `orient.frame` fits the whole complex and they fit the MHC groove alone.
- The gizmo is a **separate render pass composited at pixel coordinates**, not a CGO projected into
  the corner. PyMOL's orthoscopic viewport does not span `dist * tan(fov/2)` — measured, it is out
  by ~25% — so the projection route puts it off-frame. Do not "fix" this by reintroducing the math.
- Two PyMOL gotchas the module already handles: `auto_zoom` is **on** by default, so loading a CGO
  or a label pseudoatom re-frames the camera onto it (the module sets it off first); and
  `cmd.get_view()[0:9]` is **column-major**, so the transpose that inverts it is `rot[i*3+j]`.
- **`residue_importance(structure, interface='tcr_peptide') -> pl.DataFrame`** — the exact
  decomposition of Φ: `chain.id, residue.index, residue.aa, region.type, n_contacts, phi`, sorted
  most-favourable first. `importance_scene(pid, canon_dir, imp, by='phi'|'n_contacts')` paints it
  (sticks on a ramp, everything else pale). **The phi ramp is centred on zero**, so blue/red mean
  favourable/unfavourable rather than less/more — a range-fitted ramp reddens the least-favourable
  residue even in an all-stabilising interface. **Each contact is attributed to BOTH residues it
  joins, so the per-residue phi sums to 2x the interface total** — an attribution, not a partition;
  a test pins that factor of two.
- **`groove_importance_scene(pid, canon_dir, importance, by='p_expected', spectrum=..., surface=False)`**
  is `groove_scene`'s framing — peptide threaded along the cleft, seen from above, helices and floor
  kept as furniture — with the peptide coloured by a per-residue value carried in on the B-factor
  column. Two differences from `importance_scene` and both matter: it **zooms the groove, not the
  coloured selection**, so several complexes rendered this way are directly comparable; and it ramps
  over the column's own range, which is right for a one-sided quantity such as
  `potts.contact_map(by="position")`. Use `importance_scene` for a signed share, which centres on
  zero.
- Chain roles after `tcren orient`: A=Vα, B=Vβ, C=peptide, D=MHCα, E=MHCβ/β2m — `CHAIN_COLOURS`.
- Notebooks: `notebooks/pymol_canonical_figures.ipynb` (static gallery) and
  `notebooks/pymol_interactive.py` (**marimo**, `pip install 'tcren[marimo]'`). marimo is reactive,
  so a name may be defined by exactly one cell — import shared symbols once in the setup cell and
  thread them through the signatures, or the export fails with `MultipleDefinitionError`.
- Docs gallery page with rendered examples: `docs/gallery.rst` + `docs/_static/gallery/*.png`.
- Everything except `render`/`probe_rotation` is pure Python and unit-tested without PyMOL
  (`tests/unit/test_pymol_viz.py`); the two that shell out are `slow` + skip when pymol is absent.

## Contact typing — `tcren.contact_types` (two schemes; v1 is frozen)

- **`scheme="v2"` (default)**: `salt_bridge, hydrogen_bond, cation_pi, stacking, aromatic, hydrophobic,
  polar, vdw, other`. `other` now means only "too far to be anything", never "unrecognised" — it fell
  from **72.3% to 13.9%** of TCR:peptide contacts on five crystals. A contact carries *several* types
  (`is_<type>` booleans are independent); `contact.type` is only the top-priority label.
- **`scheme="v1"`** is the old five-type residue-level scheme, kept byte-for-byte. `recognition.py`
  **pins it** — the frozen classifiers were fitted on its `ct_*` counts. Do not change v1.
- What v2 fixed: apolarity is decided per *atom* (a carbon with no bonded N/O), not per residue —
  v1's residue set excluded Tyr, the commonest TCR interface residue. H-bonds reach 3.9 Å with
  donor/acceptor typing (two carbonyl oxygens are no longer an H-bond). `stacking.ring_stacking` is
  finally joined via `stacked_pairs(structure)`.
- `residue_pair_types(structure, interface)` types from **every** atom pair, not the closest one —
  a salt bridge whose nearest contact is two carbons is invisible otherwise. `type_weights(typed)`
  gives 0/1 weights that drop pure-proximity contacts (`tcren score --drop-untyped`).
- **Hydrogens are now filtered** in `all_atom_contacts`. This changes contacts and energies for
  H-bearing depositions (5jhd: +7 of 28 contacts, −58.5% `Phi_tcr_pep`) and breaks legacy-oracle parity
  on 5jhd/7qpj, recorded as a subset relation in the regression test.

## Feature families — `tcren features` (descriptors) vs `tcren assess` (scores)

- **Two commands, two jobs.** `tcren features` reads structures and writes descriptors;
  `tcren assess` turns a descriptor table into the score set. The feature pass is the expensive
  half, so run it once and re-score for free:

  ```bash
  tcren features -s <in> -i placement,interface,topology,energetics -o feats.tsv
  tcren assess --features feats.tsv -o scores.tsv
  ```

- **The score set (3.0.0), and every one of them is defined for a SINGLE structure.** The
  transform, the class means and the covariance are all frozen on a hold-out that ships with the
  package, so nothing is estimated from the rows being scored. Higher is better throughout.

  | column | tier | what it answers |
  |---|---|---|
  | `pose_score` | 1 | is this the kind of interface real complexes make? No negative, no label. **The bad-pose channel.** |
  | `binder_score` | 2 | log-odds that the complex is a genuine recognition interface |
  | `channel_*` | 2 | the same log-odds marginalized to one family — **which part of the structure says so** |
  | `peptide_score` | 0 | the poly-alanine-referenced energy. Nothing fitted. Ranks PEPTIDES for a fixed receptor |
  | `confidence_residual` | 1 | reported ipTM minus what the coordinates say it should have been |
  | `binder_iptm` | 2 | `binder_score + logit(ipTM)`; two log-odds added, no coefficient. The recommended read when a confidence exists |

  Five channels, named in physics and geometry terms: `placement` (where the receptor sits),
  `interface` (how much it makes, of what chemistry), `shape` (the footprint free of its size),
  `energetics` (contact chemistry in kT), `mechanics` (breakable springs). A marginal of a
  Gaussian is a sub-block of Sigma — exact, no re-fit — so a channel costs an index. **A channel
  beats the whole model where the whole model dilutes it**: on template-free cohorts
  `channel_shape` reads 0.637 against the full posterior's 0.615.

- **Pass `--peptide` when the PEPTIDE is what varies** (a combinatorial library, a mutational
  scan). Otherwise the five descriptors computed without the receptor are marginalized out,
  because they are constant across every structure of one epitope on one allele and a model
  reading them reaches the cohort's name without reading an interface.

- **Which read-out answers which question.** There are **two ranking tasks**, and magnitude is a
  level rather than a third task. Read down this list and stop at the first row that matches:

  | the question you are actually asking | what to read | how to get it |
  |---|---|---|
  | *Which peptide does this fixed receptor read?* | `peptide_score` | `tcren assess --peptide`, or `tcren ddg` / `tcren cpl` directly |
  | *Which receptor recognises this fixed epitope?* | `binder_iptm` if the table carries ipTM, else `binder_score`, reported beside `S` | `tcren assess` |
  | *Is this pose real at all — I have no label and no cohort?* | `pose_score` | `tcren assess` |
  | *Is the generator's confidence warranted?* | `confidence_residual`, then `p_nonbinder_af` | `tcren assess` on a table with an `iptm` column |
  | *Which part of the structure says so?* | `channel_*` | `tcren assess` |
  | *Which residue carries it?* | `residue_deltas` `delta` | `tcren explain --score <read-out>` |
  | *How fast does it come apart?* | the mechanics read-outs — a bound structure reaches the **off-rate**, not the equilibrium constant | `tcren mechanics` |

  Two traps this ordering exists to avoid. `peptide_score` **reads below chance on a receptor
  benchmark**, and that is a property of its poly-alanine reference frame rather than a fault — it
  is the wrong instrument for the second row, not a broken one. And the two composite arms
  (`S` + `binder_score`) z-sum quantities on different natural scales, so they need a **cohort** to
  set that scale; every other read-out above is defined for one structure on its own. Report the
  two modes apart.

- **The fit-free predecessor tier is still shipped**, and `tcren recognize --features feats.tsv`
  emits it: `complex.id` + **`Q`** (interface quality), **`T`** (the shape block) and **`S`**
  (their composition with the energy block). `S` leads the functionally validated receptor screen
  on its own and **composes** with `binder_score` rather than being replaced by it.
- **What the composite is called, and what is gone.** `S` is the current name: 2.27.0 renamed
  `reliability.s_free` → **`reliability.s_score`**, the emitted column `S_free` → **`S`**, and the
  frozen key `<set>|S_nat` → `<set>|S`, **outright, with no alias** — a caller on the old name fails
  loudly rather than reading a renamed quantity. `P_native` with `cohort.p_native` and the
  per-channel posteriors `G`/`T`/`E` went at **2.26.0**; `reliability.p_binder` with
  `available_links`, and `reliability.correct_confidence` with `available_corrections` and the
  `tcren diagnose` command, went at **2.28.0**, because the Platt links and the four correction
  coefficients were fitted out of fold. None of them returns.

- **`tcren fit-holdout` reproduces the frozen model to a relative 1e-5** — not bit for bit — from
  the shipped manifest plus a descriptor table you regenerate with `tcren features`. The 19 MB
  descriptor table does not ship; the 8,292-row manifest does. This is the contract that `P_native`
  could not offer, and 1e-5 is the bound the optimiser earns rather than a round number (see the
  reliability section for the measured spreads).

- **Six families, split by invariance** (`tcren.recognition.DESCRIPTORS`, `FAMILIES`); a
  *family* is a slice of the descriptor table and a **channel is one of the five score channels**
  in `tcren.score.CHANNELS`, with `tcren.score.CHANNEL_OF` the map between them — `topology` →
  `shape`, `kinetics` → `mechanics`, and `potts` joining `energetics` because both are contact
  energies in kT and splitting them helps nobody. Families:
  `placement` (groove-frame pose — angles, TCRdock params, ride height/shift/offset, CDR3 frames;
  frame-**dependent**), `interface` (contact size + chemistry), `topology` (the *shape* of the
  contact set, size-free), `energetics` (Φ and ΔΦ), `potts`, `kinetics` (spring network; off by
  default). **164 catalogued descriptors as of 2.30.0** — `topology` 70, `placement` 31,
  `interface` 26, `kinetics` 17, `energetics` 15, `potts` 5. `topology` is the largest since it
  absorbed the contact-graph, matrix and surface blocks.
- **The gap is in `topology`, and it is the one channel that measures space rather than incidence.**
  `tcren.topology.surface.surface_map(s, side=...)` rasterises each face as a height field on one
  shared groove-frame grid, so `gap = h_tcr - h_pmhc` cell by cell. The sign is counter-intuitive
  and is the point: the median gap is **-1.7 Å** and **71 % of cells interdigitate**, the faces
  interlocking rather than stacking. Reported pooled (`sc_gap_mean`, `sc_gap_sd`), integrated with
  the signs apart (`sc_gap_vol` void vs `sc_interlock`, both Å³) and resolved by sign
  (`sc_interlock_frac`, `sc_gap_depth`, `sc_gap_height`, `sc_gap_asym`). Measured on 4,907 labelled
  benchmark structures, `sc_gap_mean`/`sc_gap_sd` carry the largest binder contrasts of any
  published interface descriptor tested (Cohen's *d* -0.651 / -0.681) at R² 0.131 / 0.255 on all
  141 pre-2.30 descriptors.
- `placement` + `interface` were one `geometry` family and `energetics` was `physics` until
  2026-08-24. Both retired names still resolve in `descriptors()`; the split is what lets the
  independence claim be stated — measured on VDJdb, topology ⟂ interface at |ρ| = 0.023 while
  topology–placement is 0.177 (0.448 on TCRvdb), because uniform coverage *is* ride height.
- **Only the requested families are computed.** `-i topology` never builds the energies.
- **Contact counts are `interface`, not `topology`** (`FOOTPRINT_SIZE_FEATURES`). A shape channel
  carrying the interface's size would correlate with the interface channel by construction.
- **A channel costs an index, not a fit.** Each of the five is an exact marginal of the one frozen
  covariance — a sub-block of Σ, closed form — so `channel_scores` re-fits nothing and the five
  land on one scale and can be read against each other. **They do not sum to `binder_score` and
  should not**: the whole model also reads the correlations *between* channels, which a
  per-channel view cannot show. A channel is only emitted when at least 3 of its transformed
  coordinates are present.
- **`P_native` and the whole latent-class layer are gone (2.26.0), and are not to be re-added.**
  `cohort.p_native`, `P_NATIVE_CHANNELS`/`_POOL`/`_ORIENT`/`_FEATURES`/`_BANNED`, the per-channel
  posteriors `G`/`T`/`E`, `GaussianBNClassifier` and `BayesianLogisticRecognizer` with
  `encode_features`, `frozen_recognizers` and `real_probability` all went in one release. It
  re-fitted on every call, **raised when a cohort had fewer rows than features**, and its value
  depended on which rows the fit was anchored on — none of which survives contact with a user
  holding one model. The score set replaces it at a strictly better standing, and every read-out
  above is defined for a single structure.

## Single-structure reliability — `tcren.score`, `tcren.reliability` / `tcren assess`

**Every read-out here is defined for ONE structure.** The score set's transform, class means and
covariance are frozen on the shipped hold-out; the fit-free tier's blocks are frozen on the
Native2026 crystals. Nothing is estimated from the rows being scored, so no score moves depending
on what was scored beside it.

```python
from tcren.score import (peptide_score, pose_score, confidence_residual, binder_score,
                         channel_scores, score_table, holdout_model, holdout_manifest,
                         ScoreModel, MODEL_FILE, CHANNELS, CHANNEL_OF)

m = holdout_model()                       # the frozen ScoreModel; MODEL_FILE is where it ships
out = score_table(features, receptor=True, iptm=iptm, model=m)   # every read-out, one frame
b   = binder_score(features, receptor=True)      # receptor=False when the PEPTIDE varies
ch  = channel_scores(features)                   # {channel: log-odds}, keys from CHANNELS
r   = confidence_residual(features, iptm)        # reported logit(ipTM) - E[logit ipTM | x]
```

- `score_table` adds **`binder_iptm`** when `iptm` is passed: `binder_score + logit(ipTM)`, two
  log-odds summed with no coefficient to fit. On the functionally validated receptor screen it
  reads 0.771 ROC-AUC against the posterior's 0.768 and ipTM's 0.795, and it is the recommended
  read when a confidence exists.
- `peptide_score` needs only `dPhi_tcr_pep` and `dPhi_pep_mhc` and fits nothing; every other
  read-out goes through `holdout_model()`, which is `lru_cache`d and raises a message naming
  `tcren fit-holdout` rather than a bare `FileNotFoundError`.
- **Two kinds of absence, handled differently and both deliberately.** A missing *column* is
  marginalized out exactly (a sub-block of Σ), which is what lets `-i placement,interface,topology`
  still produce a binder score; a *row* missing a column that is present is left NaN rather than
  imputed, because a filled-in value is a fabricated correlation. Below `ScoreModel.MIN_COORDINATES`
  = 20 usable coordinates the model raises and names the `tcren features` call that fixes it.
- `holdout_manifest()` returns the 8,292 structures the fit used with dataset, epitope, label and
  ipTM. The descriptors do not ship (8,292 × 147 is 19 MB); `tcren fetch-data` → `tcren features` →
  `tcren fit-holdout` returns the shipped arrays to **a relative 1e-5, not bit for bit**. That
  bound is what the optimiser earns rather than a round number: each of the 55 Yeo-Johnson lambdas
  comes from `scipy.stats.yeojohnson_normmax`, whose Brent search stops at its own tolerance of
  about 1.5e-8, and every array fitted through the transform inherits it. Measured between macOS
  arm64 and Linux x86_64: one lambda differs by 8.4e-8 absolute, the class means by at most 6.0e-9
  on 18 of the 149 coordinates, the covariances by at most 2.5e-8 on 50 of the 22,201 entries. A
  test that asserted rtol 1e-9 was green locally and red on CI, which is how this was found.

**The fit-free predecessor tier**, still shipped and still reported:

```python
from tcren.potts import score_sites, available_pairs, PottsModel
from tcren.reliability import s_score, t_score, af_band

potts_scores = score_sites(available_pairs(structure), PottsModel.bundled())  # writes neg_energy
v = s_score(feature_table, energy=potts_scores["neg_energy"])  # energy optional
t = t_score(feature_table)                                     # the shape block on its own
b = af_band(meta["iptm"], reference="binder_bm|ipTM")          # the generator diagnostic
```

`S = Q/sd_Q + T/sd_T + (Pi - mu)/sd_Pi`. Three fit-free directional blocks `z(x)' C^-1 s` over
the Native2026 crystals, each divided by its own native spread.

**The outer transform is a DIVIDE, not a z.** A block score's native mean is 0 by construction, so
re-centring is a no-op; its variance is `s' C^-1 s`, which is not 1 — measured 1.43 (`Q`), 1.61
(`T`), 14.13 (`Pi`). Without the division the energy would carry ten times the weight of the
geometry. Equal weight *in native-sd units* is the claim.

`Pi` is `neg_energy` from `tcren.potts`, the interface energy read against the partition function
rather than a poly-alanine reference. It is the frozen choice because it is the **least redundant
with `Q`** of the five ways of spending `-E = log Z + L`: native Pearson +0.33, against +0.75 for
the contact count and +0.51 for `log Z`. Pass `energy=None` and the two-block form is returned —
still defined, and reported as such rather than imputed.

`t_score` is the block that **survives without a template**: on the balanced VDJdb panel
`T` loses 0.06 ROC-AUC where the epitope has no solved complex, against `Q`'s 0.24. It needs
`q_score(..., signs=T_SIGNS)`, because the footprint's connected-component fraction at 7 A runs the
other way from the rest.

**There is no calibrated probability any more, and that is on purpose.** `p_binder` and its
`available_links()` were removed at 2.28.0 with the `calibration` section of
`data/reliability_moments.json`: the Platt links were fitted out of fold — leave-one-epitope-out on
the 22-cohort VDJdb panel, within-epitope 5-fold on TCRvdb — and shipped as fold means. `moments()`
asserts in its own test that `calibration` and `corrections` are absent, so neither can return
quietly. `available_bands()` lists what still ships: `binder_bm|ipTM`, `binder_bm|pLDDT`,
`tcrvdb|ipTM`, `tcrvdb|pLDDT`. A band table is a quantile binning of the generator's confidence
carrying the observed non-binder fraction with a Wilson interval — a measurement of the benchmark,
not a model of it.

```
tcren features -s models/ -i placement,interface,topology,energetics -o feats.tsv
tcren potts   score -s models/ -o potts.tsv        # writes neg_energy, keyed on pdb.id
# join potts.tsv's neg_energy onto feats.tsv on the structure id, then:
tcren assess --features joined.tsv -o assessed.tsv
```

**Gotcha.** `tcren features` does not compute `neg_energy` — it comes from `tcren potts score` and
has to be joined. Skip the join and `assess` emits the two-block `Q + T` form; it says so in its
report, so read the first line rather than assuming three blocks.

**The forced-pose flag — why a confident model can be wrong in a readable way.** A generator pushed
into a confident but wrong pose does not build a random interface. To seat the chains it picks
residue pairs it believes are favourable, so the recognition energy comes out *good* — often better
than a genuine complex of the same epitope. The energy **inverts** under forcing instead of
degrading: on a 24-cohort forced-pose panel (1,707 structures) it reads macro ROC-AUC **0.4952** and
is below 0.5 in **15 of 24 cohorts**, where ipTM reads 0.6093.

```python
from tcren.reliability import inversion_flag, screening_yield
f = inversion_flag(feature_table, energy=potts_scores["neg_energy"])  # NaN without the energy
y = screening_yield(v, budget=0.10, prevalence=0.48)   # what testing the top 10% implies
```

`inversion_flag` is the energy block minus the mean of the two shape blocks, in native-sd units.
Large positive = the energy is vouching for a structure the footprint does not, which is the pattern
to distrust; a generator fakes favourable contacts far more easily than a well-formed footprint.
It ranks and triages; nothing in the package calibrates it to a probability.

`screening_yield` returns the cut only: how many structures, at what threshold and percentile, plus
`expected_hits` under a stated prevalence. **Enrichment over random is not returned** — it needs
labels the function does not have, and a NaN there would read like a measurement.

**`artefact_directions(binder_coords, crystal_coords, tolerance=0.05)` — which directions belong to
the generator, read with no binder label.** A direction real crystals break *harder* than modelled
binders do is the generator's own regularity rather than the interface: the tightest band breaks
4.55× against decoys' 3.18×, and scores ROC-AUC 0.504 over the 16 template-free cohorts — a coin —
where the loose band reads 0.606. The shipped model cannot see it, because Ledoit-Wolf shrinkage
floors its smallest direction at s.d. `ARTEFACT_SD` = 0.0797, above that band; measuring the bands
needs `Joint.fit(shrink=False)`, which nothing else should use.

**`tcren assess` emits four blocks, one table.** The score set (`pose_score`, `binder_score`, the
`channel_*` columns, `peptide_score`, and with an ipTM column `confidence_residual` and
`binder_iptm`); the fit-free tier (`S`, plus `inversion_flag` when the energy column is present);
triage on the recommended score (`rank`, `percentile`, and the mean over the top `--budget`
fraction); and the generator diagnostic (`af_band`, `p_nonbinder_af`, `s_roc_in_band`). The last is
the actionable one: on the balanced VDJdb panel the **top ipTM decile is 26.2 %** [18.7, 35.5]
**non-binders**, and is also the band where `S` reads highest. Flags: `--peptide`, `--model`,
`--iptm-column`, `--band`, `--budget`, `--list-bands`. **`--link` and `--list-links` went with
`p_binder` at 2.28.0**, and there is no `tcren diagnose` — `confidence_residual` replaced the
fitted confidence correction at tier 1, reading no label anywhere.

## Per-residue attribution — `tcren.score.residue_deltas` / `tcren explain` (3.2.0)

The five channels say which *part of the structure* a read-out comes from; this says which
**residue**, for any read-out and any channel.

```python
from tcren.score import residue_deltas
d = residue_deltas("model.pdb", score="binder")          # or "pose", "peptide", or a channel name
```

```bash
tcren explain -s model.pdb --score binder -o deltas.tsv
tcren explain -s model.pdb --score energetics --chain-types PEPTIDE -o pep.tsv
```

One row per interface residue: `chain.id, residue.index, residue.pdb, residue.aa, score.full,
score.without, delta`, where **`delta` = L(x) − L(x without residue i)** — the score of the complex
minus the score of the same complex with that residue's atoms removed. Positive means the residue
carries the score. The table feeds `viz.pymol.importance_scene` unchanged, so one implementation
colours a structure by `binder_score`, `pose_score` or any single channel; before this, only the
peptide read-out had a colouring at all.

- **Why leave-one-out and not a decomposition.** `peptide_score` is a sum over contacts and
  `energetics.scoring.position_profile` already splits it **exactly** — prefer that one there. The
  other four are functions of ~149 whole-structure scalars through a full covariance per class, so
  no residue owns a share of one and there is nothing to split. The leave-one-out difference is
  defined for all of them.
- **Chain typing and the MHC call run once for the whole structure.** Every masked copy is
  featurised with `annotate=False` against that one annotation, so **no mmseqs call is made per
  residue.** Cost is one descriptor pass per interface residue, about 2 s each, so a complex is a
  few minutes.
- Only residues touching a requested interface are masked. This is a restricted attribution;
  noncontact residues can affect whole-chain geometry descriptors. `chain_types` filters reported
  residues without changing the score mode. Python callers may pass `receptor=False` explicitly
  for the full peptide-presentation score; the default and CLI use receptor mode.
- A masked row whose descriptors go undefined is reported as a **null** delta, never imputed:
  `ScoreModel.coordinates` drops such a row on purpose, and that is surfaced rather than filled in.
- `importance_scene` handles the null: since 3.3.0 a non-finite value is dropped and drawn as pale
  cartoon rather than formatted into `cmd.alter` as `b=nan`, which PyMOL evaluates as Python and
  raises `NameError` on — with its own error printer swallowing the reason.

## Footprint shape — `tcren.footprint` / `tcren footprint`

- **Reach for it when the energy is at chance but the pose still looks wrong.** It reads the same
  contact map as a *shape*, not a sum: how evenly the six CDR loops spread their contacts, whether
  the germline/CDR3 division of labour holds, and whether the footprint is one connected patch.
  No potential, no reference structure, no fitted parameter.
  CLI: `tcren features -s <in> -i topology -o <out>`, then `tcren recognize --features <out>`
  for the shape block `T`. (`tcren footprint` is the same code path, now hidden and superseded;
  since 2.26.0 its `--score` emits the fit-free `reliability.t_score` in place of the removed
  `fp_score` z-sum, so there is no reason to prefer it.)
- **The MHC pass must run AFTER chain typing, and it is not optional.** `classify_chains` leaves an
  MHC chain typed generically `"MHC"`; `interface("tcr_mhc")` matches the supertype `annotate_mhc`
  assigns. Skip it and six of the twelve cells are structurally unreachable with no error —
  `p_germ_mhc` reads 0.06 instead of 0.78 and `H_cell` is normalised by ln 12 over a partition half
  of which can never be occupied. `cell_counts` now warns; that bug shipped once.
- **NaN must become null before any polars aggregation.** polars propagates NaN through
  `mean`/`std`, so one contact-free structure turned a whole z-scored channel to NaN, `fill_nan(0)`
  flattened it, and `fp_score` silently became its other channel alone (0.815 → 0.691 on TCRvdb).
- **No canonical orientation is needed** — every feature is invariant under rigid motion, which
  `test_footprint.py::test_every_feature_is_invariant_under_rigid_motion` pins by rotating and
  translating a complex and demanding an identical row. Only chain typing + CDR markup. MHC *region*
  markup is **not** used, so the "MHC needs a second pass or it silently empties" trap does not apply.
- **Annotation goes through `_iter_typed`/`iter_annotated_set`** — one mmseqs call per organism for
  the whole set. Do not call `classify_chains` per structure here (that is what `tcren surface`
  does, and it is an order of magnitude slower over a cohort).
- `footprint_features(s) -> dict` (29 features at the default two radii), `footprint_batch(paths_or_structures) -> pl.DataFrame`.
  For the shape score use `reliability.t_score(table)` (`T`), which fits nothing at call time;
  the old `footprint_score` z-sum was removed at 2.12 and the latent-class `T` posterior that
  briefly replaced it went with `P_native` at 2.26.0.
- **Coverage**: cells are the 6 CDR loops × {peptide, MHC} (12) or with the peptide split into
  thirds (24). `H_cell` is the normalised Shannon entropy, `D1`/`D2` the Hill numbers — `D2`, the
  *effective number of engaged cells*, separates better because it discounts the rare cells a decoy
  populates weakly. Refining the **peptide** side helps; refining the MHC side into helices/floor
  does not, which is why that partition is not offered.
- **Topology**: `_flag_betti` builds the flag complex on the contacted pMHC Cα at a radius —
  `fp_b0_*` patches, `fp_b1_*` holes (via a GF(2) rank of the triangle boundary), `fp_chi_*`.
  `h0_pers_ent` is the H₀ persistence entropy, whose barcode **is** the MST edge lengths, so no
  filtration library is needed. `b0` is most informative at 7 Å and `b1` at 8 Å; both ship.
- **The bipartite contact graph's cyclomatic number `E − V + C` is deliberately not offered.** With
  ~30 contacts among ~30 residues it is dominated by `E` and just tracks interface size; the patch
  count is scale-free instead. If someone asks for "the Betti number of the interface", this is the
  distinction to make.
- `n_loop_contacts` counts only what the partition sees — the six CDR loops. Framework contacts are
  excluded by construction, so it is smaller than the full interface count. The topology features
  are **not** restricted this way: they use every contacted pMHC residue.
- **It is not `n_contacts`, and the two are counted differently.** `n_contacts` belongs to the
  `potts` family — the available pairs that engaged, 29 on 1ao7 against this module's 66 — and it
  is the column the frozen Potts moments in `data/reliability_moments.json` are standardized on.
  Through 2.19.0 both passes wrote the one name and the later pass won, so a table carrying the
  loop tally under `n_contacts` is off by a factor of two to four with no error and no NaN. The
  guard `reliability._check_potts_contacts` was written for `correct_confidence` and **lost its
  only caller when that went at 2.28.0**; what refuses a stale table now is the provenance digest
  (`tcren.provenance.check`), which rejects any feature table written under a different catalogue
  and names the command that regenerates it.

## Surface topology — `tcren.surface` / `tcren surface`

- **Reach for it when the question is about the pMHC alone**, with no TCR in the structure or before
  one is docked: how featureless is this epitope, which epitopes present a similar face, is the
  peptide bulged. A contact potential cannot answer any of these — it needs an interface that exists.
  CLI: `tcren surface -s <in> -o <out> [--compare] [--svg <dir>] [--cells <csv>] [--channel h]`.
- `surface_map(structure) -> SurfaceMap` — height + hydropathy + charge on a 64×32 raster over the
  groove; `surface_stats` gives `relief`, `peak_to_valley`, `frac_above_ridge`, `phobic_centre`;
  `surface_distance`/`surface_tree` compare epitopes (SURFMAP Manhattan distance + linkage).
- **The frame is refit per structure** — z from the groove-floor plane normal, **y from the peptide**
  (the floor's own principal axis is NOT the groove axis), origin on the peptide centroid. So maps
  compare without prealigning inputs, and `force_pca`/canonicalisation are not prerequisites.
- Heights come from **ray casting in the groove frame**, not Shrake-Rupley points: sphere sampling is
  fixed in global axes, so the map wobbled 1.35 Å (median cell) under a rigid rotation of the input.
- Validated on all 374 Canonical2026: `frac_above_ridge` 0.054 (8-mers) → 0.569 (13-mers); same-epitope
  maps closer than cross-epitope at P = 0.917.
- **Both faces, one frame.** `surface_map(structure, side="tcr")` rasters the receptor's V-domain
  *underside* instead of the groove — same frame, same grid, lowest point per cell instead of
  highest — and `surface_complementarity(pmhc_map, tcr_map)` compares them cell for cell:
  `shape_r` (positive = complementary), `charge_r`/`charge_product` (**negative** = complementary,
  plus meeting minus), `phobic_r`/`phobic_product` (positive), `gap_mean`/`gap_sd`, and the
  per-channel map distances `d_h`/`d_charge`/`d_phobic`. CLI: `--side tcr`, `--complementarity <csv>`.
- **The two limits on that comparison are calibrated, not assumed** (60 Native2026 crystals).
  `DEFAULT_EXTENT` is sized for a class-II 15-mer, far wider than any receptor footprint: over the
  full extent a TCR projection reaches only **0.741** of occupied pMHC cells *at any Z cutoff*, and
  the shortfall is nearly all at the far groove end (0.348 beyond y = +15 Å against 0.987 near
  y = 0). So `COMPARE_WINDOW = (12, 12)` Å crops first, and `MAX_GAP = 10` Å is where the Z curve
  has gone flat (0.895 at 4 Å → 0.951 at 10 Å → 0.962 uncapped). Peptide-owned cells are covered at
  **0.917** even over the full extent, so cropping discards no epitope surface.
- **The gap is mostly negative and that is correct**: median `h_tcr − h_pmhc` = −1.7 Å with 71% of
  cells interdigitated, because the two faces interlock rather than stack. The Z cutoff is one-sided.
- Figures: `viz.surface2d.render_surface_map(smap, channel)` → SVG string (hand-built, zero deps).

## Parsing speed — the PDB fast path

- `parse_structure` takes a **vectorised ATOM-record path** for `.pdb` at `model=0`: the records are
  sliced as one uint8 array and residue breaks found with `np.flatnonzero`. **3.3x** end to end
  (22.1 -> 6.6 ms per crystal), and it is what makes a dataset-scale pass tractable — Biopython was
  86% of the wall clock (19.0 s of 30.3 s over Native2026, building 2.1 M Atom objects to discard).
- **Exact or it does not run.** `_parse_pdb_fast` returns `None` — falling back to Biopython — on a
  blank element column, a short ATOM line, or an unparseable coordinate. Coordinates are read
  **through float32**, matching Biopython bit for bit rather than being more accurate: a PDB's three
  decimals come back as `59.42599869`, and every previously computed number depends on that.
  Verified by re-deriving all three potentials of `pot_realmock.py` from structures: max |delta| =
  0.0000000000 over 1,140 cells.
- The remaining cost is ~655 k `Atom` constructions per 120 structures. That is the data model, not
  the language — **a C++ extension would not help without making residues array-backed**, which is a
  redesign, not an optimisation.

## Rotamer-averaged contacts — `tcren.rotamers`

- `contact_probabilities(structure, interface) -> df[..., p]` and `soft_energy(structure, potential)`.
  Replaces the 0/1 contact indicator with `p_ij` averaged over χ rotamers, Boltzmann-weighted under
  **DOPE** (never under the potential being scored — that would be circular).
- Why: under a deliberately wrong χ1 the hard contact set keeps Jaccard 0.66 and the energy moves by
  |ΔΦ| = 0.524; the averaged map keeps 0.95 and moves 0.054, against energies of magnitude 0.4–2.2.
- Rotating everything past Cβ about Cα–Cβ **is** χ1 exactly. `max_chi=2` default (3^n rotamers).
  ~0.24 s/structure — see `refine/CPP_REWRITE.md` for when this needs to be C++ (MC loops, not one-shot).

## Backbone dynamics / peptide conformational stability — `tcren.dynamics` (the Sewell hypothesis)

- **Library only — there is no `tcren stability` command.** Call `peptide_stability` /
  `stability_table`, both exported from the top-level `tcren` namespace. Batch driver, with the CPL
  design and the AUC bookkeeping: `scripts/sewell_stability.py`.
- `peptide_stability(structure, intra_weight=1.0) -> Stability`: flexible-backbone Metropolis MC of
  the peptide's φ/ψ against DOPE. The readout is **not a better pose** — it is `rmsf` (ensemble
  spread), `drift` (how far the mean moves from the input) and `energy_gap`.
- **Why it exists**: a contact potential scores whichever conformation it is handed and cannot tell a
  peptide held in the TCR-facing conformation by its own side chains from one merely drawn there.
  Sewell (`suggestions/sewell.txt`): intra-peptide interactions stabilise the productive bulge, and
  "poor binders could perhaps still make many contacts but fail to stabilise" it.
- **`intra_weight` is the switch.** `1.0` includes the peptide's DOPE contacts with itself, `0.0`
  removes them; `stability_table` runs both at the same seed and reports `delta_rmsf`.
- Moves are torsional and exact — perturb one φ/ψ, rotate everything downstream. Both peptide termini
  are free so no loop closure is needed. `backbone_torsions` is the tree; **φ moves the side chain**
  (Cβ hangs off Cα, off the N–Cα axis), **ψ does not** (Cβ is on the Cα side of Cα–C).
- **Noise**: `rmsf` is the least noisy readout (CV 0.115 over seeds), `drift` 0.25, `energy_gap` 1.56
  — the last is a max over the trajectory, so do not use it for per-structure comparisons. Aggregate
  over structures; ~0.9 s per run at 4000 steps.
- **Not MD**: no solvent, no force field, no time. DOPE is knowledge-based, so `rmsf` compares
  *between* structures run with the same settings, not against an MD RMSF in Å.

## Side-chain repack — `tcren.repack` (native `_relax.repack`)

- `repack(structure, chains=("PEPTIDE",), max_chi=2) -> (structure, report)`. Places each side chain
  in the χ conformer DOPE prefers; `report` carries `n_conformers`, `energy`, `p_best` per residue.
- **Same input, same atoms, measured**: side-chain RMSD 4.131 → **2.364 Å in 6 ms**; OpenMM returns
  4.133 Å (unchanged) in 3103 ms, because a local minimiser cannot cross a torsional barrier.
  8/8 structures improved. `tcren refine --repack` runs it after the rigid-body MC.
- **It rotates side chains a model has; it does not build missing ones.** `substitute_peptide` strips
  past Cβ, so that path still returns 44 of 77 heavy atoms — side-chain *construction* is a separate,
  open roadmap row (`refine/CPP_REWRITE.md`).
- Conventions the kernel relies on: the input conformer is index 0 of every enumeration (so each local
  frozen-environment objective cannot worsen; this is not a guarantee for the assembled pose), and weighting is **mean field** (each residue against its
  neighbours at their input conformation — coupled side chains are not resolved).

## Peptide position — `tcren.scoring`

- `peptide_positions(cm, structure)` adds `peptide.pos` (1-based), `peptide.aa`, `peptide.role`
  (anchor/tcr_facing). **Pass the structure**: the sequence comes from it, not from the contacts —
  the class-II 9-mer register heuristic needs the whole peptide (4ozg's gliadin core resolves to
  P1/P4/P6/P9 = positions 2/5/7/10, the published register).
- `position_weights(ann, "uniform"|"central"|"tcr_facing")` → `score_peptides(..., weights=)`.
  `position_profile` decomposes Φ along the peptide (sums exactly to the total); `central_strain`
  isolates the middle band.
- All scoring reweighting goes through one hook: `pipeline._contact_weights`.

## Wrong-TCR decoys — `tcren.shuffle` (`tcren shuffle`)

- `make_decoys(structures, n_per=10, within_class=True, seed=0)` / `graft_tcr(pmhc_source, tcr_source)` /
  `run_shuffle(dir, out, n=10, ...)`. Keep each **oriented** complex's pMHC intact, graft on a **different**
  complex's TCR (within-MHC-class derangement, no real pairing) → wrong-TCR-on-real-pMHC negatives. Real
  (label 1) vs decoy (label 0) trains a **label-free** TCR-recognition classifier; peptide:MHC energy is
  invariant under the graft (built-in control), TCR:peptide/TCR:MHC contacts are new.
- **Direct chain replacement, NOT `orient.substitute_tcr`.** Inputs must be canonically oriented (one common
  frame); the graft copies chains with no per-pair alignment, so each grafted TCR keeps its native MHC–TCR
  docking angle → the decoy set spans the real docking-angle variance (substitute_tcr would collapse every
  donor onto the host's MHC pose). CLI: `tcren orient -s natives/ -o oriented/ && tcren shuffle -s oriented/ -o shuffled/ --n 10`.
- **Finding (2026-07-05):** real-vs-shuffled is learnable at **AUC 0.876** (RF; `Phi_tcr_pep` the top feature).
  But a shuffled-trained (crystal) model does **not** transfer to AF-modeled TCRvdb (0.55–0.62 vs AF 0.79) —
  crystal→AF distribution shift. Use it as a label-free recognition prior / supplementary benchmark, not as a
  drop-in TCRvdb scorer.
- **The two shuffle-trained recognizers are gone (2.26.0), with their shipped weights.**
  `GaussianBNClassifier`, `BayesianLogisticRecognizer`, `encode_features`, `frozen_recognizers`,
  `real_probability`, the `p_real` / `p_real_bn` columns and
  `src/tcren/data/shuffle_{bn,logistic}.json.gz` were removed in one release: their coefficients
  were frozen against training sets nobody could reconstruct. **What they measured stands and is
  the reason to keep the decoy set**: real-vs-shuffled 5-fold CV ROC-AUC **0.885** for the
  distribution-aware logistic, 0.865 BN, 0.870 raw logistic, and on TCRvdb a *supervised* refit
  with the same encoding reads **0.860** pooled against TCRmodel2 confidence's 0.794. The frozen
  crystal-trained transfer does **not** carry to AF-modelled TCRvdb (0.53). Decoys are regenerable
  (`tcren shuffle --seed 0 --n 10`); the manifest is committed, the full PDBs live on HF (351 MB).
- **`tcren recognize` / `recognition_features` (2026-07-06):** `recognition.recognition_features(struct)`
  ports the manuscript's descriptor extractor into tcren (docking geometry + TCRen/MJ Φ & poly-Ala
  ΔΦ + contact tallies + biopython ΔSASA `burial` + `mhc_class_bin`) — verified **byte-exact** vs
  `canonical2026_features.csv` (burial max diff 4e-11). Uses `import_structure` (C-gene trimmed);
  **no `_geom` C-ext needed** (only arda for annotation). CLI `tcren recognize -s pdbs/ -o out.tsv`
  writes one TSV row per PDB, **descriptors only — no fitted composite**; `--features` instead
  turns an existing `tcren features` table into `Q`, `T` and `S`. The `--scores` and
  `--features-only` flags went at 2.26.0 with the models they switched. koff joins the same rows
  under `--mechanics` (2026-07-28), and only the alanine ΔΔG, which is per-residue rather than
  per-structure, stays its own command `tcren ddg`.
- **`--full` feature table (2026-07-13, audited 2026-07-28):** `recognition_features(struct, full=True)` /
  `tcren recognize --full` append the **18 CDR3-frame** descriptors
  (`cdr3{a,b}_{reach,ou,ow,on,au,aw,an,topep,ext}`, FramePose groove-frame projection — the `cdr3b_*`
  strain signal). Tuples: `RECOGNITION_FEATURES` (**40**), `CDR3_FRAME_FEATURES` (18),
  `FULL_FEATURES` (**58**). The 34 / 52 recorded here were the pre-2.26.0 counts.
- **Descriptor audit (2026-07-28):** every energy column is `Phi_*` (`e_cdr12`/`e_cdr3a`/`e_cdr3b` →
  `Phi_cdr12`/`Phi_cdr3a`/`Phi_cdr3b`); the duplicate `e_tcr_mhc` and `ct_tp_hydrogen_bond` columns
  are gone (they equalled `Phi_tcr_mhc` and `n_hbond`); the **12 matrix-swap** columns
  (`{tcren,mj,d}_{tp,cdr12,cdr3a,cdr3b}`) were **removed** — `tcren_*` duplicated `Phi_*`, and MJ is
  not the potential used on TCR:peptide. New: `crossing_signed` (signed scanning angle, carries
  docking polarity) and `DESCRIPTORS` / `descriptors(family, tcr_only=)` — the catalogue giving each
  column's family (`placement`/`interface`/`topology`/`energetics`/`potts`/`kinetics`, with
  `geometry` and `physics` surviving as aliases) and whether the receptor enters it. **The `score`
  family and the `score` invariance class went at 2.26.0 — the catalogue is descriptors only**, and
  `descriptors("score")` now raises. Five columns do not involve the receptor —
  `Phi_pep_mhc`, `dPhi_pep_mhc`, `mhc_class_bin`, `Phi_pep_int` and `n_pep_int`; they carry cohort
  identity, so receptor questions must use `tcr_only=True`.
- **`-t/--threads` on `tcren scoring` and `tcren recognize` (2026-07-26):** both accept a file, a
  directory, a `.tar.gz`, a quoted glob or a `.txt` manifest; `-t N` runs N concurrent workers (`-t 0`
  = all cores). The cohort-relative `Q` is still computed over the **whole** set, never per batch.
  `scoring` gains ~7.6x on 8 threads; `recognize` less (its cost is Python
  featurisation, not mmseqs), so batch its annotation rather than expecting linear scaling.

- **`cohort.q_coupled` / `cohort.coupling` — DEPRECATED at 2.12 (2026-07-26).** Both survive as
  functions; what replaced them is the score set's `binder_score`, whose weights between geometry
  and energy come from one covariance rather than from a cohort correlation:
  `q_coupled(q, energy)` = `¼[1+erf(z(Q)/√2)]·[1+erf(r·z(ΔΦ)/√2)]` with `r = coupling(q, energy)`, the
  cohort correlation between the geometry and energy channels. Two Gaussian tail probabilities multiplied
  — binding needs both an interface and favourable residues in it — with the energy admitted in
  proportion to `r`, because `E[z(Q)|z(ΔΦ)] = r·z(ΔΦ)` is exactly the part of it that is evidence about
  interface nativeness. **No threshold, no softness constant, no fitted coefficient.** A forced cohort
  disarms itself: `r < 0` flips the energy's sign, `r ≈ 0` collapses its factor to the constant ½ and
  leaves geometry alone. `coupling()` alone is a label-free forced-pose diagnostic (TCRvdb: −0.25 on the
  template-forced GLCTLVAML poses, +0.48 on the clean YLQPRTFLL ones). On TCRvdb it reaches macro
  ROC 0.799 / PR 0.817 / P@10%recall 0.949, ahead of every TCRmodel2 confidence with no generative term.
  For receptor ranking pass the **TCR**-referenced ΔΦ (the peptide is fixed, so the peptide reference
  carries nothing); for peptide ranking pass `reference_delta`'s peptide-referenced ΔΦ.

- **`cohort.phi_score` — the contact-energy channel (2026-07-24, renamed from `f_score` at 2.26.0):**
  `phi_score(table)` = `z(-(Phi_tcr_pep + Phi_tcr_mhc))`, binder-oriented (`cohort.PHI_TERMS`, the
  default `terms=` argument). Φ reads contact chemistry but is **pose-conditional**: it inverts on
  forced poses, which is the whole reason the shape channel exists. Do not hand-combine
  `phi_score` with `Q` — the score set's `energetics` channel is what weights them, and it does so
  from the covariance rather than from a choice.

## MHC allele reference — `tcren build-mhc-ref`, built on demand

- **Not bundled in the wheel.** The curated allele FASTA is built from IMGT/HLA + UniProt mouse by
  `tcren build-mhc-ref`, which must be run **once after a `pip install`**; every command that
  annotates a structure needs it. It lands in `database/mhc/` under `paths.tcren_home()`.
- `paths.tcren_home()` is the single root for on-disk reference data: `$TCREN_HOME` when set,
  else the source checkout (recognised by its `pyproject.toml`), else `$XDG_CACHE_HOME/tcren`
  (`~/.cache/tcren`). Added at 2.12.1, because deriving the roots from `Path(__file__).parents[2]`
  resolved to `site-packages`' parent in an installed wheel and every annotate failed.
  `paths.data_dir()` is `data/` under it, overridable on its own with `$TCREN_DATA_DIR`.

## MHC mapping speed — `mhc.reference.reference_db()`

- `easy_search(query, reference_fasta())` rebuilt the 28k-allele target DB + k-mer index on EVERY
  call (~4.5 s). `reference_db()` runs `createdb` + **`createindex`** once into gitignored
  `data/mhc_cache` (the index is the real cost); `map_mhc` / `annotate_mhc_batch` pass that DB →
  repeated single-structure searches drop to ~0.9 s (5×). Rebuilds if the FASTA is newer.

## Annotation CLI — one `annotate`, no separate `mhc`

- `tcren annotate -s … [--regions all|tcr|mhc|peptide] [--pseudo]` emits ONE per-residue markup
  covering TCR (CDR/FR), MHC groove (HELIX/FLOOR) and peptide. `--regions` filters by chain class;
  `--pseudo` adds `MPS` rows. The old `tcren mhc` command was removed — its allele/class info is
  available in the library via `mhc.map_mhc` / `mhc.annotate_mhc`.

## MHC pseudosequence (MPS) — `tcren.mhc.annotate_pseudo`

- Marks the NetMHCpan 34-residue groove pseudosequence on an annotated structure (region `MPS`).
  Committed FASTAs `src/tcren/data/{mhci,mhcii}_pseudo.fa` (built by `scripts/build_pseudo_fasta.py`
  from NetMHCpan tables; unique seqs, header `<allele>|n=<count>`).
- The 34 positions are **scattered**, so mmseqs/local search can't find them (no shared k-mer).
  Instead each candidate 34-mer is threaded through the chain with a **fitting alignment** (free
  chain gaps; positions are N→C ordered) — ~0.1 s over all ~4k seqs, no prebuilt index. One best
  hit is chosen; class I marks MHCa only (never β2m), class II splits across MHCa+MHCb.
- Validation notebook `notebooks/mhc_pseudosequence_mps.ipynb`: MPS residues vs. 5 Å peptide
  contacts (~half are direct contacts; the rest line the groove toward the TCR).

## Structure output format — `--mmCIF` / `--compress`

- Every command that writes a structure (`orient`, `superimpose`) outputs plain `.pdb` by
  default; `--mmCIF` switches to `.cif`, `--compress` adds a trailing `.gz`. In the library:
  `structure_output_path(dir, id, mmcif=…, compress=…)` + `write_structure(s, path)` (dispatches
  PDB/mmCIF by suffix; a minimal `_atom_site` mmCIF loop that round-trips through Biopython).

## Fetching recent structures — `tcren fetch-recent` / `tcren.recent`

- `tcren fetch-recent [--discover --after YYYY-MM-DD]` → `data/pdb_recent/` (gitignored):
  downloads PDB ids (Native2026 seed; `--discover` adds an RCSB full-text TCR:pMHC search) as
  **mmCIF `.cif.gz`** (the PDB deprecates split `.pdb`; handles **extended >4-char ids**), then
  keeps only complexes with all **5 required chains** (MHCa + b2m/MHCb + peptide + TCR pair),
  validated by one batched annotation pass. `tcren.recent.{fetch_ids,discover_similar,native2026_ids}`.

## Paper-reproduction module (`tcren.paper`)

```python
from tcren.paper import (
    bootstrap, fetch_hf_structures, fetch_vdjdb, fetch_pdb_dates,
    copy_external_inputs, copy_legacy_results,
    contact_table,            # mir extract_contact_map replacement (per structure)
    annotate_structure_set,   # batched TCR annotation over a folder -> (contacts, markup)
    mhc_annotation,           # per-structure MHC allele + class over a folder
    compare,
)
```

- Notebooks live in `notebooks/natcompsci2022/`. HF structure sets are fetched (gitignored) into
  per-set folders directly under `notebooks/data/` — **`notebooks/data/Native2022`,
  `notebooks/data/Native2026`, `notebooks/data/PolyV2022`, `notebooks/data/Bobisse`,
  `notebooks/data/Bigot`** (no `structures/` wrapper). `Canonical2026` is Native2026 after
  `tcren orient`. All structures are gzipped (`*.pdb.gz`).
- Non-structure inputs + 2022 comparison baselines are **committed** under
  `notebooks/natcompsci2022/data_legacy/` (vdjdb, Birnbaum, MJ/Keskin, IEDB, epitope lists,
  `TCRpMHCmodels.tar.gz`, PDB dates, mir/R oracle) — never a pipeline input. `results_new/` is computed.
- Root `data/` — resolved through `paths.data_dir()`, i.e. `data/` under `paths.tcren_home()`,
  which is the checkout in a dev install and `~/.cache/tcren` from a wheel — holds the library
  dataset (gitignored structures): `Native2026` (orientation
  references), `Canonical2026` (the default `superimpose` database), `PDB_date.tsv`,
  `TCRen_potential.csv`. `Canonical2026`'s `orient_metadata.json` moved into
  `src/tcren/data/` — it must ship in the wheel, since `fetch-data` downloads only structures and
  an installed user has no repo `data/`. `setup.sh` runs `tcren fetch-data` at install to
  populate `Native2026` + `Canonical2026` from HF (or lazily on first `superimpose`/`orient`).
  Orientation references load 1ao7/1fyt from `data/Native2026` via `tcren.paths`. The numerical
  regression oracle (legacy mir/R outputs: `contact_maps_PDB.csv`, `tcren_am/tcren.txt`, the
  `example/` set) lives under `tests/assets/oracle/`; the legacy R/Java pipeline was deleted.

## Directed binder-to-target contact models

Residue counts, potentials, backgrounds and contextual likelihoods retain ordered
binder → target roles. Derivation never adds reverse residue counts or averages
a potential with its transpose. TCRen's receptor and target marginals are distinct;
published MJ tables retain their original entries without altering TCRen to match
their intrinsic symmetry. Only individual atom-pair chemistry may use exchangeable
whole amino-acid:atom tuples. Angular keys retain binder/target residue order and
the endpoint whose theta/phi is measured. New resources use contact-marginals-2;
version-1 angular resources are rejected by active scoring.

Potential.scale() and offset() use every finite defined directed cell. Frozen
historical interface moments are unchanged. An intra-peptide term requires an
intrinsically symmetric source potential; it rejects asymmetric TCRen instead
of averaging its entries.

## Geometry: contacts, region pairs, docking angle

```python
from tcren.project2d import region_pair_contacts, region_pair_summary  # needs chain-typed + MHC-annotated structure
from tcren.docking import docking_angles
```

- `region_pair_summary(s, kind="closest"|"cb"|"ca")` — inter-chain contact counts for **every**
  region pair (CDR↔peptide, peptide↔MHC, TCR↔MHC, intra-TCR …), not just one interface. Three
  contact definitions: `closest` (5 Å closest heavy-atom pair — the original TCRen definition, the
  only kind that carries a `contact_type` bond classification), `cb` (8 Å Cβ, Cα fallback for Gly),
  `ca` (12 Å Cα). Region-pair labels are ordered canonically (direction-independent).
- Bond types come from the heavy-atom heuristic `project2d.classify_contact` (salt_bridge /
  hydrogen_bond / aromatic / hydrophobic / polar). The external `biotite.structure.hbond`
  (Baker-Hubbard) needs **explicit hydrogens** — it returns 0 on X-ray crystals (no H), so it is
  only useful on protonated / NMR / MD structures. Use the heuristic for crystal structures.
- `docking_angles(s)` — TCR crossing + incident angle from a groove frame built from the peptide
  principal axis + peptide→TCR normal (NOT the whole-complex PCA basis, which the Vα–Vβ spread
  contaminates). ~20–70° crossing for αβ; requires a peptide chain (γδ without peptide raises).

## Gotchas

- nbconvert: pass `--ExecutePreprocessor.kernel_name=python3` (or `=tcren-nb`) or cells silently don't run.
- MHC allele strings from the mapper carry full resolution (e.g. `HLA-A*02:608N`); for
  IEDB-style matching, truncate to 2-field group (`HLA-A*02`).

### MHC sequence ambiguity and assay joins

`MhcCall.allele` is a deterministic lexicographic representative of co-best returned
MMseqs hits, not a resolved genotype. `allele_candidates` preserves their distinct
alleles (bounded by the hit limit). `scoring` exports a semicolon-separated
`allele_candidates`; `annotation.batch.mhc_annotation` uses `mhc.allele_candidates`.
A cropped groove cannot distinguish alleles differing outside its sequence, or infer
null expression from a tied hit. Keep the experimental restriction separate. Validate
coordinate peptide, TCR junctions and MHC evidence before attaching experimental
labels; a matching filename/hash is insufficient. ARDA owns TCR FR/CDR coordinates;
`annotate_mhc_batch` owns MHC class/role and canonical groove-region projection.

## Deterministic literature-topology reductions

`features -i topology` / `topology.literature` sorts physical residue keys before
participation means and target-chain/CDR-loop keys before contact-order means.
The definitions, length normalisation and missingness remain unchanged. Class II
keeps each groove chain's sequence axis. Input edge order must not alter the four
`co_pep/co_mhc/partcoef_tcr/partcoef_pmhc` values in one numerical environment.
No rounding/tie-merging; recalculated ranks can differ from historical unordered
means. Retain all matched metric losses alongside gains. See `docs/features.rst`.

Spatial co-localization: `export --quantity pmhc-cells -s MANIFEST -t N -o cells.tsv`
uses the existing batched annotation/region APIs and receptor-independent groove
frame, then writes exact occupied grid cells and `cells.tsv.summary.tsv`.
`topology.pmhc.spatial_features(cells)` reduces those exports with vectorized
Polars: population height/KD and positive ridge-excess/KD covariances, marginal
means, E(excess*KD), common counts and conditional above-ridge chemistry. Existing
90th-percentile ridge and smoothed KD map defaults are retained; no fitted weight
or new threshold. Covariance units are Å times KD index; missing ridge remains
undefined and no-above conditional hydropathy is missing. Positive covariance
means relatively more hydrophobic values occur at higher cells, not affinity.
Use the summary with `associate` on original matched score/ipTM rows and separate
class/template strata. Scalar or spatial association cannot alone identify a
forced pose, apo strain, or loop-formation kinetics. See docs/associations.rst.

External pose truth must verify physical chain roles, exact peptide and shared
observed receptor/MHC sequence spans, construct substitutions, reference-coordinate
hashes and per-interface definitions. Do not treat a six-interface DockQ average
(including TCR internal and peptide:MHC) as a pure receptor-placement metric.
Keep recognition interfaces separate; a noncognate or substituted reference is a
reference-deviation diagnostic, not structural truth for the candidate's own
complex. Retain all exclusions and matched per-feature/ipTM rows, class/template
and clone separation. A singleton reference check supplies no detector accuracy.
External DockQ 2.1.3 requires an isolated NumPy<2 environment; do not alter the
locked TCRen wheel to accommodate it. Benchmark owns this validation, with the
library export unchanged. See audit/2026-10-16/CPL_POSE_TRUTH.md.

### Cognate saved-pose validation

Reuse `features -i placement,interface,topology --no-metadata -s MANIFEST -t N`
for receptor-relative diagnostics and `export --quantity pmhc-cells` for pMHC
context; each pass annotates its entire set in one batch. All currently needed
quantities are already exposed. Keep original per-rank overall ipTM and
TCR-pMHC-specific ipTM separate: TCRmodel2 documents a two-partner calculation,
so a generic internal-interface bias cannot be assumed for its specific field.
`clash_score` includes peptide overlaps with both receptor and MHC, in summed Å;
it is not receptor-only. Better reference agreement can coexist with worse
side-chain packing. Compare raw signed descriptor associations on the same finite
feature/truth/both-confidence rows; do not choose signs or fit a replacement.
Saved ranks of one complex are correlated models, not independent target samples.
Exact-reference template inclusion and broad template-coverage strata are distinct
metadata. See `audit/2026-10-16/COGNATE_POSE_VALIDATION.md`.

### Structural class-II core export

Use `export --register-reference-id 1fyt --register-p1 3` with the verified reference
included in the same input batch. `.registers.tsv` ranks every complete forward/reverse
9-residue core by Cα RMSD after **MHC-only** groove alignment; exports physical
P1/P4/P6/P9 identities, runner-up margin and groove/core residuals in Å (plus the exact rigid transform). No motif,
fit, threshold, or extra annotation search. Inspect source evidence before promoting
register minima. `anchor_indices` is still the historical sequence prediction.
`--peptide-selection TSV` supplies source-verified `pdb.id, chain.id, first_pdb,
last_pdb` (tab-separated, inclusive author numbers) for engineered fused antigens;
auxiliary peptide chains become unassigned. Coordinates/indices stay unchanged.
Selection and core exports are hashed. API: `refine.structural_register` module.

### Physical roles before DockQ correspondence

Unrestricted external DockQ with permissive sequence mismatches can assign a short
peptide to an MHC chain and exchange their biological roles. A reproduced scalar
then authenticates numerical correspondence only. For pose validation, constrain
TRA/TRB/peptide/MHC roles explicitly and retain observed-span/construct checks,
physical copies, per-interface values, original totals and missingness separately.
5NMD is unbound 868 (no pMHC); 4MNQ is affinity-enhanced ILA1, not SB27.
Neither supplies identical-complex SB27 truth. Existing batched CLI descriptors
remain sufficient; do not add a fitted correction for a reference-mapping error.

Permissive historical mappings can also combine crystallographic copies and
change the number of nonempty interfaces in the total. Retain that count; do not
assume every historical multimer total averages the same six physical interfaces.


## Pose accuracy versus forced deformation (17 October 2026)

A forced-pose hypothesis concerns implausible peptide/CDR deformation accommodating a non-binding pair. Receptor-relative placement, DockQ disagreement and score/confidence inversion are distinct measurements and cannot label that hypothesis as truth. Existing `md quality`, coordinate diagnostics and local frame descriptors are structural diagnostics; no detector accuracy follows without independently verified outcomes and appropriate construct correspondence. Natural induced fit can rearrange loops in genuine binders. Compare physical biological roles and observed receptor/peptide/MHC constructs before accepting pose truth; keep saved ranks within one target rather than counting them as independent targets. Fit lines used to display cohort associations do not become fitted scoring models.

## Separate contact topology and effective rank

`tcren export -s structures.txt --quantity contact-topology -o graphs.tsv -t 1`
uses the existing batched annotation and ContactMap. It returns three rows per
complex: CDR:peptide, CDR:MHC and peptide:MHC (both groove chains for class II).
Vertices are contacted physical residues; edges are unique residue pairs at the
specified heavy-atom cutoff. `beta0` counts graph components and
`beta1 = edges - vertices + beta0` counts independent graph cycles. Empty graphs
have zero counts and `available=false`; they are not zero-valued binder scores.
These are not the spatial flag-complex Betti numbers in `footprint`.

`m_erank_tp` and `m_erank_tm` are different quantities: singular-value diversity
of CDR C-alpha proximity to peptide and MHC helices, respectively. With
`K_uv=exp(-d_uv²/(2 sigma²))`, sigma=5 Å, singular values s and
q=min(number of resolved residues on each side), the normalized order-2 Hill
rank is `(sum(s))² / (q sum(s²))`. Its ideal range is 1/q to 1; underspecified
matrices return NaN. It is not a fitted score, contact count or graph Betti number.
Docking angles, CDR3 apex/reach and peptide protrusion remain distinct descriptors;
none alone validates unnatural deformation or nonbinding.

`tcren export -s structures.txt --quantity backbone -o backbone.tsv -t 1`
exports intrinsic C-alpha bond lengths (Å), pseudo bond angles and signed
pseudo torsions (radians) separately for peptide and each annotated CDR1/2/3
on TRA/TRB. It reuses the geometry primitives underlying `cdr3_internal_coords`.
Every window retains observed indices, author residue numbers, insertion codes
and amino acids, plus the export's complete residue/annotation sidecars.
Missing C-alpha/observed-index gaps are not bridged; these are observed-coordinate
windows, not proof that the chain's missing SEQRES residues were resolved.
Correspondence must be audited before comparing windows; torsion differences
must be circular. These raw measurements have no strain threshold, fitted
reference density or nonbinding label. Natural induced fit remains possible.

## Local topology/topography of individual contacts

`export --quantity contact-context --potential tcren2` retains the complete
TCR:peptide contact table and adds all-TCR bipartite endpoint degrees, four-cycle
participation, closest peptide atom height, projected receptor-free pMHC surface
ownership/height/KD hydropathy, excess above the existing helix rim and exact
`pair_energy`. Physical identities/region markup and existing batch/provenance
sidecars apply. Graph support and surface protrusion are distinct quantities;
this is explicitly not chemical `--drop-untyped` filtering. No mask is silently
applied. Missing ridge/out-of-grid heights remain undefined. Zero retained
contacts is distinct from missing topography. Use fixed, declared masks and
matched score/ipTM rows; do not tune cutoffs or score directions on assay labels.
See `docs/fit-free.rst` for definitions and the projection limitation.

`contact-summary -i contacts.tsv -o contacts.summary.tsv` owns table-only canonical
pair-count reductions for all/cycle/above-ridge/joint contact selections. It checks
the stamped source and complete accounting, keeps zero-contact models, exposes
missing entries and undefined surface selections, and records input hashes. No
structures, MMseqs searches, label inputs, tolerance rounding or fitting. Use these
CLI sums downstream; R owns tidy summaries/macros rather than replacement scoring.

## Atom-resolved context exposure

`export -s MANIFEST --quantity contact-exposure -o OUTPUT` exports all heavy atoms
of annotated alpha-beta TCR and intact pMHC, including noncontacts, without a
potential. `isolated_sasa_A2`, `bound_sasa_A2`, `partner_burial_A2` retain physical
atom keys and binder/target roles. Target isolation retains MHC/B2M; peptide-alone
area and `mhc_shielding_A2` are separate peptide-only fields. Residue areas are
sums over atoms. Uses existing Bio.PDB Shrake-Rupley, 100 points, 1.4 A probe,
supplied-coordinate frame. Coordinate stripping is not experimental apo geometry.
This is an exposure export, not a context-normalized score or water free energy.
Use bounded manifest chunks: the general export collects its output tables.

`export --quantity surface-contacts` uses the native SASA sphere and emits
binder→target atom pairs with separate endpoint allocated areas; the companion
`.atoms.parquet` retains exposed noncontacts. Area is conserved per endpoint;
nearest expanded-sphere signed gap owns each patch, exact ties split equally.
`contact-penalties --surface-context` adds raw/normalized area-weighted sums
using isolated CDR and intact-pMHC peptide exposed compositions. Original
scores remain unchanged. `contact-compare` compares both new sums on the same
finite rows; the nullable per-area diagnostic is not a comparison arm. See
`docs/contact-penalties.rst` for equations, missing-support accounting and units.


`contact-penalties --surface-context --surface-ownership accessible_boundary`
selects nearest sampled isolated opposing boundary ownership. The default
`deepest_sphere` remains unchanged. This option conserves each endpoint's area
and assigns only to atoms with positive sampled isolated exposure; exact ties
split among distinct atom owners. Python API: `contact_surface(structure,
ownership="accessible_boundary")`. The native sampler's optional `samples=True`
exports exposed sample coordinates, atom indices, areas and lost/bound flags.
SciPy batched nearest-neighbour queries use one worker. Record the chosen method,
keep runs separate, and do not claim an accuracy gain from geometric invariants.

## Frozen candidate-pool scoring

`score -s templates.txt -c candidates.tsv -p karnaukhov2022|tcren2 -t 0 -o scores.csv`
uses parent-only batch annotation followed by coarse template workers, sharing
one potential/options tuple per chunk. Default0 uses the allocation;1 is serial.
Candidate TSV uses `complex.id,peptide`; metadata/labels are benchmark-owned.
Preserve the exact candidate pool and duplicate multiplicities when comparing
potentials. `rank` generates another random background and is not a substitute.
Sequence threading changes residue identities on template contacts; it does not
provide physical mutant sidechain atoms for the atom-count/surface correction.

`score-receptors -s templates.txt -c receptors.tsv -p karnaukhov2022 -t 0 -o scores.csv`
threads paired anchor-excluded CDR3s onto fixed CDR3:peptide contacts using
`score_peptides(..., substituted_side="from", require_same_length=False)`.
TSV columns: `complex.id,candidate.id,cdr3a,cdr3b`; IDs unique per template,
distinct IDs with identical sequences retained. Standard uppercase amino acids,
exact observed template CDR3 lengths, one physical TRA/TRB CDR3 and exactly one
observed standard-amino-acid peptide are required. Cutoff must be finite and
positive; missing templates fail. One fresh ContactMap per template, parent batched annotation,
coarse bounded workers and serial native kernels; `-t 0` uses allocation.
Output `score.tra,score.trb,score` retains directed potential sums and original
missing-cell omission. Default CDR3 mode excludes framework/CDR1/2 contacts.
No coordinates, anchor synthesis, peptide-only weighting/reference options or fitting.

`--receptor-region junction` requires `junctiona,junctionb` instead of `cdr3a,cdr3b`.
Junctions include the actual observed Cys104 and F/W118: recover FR3 last/FR4 first,
requiring unambiguous annotations and immediate adjacency to contiguous observed
CDR3. Missing/gapped/invalid anchors fail. Candidate junctions start C, end F/W and
match the exact observed junction length. Include physical anchor contacts and
index the filtered view from junction start; never mutate structures/contact maps.
Default `--receptor-region cdr3` behavior and columns remain unchanged.
