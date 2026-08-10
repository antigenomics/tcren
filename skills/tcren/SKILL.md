---
name: tcren
description: tcren — TCR-pMHC contact potential (TCRen) pipeline; conventions and public API
---

# tcren Skills Guide

`tcren` reproduces and extends the TCRen contact-energy potential (Nat Comput Sci 2022)
on a pure-Python pipeline (structure parsing → contacts → TCR/MHC annotation → potential
derivation → epitope-ranking benchmarks). Annotation uses the `arda` package
(mmseqs2-backed), a runtime dependency published to PyPI as `arda-mapper` (imports as `arda`,
`>=2.5.7`) — no separate checkout and **no `ARDA_HOME`** (arda auto-fetches its reference and a
static mmseqs binary on first use). Repo-local `.venv` via uv, no conda (`bash setup.sh`).

## Batch annotation — never loop (mmseqs2 is the parallel layer)

**All structure annotation (TCR chain typing AND MHC allele mapping) must gather every
sequence first, make ONE batched mmseqs2 call, then map the output back for downstream
per-structure analysis.** mmseqs2 parallelises internally across threads — that is the
parallel layer; Python orchestration is a single call.

- Each per-structure annotate call pays a fixed ~825ms mmseqs2 process+index-load cost;
  a batch of 300 sequences costs the same ~930ms total.
- A `ProcessPoolExecutor(fork)` over structures **deadlocks** (fork after mmseqs2/BLAS
  spawn threads). A `ThreadPoolExecutor` runs but still pays the fixed cost N times.
- `paper/helpers.py::_batch_annotate` does TCR annotation for a whole dataset in 2 arda
  calls (human + mouse). MHC annotation uses the same pattern: `mhc.annotate_mhc_batch(structures)`
  — ONE mmseqs search over every candidate MHC chain, sliced back per structure.

Reference: `arda.annotate_sequences([(id, seq), ...])` — one call, threads internally.

## Threading model — annotation batched, threads only for structural ops

- **Annotation (TCR + MHC) is never Python-threaded and never per-structure.** It is one
  batched mmseqs2 call; mmseqs2 is the parallel layer (do NOT pass it a thread count). No
  `ProcessPoolExecutor`/`workers`.
- **Use threads ONLY for the embarrassingly-parallel, mmseqs-free stages:** structural
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
- **Φ and ΔΦ with the TP / TM / PM breakdown come from this one command.** Columns `F_tcr_pep`,
  `F_tcr_mhc`, `F_pep_mhc`, `F_total`; `run_pipeline(…, reference_aa="A")` / `tcren scoring --delta`
  adds `dF_*` and `dF_total` (each interface's `tcren.ddg.reference_delta`; `dF_tcr_mhc` ≡ 0, the
  peptide is not in that interface). Use ΔΦ, not Φ, when each candidate carries its **own**
  generated pose. **Column names are shared with `tcren recognize`** — one vocabulary, tables join
  on `pdb.id`. (Manuscript notation: Φ for the energy, φ for a potential matrix entry, F for the
  binder-direction channel `−Φ_TCR:pep` in `cohort.f_score`.)
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
  `src/_geom/geom.cpp` (interface geometry for `tcren binder`; also `interface_clashes` and
  `contact_stability`, the clash + fragility kernels). Adding a sixth = same pattern; a new *function*
  in an existing module (as clashes/stability were added to `_geom`) needs no `CMakeLists.txt` change.
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
- `contact_stability(structure, cutoff=5.0, delta=1.0) -> StabilityReport` (`tcren.stability`) — TCR:peptide
  contact fragility read straight off the contact map (native `_geom.contact_stability`, numpy reference
  behind it). Per contact `margin = cutoff − dmin`; report has `n_contacts`, `mean_margin`, `frac_robust`,
  `frac_marg_lt1`, `exp_lost` (expected contacts lost under a `delta`-Å shift). A coordinate-only interface
  positional-confidence readout. Self-check `python -m tcren.stability`.
- `tcren recognize` appends both as output columns: `n_clashes`, `clash_score`, `exp_lost`, `mean_margin`,
  `frac_robust` (extra columns, **not** part of the 35 `RECOGNITION_FEATURES` the models consume).
- `check_register(model, reference=None) -> RegisterReport` — always reports the clash burden; with a
  correctly-registered `reference` (crystal / trusted pose) it adds the **anchor-Cα RMSD** in the
  MHC-groove frame (`peptide_rmsd`) → `wrong_register` (True/False; `None` without a reference).
  **A heavy clash burden alone is NOT a register call** — AF swap models are routinely clashy; register
  needs the reference. (The ila1 CPL lesson: forced pose → raw TCRen ROC 0.35, recovered to ≈0.77 on the crystal register.)
- `fix_register(model, template, engine="ccd") -> ModelResult` — re-threads the model's peptide sequence
  onto `template`'s correctly-registered backbone and re-refines via `model_peptide` (the
  FlexPepDock-functional path; needs the `_refine`/`_fold` kernel). Equal peptide length required.
  Template-free re-docking (CCD to canonical anchor targets) is the future extension.

## Poly-alanine reference score — `tcren.ddg.reference_delta` (geometry-normalized TCRen)

- `reference_delta(cm, peptide, pot, interface="tcr_peptide", reference_aa="A") -> float` = ΔΦ = Φ(peptide)
  − Φ(all-`reference_aa` peptide). It is the **full-peptide alanine scan** (== sum of `alanine_scan().ddG`)
  and subtracts the pose's identity-independent geometry baseline Φ(polyAla).
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

## Intra-peptide term — `tcren.intra_peptide_energy` (the contacts a chain makes with ITSELF)

- Every Φ in the package sums over contacts between two **different** chains, so a peptide held in its
  bound conformation by its own side chains scores the same as one that is not. This term is that
  omission, and it is **off everywhere by default** — `scope="inter"`, `peptide_internal=False`,
  `intra_weight=0.0` all reproduce the previous output byte-for-byte.
- `all_atom_contacts(..., scope="inter"|"intra"|"all")`; `peptide_internal_contacts(structure,
  cutoff=4.0, min_seq_sep=3)` is the intra case plus the separation filter. **4 Å, not 5** (the vdW
  convention for intra-chain contacts, vs the interface 5 Å), and `|i−j| ≥ 3` because sequence
  neighbours touch by covalent geometry: over `tests/assets/pdb` the totals are 11 contacts at
  `|i−j| ≥ 3` vs 134 at `|i−j| ≥ 2`, the jump being `i`/`i+2` pairs of an extended chain.
- `intra_peptide_energy(cm, pot, peptide=None, contact_weight="residue") -> float`; needs
  `ContactMap.from_structure(..., peptide_internal=True)` (which stores the pairs *beside* `contacts`,
  never in it). `peptide=` threads a candidate onto the structure's peptide positions.
- **The potential is symmetrised**, `(F + Fᵀ)/2`. An intra-chain pair has no from/to orientation, and
  which residue lands on which side is the contact table's canonical `(chain.id, residue.index)`
  ordering, not chemistry. Matters for TCRen (directed), no-op for MJ. **Default MJ**: TCRen is derived
  from TCR↔peptide contacts and says nothing about a chain's contacts with itself.
- Wired in at three layers: `score_peptides(..., intra_weight=w, intra_potential=)` /
  `tcren score --intra-weight` (score = Φ + w·E_intra); `run_pipeline(..., intra_weight=w)` /
  `tcren scoring --intra-weight` (reports `F_pep_int` raw, folds w·it into `F_total`; potential
  overridable via `potentials={"peptide_internal": …}`); `tcren recognize --full` (`F_pep_int`,
  `n_pep_int`, catalogued `involves_tcr=False` — the peptide's own conformation is cohort identity).
- **Expect a sparse term.** A canonical extended class-I 9-mer makes 0–1 internal contacts, so it only
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
- **For a cohort, prefer `tcren recognize --scores --mechanics -t 0`** over running `recognize` and
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
- Self-check (no PDB): `conda run -n tcren-fold python -m tcren.mechanics`.

## Docking geometry — `tcren.orient.docking` + `tcren.orient.tcrdock_geometry`

- **Two interpretable angles** (existing): `docking_angles(structure) -> DockingAngles(crossing_angle,
  crossing_angle_signed, incident_angle, ...)`. `crossing_angle` = the groove-plane "scanning" angle,
  `incident_angle` = the tilt. Computed from the Vα→Vβ axis in the groove frame; no reference DB.
- **Full rigid-body pose** (new, `tcrdock_geometry.py`): `docking_geometry(structure) -> DockingGeometry(d,
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
- **`CANONICAL_AXES` is the naming of `tcren.orient.frame`, and the gizmo prints it**: x=`width`
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
- Chain roles after `tcren orient`: A=Vα, B=Vβ, C=peptide, D=MHCα, E=MHCβ/β2m — `CHAIN_COLOURS`.
- Notebooks: `notebooks/pymol_canonical_figures.ipynb` (static gallery) and
  `notebooks/pymol_interactive.py` (**marimo**, `pip install 'tcren[marimo]'`). marimo is reactive,
  so a name may be defined by exactly one cell — import shared symbols once in the setup cell and
  thread them through the signatures, or the export fails with `MultipleDefinitionError`.
- Docs gallery page with rendered examples: `docs/gallery.rst` + `docs/_static/gallery/*.png`.
- Everything except `render`/`probe_rotation` is pure Python and unit-tested without PyMOL
  (`tests/unit/test_pymol_viz.py`); the two that shell out are `slow` + skip when pymol is absent.

## Contact typing — `tcren.contact_types` (DSSP-style, dep-light)

- `contact_type_counts(cm, interface='tcr_peptide', tcr_regions='all') -> {n_<type>, pairs_<type>}` and
  `classify_contacts(interface_df) -> df + 'contact.type'`. Types by priority: `salt_bridge`, `hydrogen_bond`,
  `aromatic`, `hydrophobic`, `other`, from heavy-atom geometry (no H, no external DSSP). `pairs_hydrogen_bond`
  is the documented, reproducible replacement for the lost ad-hoc `n_hbond` (tracks it at r≈0.68).

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
- **Finding (2026-07-05):** real-vs-shuffled is learnable at **AUC 0.876** (RF; F_tcr_pep the top feature).
  But a shuffled-trained (crystal) model does **not** transfer to AF-modeled TCRvdb (0.55–0.62 vs AF 0.79) —
  crystal→AF distribution shift. Use it as a label-free recognition prior / supplementary benchmark, not as a
  drop-in TCRvdb scorer.
- **BN classifier** (`tcren.recognition.GaussianBNClassifier`): pure-numpy conditional-linear-Gaussian Bayes
  net — DAG over standardized features (BIC hill-climb on within-class-centred data) + class `y` and MHC class
  as discrete parent nodes; classifies by the Gaussian log-likelihood ratio. `fit/predict_proba`, gzip-JSON
  `save/load`, `to_dot` (graphviz). Trained model shipped at `src/tcren/data/shuffle_bn.json.gz`; the
  reproducible appendix (train+eval, gnuplot ROC/PR + balanced metrics, graphviz BN, marginals) is
  `shuffle_bn/` in the technical appendix (`make`) — moved to the manuscript repo 2026-07-28,
  `2026-tcren/archive/tcren-appendix/`; still in this repo's git history. Decoys are regenerable (`tcren shuffle --seed 0 --n 10`); manifest committed,
  full PDBs belong on HF (351 MB).
- **Distribution-aware logistic** (`tcren.recognition.BayesianLogisticRecognizer` + `encode_features`): a
  *discriminative* alternative to the BN. `encode_features` maps each feature by its natural family — circular
  `dock_torsion` → (cos, sin) von-Mises stats, `chain_balance` → logit, counts/continuous linear, drops the
  duplicate `n_hbond` — then a Bayesian logistic (PyMC NUTS, weakly-informative or horseshoe prior) is fit and
  its posterior mean frozen into `src/tcren/data/shuffle_logistic.json.gz` (dep-light numpy `predict_proba`).
  Real-vs-shuffled 5-fold CV **ROC-AUC 0.885** (matches RF, > BN 0.865 / raw-logistic 0.870). On TCRvdb a
  *supervised* refit with the encoding gives **0.860** pooled (> AF 0.794, raw-feature logistic 0.855); the
  frozen real-vs-shuffled transfer does NOT carry (0.53, crystal→AF shift, same as the BN). Appendix
  `logistic_stan/` (`make PY=<pymc-venv>`; ROC/PR + posterior-forest gnuplot, encoding table).
- **`tcren recognize` / `recognition_features` (2026-07-06):** `recognition.recognition_features(struct)`
  ports the manuscript's 35-descriptor extractor into tcren (docking geometry + TCRen/MJ F & poly-Ala dF +
  contact tallies + biopython ΔSASA `burial` + `mhc_class_bin`) — verified **byte-exact** vs
  `canonical2026_features.csv` (burial max diff 4e-11). Uses `import_structure` (C-gene trimmed) to match
  training; **no `_geom` C-ext needed** (only arda for annotation). `frozen_recognizers()` loads both
  shipped models; `real_probability(rows)` → `{"logistic","bn"}` P(real). CLI `tcren recognize -s pdbs/ -o
  out.tsv` writes one TSV row/PDB = 35 descriptors + `p_real` + `p_real_bn` (`--features-only` skips models).
  The user-facing "one TSV for a/b/d" answer; koff joins it under `--mechanics` (2026-07-28), and only
  ddF (ala), which is per-residue rather than per-structure, stays its own command `tcren ddg`.
- **`--full` feature table (2026-07-13, audited 2026-07-28):** `recognition_features(struct, full=True)` /
  `tcren recognize --full` append the **18 CDR3-frame** descriptors
  (`cdr3{a,b}_{reach,ou,ow,on,au,aw,an,topep,ext}`, FramePose groove-frame projection — the `cdr3b_*`
  strain signal) → 52 features total. Tuples: `RECOGNITION_FEATURES` (34), `CDR3_FRAME_FEATURES` (18),
  `FULL_FEATURES` (52).
- **Descriptor audit (2026-07-28):** every energy column is `F_*` (`e_cdr12`/`e_cdr3a`/`e_cdr3b` →
  `F_cdr12`/`F_cdr3a`/`F_cdr3b`); the duplicate `e_tcr_mhc` and `ct_tp_hydrogen_bond` columns are gone
  (they equalled `F_tcr_mhc` and `n_hbond`); the **12 matrix-swap** columns
  (`{tcren,mj,d}_{tp,cdr12,cdr3a,cdr3b}`) were **removed** — `tcren_*` duplicated `F_*`, and MJ is not the
  potential used on TCR:peptide. New: `crossing_signed` (signed scanning angle, carries docking polarity)
  and `DESCRIPTORS` / `descriptors(family, tcr_only=)` — the catalogue giving each column's family
  (`geometry`/`physics`/`kinetics`/`score`) and whether the receptor enters it. Only `F_pep_mhc`,
  `dF_pep_mhc` and `mhc_class_bin` do not; they carry cohort identity, so receptor questions must use
  `tcr_only=True`. Frozen recognizers verified **bit-identical** through `_FROZEN_ALIASES`.
- **`--scores` good-results scores (2026-07-13):** `tcren recognize --scores` (implies `--full`) also emits
  `p_bind` (`binder.binder_score`, TCRvdb denoised AUC 0.928) and `p_forced`
  (`recognition.forced_pose_score` / `FORCED_POSE_MODEL`, a frozen 6-feature strain logistic: crystal-natural
  vs AF-forced, 5-fold AUC 0.762 — the "too-good-to-be-true" hallucination flag). Cohort-relative z-combos
  (Φ_bind, the crystal<af_real<af_decoy strain gradient) stay analysis-side, computed from these features.
  Full column reference: `docs/features.rst`.
- **`kit_score` — AF×tcren synergy (2026-07-13):** `recognition.kit_score(p_bind, iptm)` = cohort-relative
  `z(p_bind)+z(iptm)`, the fixed no-fit combination of the tcren binder score with the AF ipTM. On TCRvdb raw
  labels it beats **either alone** at precision (macro-PR 0.847 vs 0.782/0.804; P@10% 0.969; Δ vs ipTM +0.041
  CI [+0.006,+0.074]) and corrects AF errors (strain flags AF false-positives 0.633; p_bind rescues AF
  false-negatives 0.732>ipTM 0.697). The "kit for AI-generated structures" decision procedure is `docs/kit.rst`.
  (No synergy on VDJdb real-vs-mock — there tcren's role is the forced-pose gradient, not discrimination.)
- **`-t/--threads` on `tcren scoring` and `tcren recognize` (2026-07-26):** both accept a file, a
  directory, a `.tar.gz`, a quoted glob or a `.txt` manifest; `-t N` runs N concurrent workers (`-t 0`
  = all cores). Cohort-relative scores (`q_bind`, `s_strain`) are still computed over the **whole** set,
  never per batch. `scoring` gains ~7.6x on 8 threads; `recognize` less (its cost is Python
  featurisation, not mmseqs), so batch its annotation rather than expecting linear scaling.

- **`cohort.q_iptm` — fit-free `z(ipTM)+z(Q)` (2026-07-24, v2.3.1):** `cohort.q_iptm(table, iptm, features=Q_FEATURES)`
  ships the geometry synergy as one call (the fit-free analog of `kit_score`, which pairs ipTM with the *fitted*
  `p_bind`). `cohort.Q_FEATURES_GEOM` is the 4 geometry-only terms (`Q_FEATURES` minus `pp_combo`) → `Q_geom`,
  robust to the forced-pose energy inversion (C27). `z(ipTM)+z(Q_geom)` beats raw-AF ipTM on both ROC and PR on
  well-modelled ("template-covered") epitopes and ties it fit-free on TCRvdb (benchmark C42). Single-line CLI:
  `tcren recognize -s pdbs/ --iptm meta.tsv -o out.tsv` joins ipTM (key col matched to `complex.id`) and appends
  `Q_geom` + `z(ipTM)+z(Q_geom)`.
- **`cohort.q_coupled` / `cohort.coupling` — the parameter-free binder score (2026-07-26):**
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

- **`cohort.f_score` / `cohort.q_f` — the contact-energy channel (2026-07-24, v2.3.2):** `f_score(table)` =
  `z(-(F_tcr_pep+F_tcr_mhc))` (binder-oriented, `cohort.F_TERMS`); `q_f(table, sign=+1)` = `z(Q_geom)+sign·z(F)`,
  the pure-tcren combiner with **no DL term**. F reads contact chemistry but is **pose-conditional**: it inverts
  on forced poses (GLC↔ila1, C27). On template-covered poses `z(Q)+z(F)` beats raw-AF ipTM on both ROC/PR
  (0.759/0.725 vs 0.692/0.693, C42); on forced poses read `z(Q)-z(F)` (`sign=-1`; GLCTLVAML 0.71 vs 0.52 for
  `+F`). `tcren recognize --cohort` emits `Q_geom`, `F_score`, `z(Q)+z(F)`, `z(Q)-z(F)` (and `z(ipTM)+z(Q)+z(F)`
  with `--iptm`) — the full fit-free panel in one line. `z(ipTM)+z(Q_geom)` stays the channel robust to the
  inversion; grade forced-ness with `s_strain` before trusting `+F`. Reproduced by
  `models/qf_panel.py` in the benchmark repo.
- **AlphaFold-synergy auto-inversion (v2.3.2):** `cohort.q_f_iptm(table, iptm, threshold=0.5)` uses AF's own
  ipTM to pick F's sign per structure — `+z(F)` on confident poses, `-z(F)` on forced (low-ipTM) ones —
  and `cohort.f_invert_by_iptm(iptm, thr)` is the boolean flag. With `--iptm`, `tcren recognize --cohort`
  emits `F_invert` + `z(Q)+z(F|iptm)` (threshold `--invert-f-thresh`) and **prints how many poses are
  forced** so the user is told F inversion is in play; without `--iptm` it warns that F is trusted
  unconditionally. The versatile AF-post-analysis path: geometry (`Q`), chemistry (`F`), and AF confidence
  (ipTM) combined, with F used only where AF's pose is trustworthy.

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
- Root `data/` holds the library dataset (gitignored structures): `Native2026` (orientation
  references), `Canonical2026` (the default `superimpose` database), `PDB_date.tsv`,
  `TCRen_potential.csv`. `Canonical2026`'s `orient_metadata.json` moved into
  `src/tcren/data/` — it must ship in the wheel, since `fetch-data` downloads only structures and
  an installed user has no repo `data/`. `setup.sh` runs `tcren fetch-data` at install to
  populate `Native2026` + `Canonical2026` from HF (or lazily on first `superimpose`/`orient`).
  Orientation references load 1ao7/1fyt from `data/Native2026` via `tcren.paths`. The numerical
  regression oracle (legacy mir/R outputs: `contact_maps_PDB.csv`, `tcren_am/tcren.txt`, the
  `example/` set) lives under `tests/assets/oracle/`; the legacy R/Java pipeline was deleted.

## Symmetric TCRen — `derive_tcren(..., symmetric=True)`

TCRen is **directed**: `from` is a TCR residue, `to` a peptide residue, so `F[a,b] != F[b,a]` and
the shipped matrix is 19x20 (classic drops the `from == "C"` row). MJ — used for the presentation
interfaces — is symmetric 20x20. `symmetric=True` folds the **raw counts** onto their transpose
(`potential.symmetrize_counts`, `N + N.T`) *before* the log-odds, so the marginals and hence the
expected term are rebuilt from the folded counts.

```python
from tcren.potential import derive_tcren, symmetrize_counts
sym = derive_tcren(contacts, variant="classic", symmetric=True)   # exactly symmetric, 20x20
```

- **Not a post-hoc average.** Averaging the finished energies leaves the asymmetric background in
  place; on Native2026 the two differ by 0.29 mean / 0.82 max. Symmetrise counts, not energies.
- **Cys is grafted, not dropped.** Free Cys is essentially absent from CDR loops — 4 of 8062
  Native2026 contacts (0.05%) are TCR-side Cys vs 32 (0.40%) peptide-side. The fold grafts the
  peptide-side observations onto the Cys row, so the symmetric matrix keeps a full 20x20 alphabet
  instead of dropping a column for lack of data. `drop_cys` is forced `False` when `symmetric`.
- **Why it is interesting:** it puts TCRen in MJ's functional form, making them directly
  comparable for the first time. They turn out to be nearly **uncorrelated** (Pearson +0.07,
  Spearman +0.10 over the 210 unordered pairs) — TCR:peptide recognition chemistry is not generic
  protein-folding contact chemistry, which is the argument for TCRen existing separately at all.
  Symmetric vs directed TCRen: r = +0.67, mean|diff| 0.29.
- Default is `symmetric=False`; the shipped `TCRen_potential.csv` is unchanged.

## Geometry: contacts, region pairs, docking angle

```python
from tcren.project2d import region_pair_contacts, region_pair_summary  # needs chain-typed + MHC-annotated structure
from tcren.orient import docking_angles
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
