---
name: tcren
description: tcren — TCR-pMHC contact potential (TCRen) pipeline; conventions and public API
---

# tcren Skills Guide

`tcren` reproduces and extends the TCRen contact-energy potential (Nat Comput Sci 2022)
on a pure-Python pipeline (structure parsing → contacts → TCR/MHC annotation → potential
derivation → epitope-ranking benchmarks). Annotation uses the `arda` package
(mmseqs2-backed), a runtime dependency published to PyPI as `arda-mapper` (imports as `arda`,
`>=2.0.3`) — no separate checkout and **no `ARDA_HOME`** (arda auto-fetches its reference into
`~/.cache/arda` on first use). Conda env `tcren` (`bash setup.sh`).

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
  `orient_metadata.json`), superposes the query groove Cα onto each, and **averages** the rigid
  transforms (chordal/SVD mean rotation + mean translation) into one consensus placement. The
  matching DB subset is batch-annotated once and cached per process.
- **`tcren orient` / `tcren.orient.run_folder(...)`** — BUILD a canonical DB from native
  complexes using the per-class derived frame (how `Canonical2026` is produced). Not for orienting
  a single new structure — use `superimpose` for that.
- **HF upload is NOT a user command.** `--push-to-hub` was removed; maintainers run
  `scripts/push_canonical_to_hub.py` instead.

## End-to-end pipeline — `tcren.run_pipeline` / `tcren pipeline`

- `run_pipeline(structure, superimpose=True, db_dir=…)` → `PipelineResult`: import → annotate
  (alleles + chains + MHC groove) → superimpose onto the canonical DB (canonical Cα) → resmarkup
  + 5 Å contacts → per-interface energies. Scores: **TCRen** for TCR↔peptide, **MJ** for TCR↔MHC
  and peptide↔MHC, plus `total` (sum of the residue-pair potential over each interface's
  contacts). CLI `tcren pipeline -s … -o scores.csv` writes one row per structure.

## Compiled extensions — `tcren._align`, `tcren._refine` (pybind11 / scikit-build-core)

- TWO C++ exts, both in `CMakeLists.txt` (`pybind11_add_module` each, `install(TARGETS _align
  _refine …)`): `src/_align/align.cpp` (MHC pseudoseq fitting alignment) and `src/_refine/refine.cpp`
  (potential-guided peptide refinement). Adding a third = same pattern.
- The MHC-pseudosequence fitting-alignment hot path is a C++ ext (`src/_align/align.cpp`,
  `CMakeLists.txt`). Build backend is `scikit-build-core` (not hatchling); `pip install -e .`
  builds it (editable.rebuild on import). Funcs: `fitting_score`, `best_hit` (GIL released over
  candidates), `align` (traceback). Scoring matches Bio.Align's fitting config EXACTLY (BLOSUM62,
  placed-gap open -11/extend -1, free target + end gaps), so `tcren.mhc.pseudo` falls back to
  Biopython transparently when the ext is absent. ~40 ms vs Bio 59 ms vs pure-Python 15 s for 4k
  candidates (a modest 1.5x — Bio's aligner is already C). `editable.rebuild = false` in
  pyproject: the ext is built once at `pip install -e .` (do NOT rebuild on import — that needs
  cmake on PATH at import time and breaks pytest/CI). CI installs `--no-deps` + explicit runtime
  deps (so arda-backed tests skip) and `pip install cmake ninja` to build the ext.

## Peptide substitution + refinement — `tcren.refine` (`tcren refine`)

- `substitute_peptide(structure, new_peptide)` — backbone-preserving identity swap on the peptide
  chain (keep N/Cα/C/O+Cβ, drop side-chain beyond Cβ); pure data-model, no atoms moved. `score_peptides`
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
- **Caveat: these track the dissociation off-rate koff / kinetic stability (Bell–Evans rupture
  resistance ~ r0.5 on ATLAS), NOT the equilibrium ΔG/Kd** — rupture reflects the dissociation
  barrier, not the well depth (physically apt for the TCR mechanosensor / catch bonds). Use them
  **between structures** (one value per complex) to rank/compare; do not pool many per-residue or
  per-spring rows from one structure as independent samples — that is pseudo-replication.
- Self-check (no PDB): `conda run -n tcren-fold python -m tcren.mechanics`.

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
  `orient_metadata.json`, `TCRen_potential.csv`. `setup.sh` runs `tcren fetch-data` at install to
  populate `Native2026` + `Canonical2026` from HF (or lazily on first `superimpose`/`orient`).
  Orientation references load 1ao7/1fyt from `data/Native2026` via `tcren.paths`. The numerical
  regression oracle (legacy mir/R outputs: `contact_maps_PDB.csv`, `tcren_am/tcren.txt`, the
  `example/` set) lives under `tests/assets/oracle/`; the legacy R/Java pipeline was deleted.

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
