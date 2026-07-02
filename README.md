<p align="center">
  <picture>
    <source media="(prefers-color-scheme: light)" srcset="assets/tcren_light.svg">
    <img alt="tcren" src="assets/tcren_dark.svg" width="340">
  </picture>
</p>

<h1 align="center">tcren — structure-based prediction of TCR–epitope recognition</h1>

<p align="center">
  <a href="https://pypi.org/project/tcren/"><img alt="PyPI" src="https://img.shields.io/pypi/v/tcren"></a>
  <a href="https://github.com/antigenomics/tcren/actions/workflows/tests.yml"><img alt="tests" src="https://github.com/antigenomics/tcren/actions/workflows/tests.yml/badge.svg"></a>
  <a href="https://docs.isalgo.dev/tcren/"><img alt="docs" src="https://github.com/antigenomics/tcren/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-GPLv3-green"></a>
</p>

**TCRen** predicts which epitopes a T-cell receptor recognises from a single TCR–peptide–MHC
structure (experimental or modelled). It extracts the TCR–peptide contact map and scores every
candidate peptide with a **residue-level statistical potential** derived from contact preferences
in TCR:pMHC crystal structures — answering not "what fancy complex can a model draw?" but "is this
binding physically plausible?".

This is a documented, tested, CLI-driven Python library. TCR chains are annotated with the sibling
[`arda`](https://github.com/antigenomics/arda); MHC chains are mapped and the groove partitioned
against a curated reference; structures are oriented into one canonical frame; and the original
contact maps, potential, and scores are reproduced numerically (validated against committed oracles
to floating-point precision).

While the original tcren focused on TCR:peptide contacts, the new version brings in features to 
score TCR:MHC and peptide:MHC interactions, required to get full picture of TCR:pMHC binding 
mechanics and estimate ddG values.

## What it does

From one TCR–peptide–MHC structure (crystal or model), each task is one command or one call:

| task | command | library |
|---|---|---|
| Score candidate epitopes for a TCR | `tcren score` | `score_peptides` |
| Percentile-rank a peptide vs background | `tcren rank` | `percentile_rank` |
| ΔΔG of mutations (alanine scan / neoantigen) | `tcren ddg` | `alanine_scan`, `neoantigen_ddg` |
| Binder vs non-binder for a TCR model | `tcren binder` | `binder_score` |
| Three-interface energy breakdown + total | `tcren pipeline` | `run_pipeline` |
| Annotate chains + region markup | `tcren annotate` | `classify_chains`, `annotate_mhc` |
| Interface contact table (5/8/12 Å) | `tcren contacts` | `ContactMap`, `multi_contacts` |
| Orient into the canonical MHC frame | `tcren superimpose` / `orient` | `superimpose`, `canonicalize_structure` |
| Substitute a peptide + refine its pose | `tcren refine` | `substitute_peptide`, `refine_peptide` |
| DOPE interface energy (ΔΔG `e_native`) | `tcren energy` | `interface_energy` |
| Interface mechanics — koff proxies (stiffness / rupture) | `tcren mechanics` | `stiffness_tensor`, `rupture`, `coupling_residues` |
| Re-derive the statistical potential | `tcren derive-potential` | `derive_tcren` |
| Steric-clash / wrong-register QC | — | `interface_clashes`, `check_register` |
| 2D complementarity map + 3D pocket/CDR view | — | `render_complementarity_map`, `view_pocket_cdr` |

## Install

```fish
pip install tcren          # from PyPI — binary wheels ship the C++ extension; pulls in arda-mapper
```

For development (editable install, conda env with the build toolchain, and the reference data
fetched into `data/`):

```fish
bash setup.sh              # creates the `tcren` conda env, installs arda + tcren, fetches data/
conda activate tcren
```

tcren ships five small **pybind11/C++ extensions**, built on install by `scikit-build-core`:
`tcren._align` (MHC-pseudosequence fitting alignment; a Biopython fallback runs if unbuilt),
`tcren._refine` (DOPE atom-level Monte-Carlo peptide refinement), `tcren._relax` (DOPE interface
energy for `tcren energy` / ΔΔG), `tcren._fold` (CCD loop closure) and `tcren._geom` (interface
geometry for `tcren binder`). TCR annotation is provided by [`arda`](https://github.com/antigenomics/arda), a runtime
dependency published to PyPI as [`arda-mapper`](https://pypi.org/project/arda-mapper/) (it imports
as `arda`); `pip`/`setup.sh` pull it automatically, and from `arda-mapper >= 2.0.3` it auto-fetches
its own reference on first use (no `ARDA_HOME` to set). `setup.sh` also runs `tcren fetch-data` to
populate `data/` with the reference structure sets (`Native2026`, `Canonical2026`) used by
`orient`/`superimpose` (set `TCREN_NO_FETCH=1` to skip).

## Command line

```fish
# Full pipeline: annotate -> superimpose -> resmarkup / canonical Cα / contacts -> per-interface
# energies (TCRen for TCR↔peptide, MJ for TCR↔MHC and peptide↔MHC) + total
tcren pipeline -s complex.pdb -o scores.csv

# Configurable per-interface potential: swap a bundled name (tcren|mj|keskin), a CSV, or
# None for any interface; default reproduces the built-in per-interface families exactly.
tcren pipeline -s complex.pdb -o scores.csv --tcr-mhc-potential keskin

# Opt-in TCR framework regions: --regions {all,cdr,cdr+fr} chooses which TCR regions
# contribute on the TCR side (cdr = CDR1-3 only; cdr+fr adds FR1-3; all = unfiltered, default).
tcren score -s complex.pdb -c candidates.txt -o ranked.csv --regions cdr+fr

# Percentile-rank the native (or candidate) peptide's TCRen energy against a random pMHC
# background — small rank_pct = the peptide scores among the best binders.
tcren rank -s complex.pdb -o rank.csv

# Fast ΔΔG of peptide point mutations (virtual-matrix path: no atoms move, no re-docking).
tcren ddg -s complex.pdb -o ddg.csv

# Binder vs non-binder P(binder) from AF-orthogonal interface geometry + the CDR1/2-vs-CDR3a
# TCRen term — ranks candidate TCRs against a fixed pMHC, beating AlphaFold/TCRmodel2 confidence
# (denoised AUC 0.928 vs 0.872) with no external tool. See tcren.binder.binder_score.
tcren binder -s complex.pdb -o binder.csv

# End-to-end candidate-epitope scoring from a structure
tcren score -s complex.pdb -c candidates.txt -o ranked.csv

# Substitute a peptide and refine its pose (knowledge-based MC scored by the DOPE atom-level
# statistical potential — independent of the TCRen/MJ scoring potentials, restrained to the input).
# Not physics relaxation — use Rosetta FlexPepDock for that.
tcren refine -s complex.pdb -o refined/ --substitute KQWLVWLFL

# Structures: any of .pdb / .cif / .pdb.gz / .cif.gz, a directory, or a .tar.gz batch
tcren contacts -s batch.tar.gz -o contacts.csv --interface tcr_peptide

# Per-residue markup: TCR (CDR/FR) + MHC groove (helix/floor) + peptide in one table.
# --regions all|tcr|mhc|peptide filters; --pseudo also marks NetMHCpan groove residues (MPS).
tcren annotate -s complex.cif.gz -o markup.csv --regions mhc --pseudo

# Superimpose structure(s) onto the canonical frame, by MHC, against the canonical database
# (data/Canonical2026, fetched at install). Detects MHC class + species and averages the
# superposition over every database structure of that class/species. Chains -> A=Vα B=Vβ
# C=peptide D=MHCα E=MHCβ/β2m. -s takes a file / directory / .tar.gz / glob; -o is a directory,
# or a single structure file (one input) whose extension must match --mmCIF/--compress; -t threads.
tcren superimpose -s complex.pdb -o oriented.pdb           # single file
tcren superimpose -s 'data/*.pdb' -o oriented/ -t 8        # glob -> directory, threaded

# Build a canonical database from native complexes (how Canonical2026 is produced). Annotation
# is one batched mmseqs call; -t threads only the structural alignment + write.
tcren orient -s data/Native2026 -o data/Canonical2026 -t 8

# Structure outputs are plain .pdb by default; add --mmCIF for .cif and --compress for .gz.
tcren superimpose -s complex.pdb -o oriented/ --mmCIF --compress   # -> oriented/<id>.cif.gz

# Fetch recent TCR-pMHC structures from RCSB -> data/pdb_recent (mmCIF .cif.gz, 5-chain validated)
tcren fetch-recent --discover --after 2024-01-01

# Build the MHC reference once (IMGT/HLA + mouse H-2; cached, not committed)
tcren build-mhc-ref

tcren info
tcren --install-completion        # shell tab-completion (bash/zsh/fish)
```

`tcren orient` and `tcren superimpose` need the reference sets in `data/` (`Native2026`,
`Canonical2026`); `setup.sh` fetches them at install via `tcren fetch-data` (re-run it any time).

## Library

```python
from tcren import run_pipeline, parse_structure, import_structure, ContactMap, score_peptides
from tcren.annotation import classify_chains
from tcren.potential import tcren

# One call: annotate -> superimpose -> contacts -> per-interface energies + total
res = run_pipeline("complex.pdb")              # res.scores, res.markup, res.contacts, res.oriented

# Oracle facade: one structure -> a bundle of ready-to-tabulate frames for the paper
# notebooks (scores, percentile rank, ΔΔG alanine scan, markup, contacts). Configurable
# per-interface potentials and TCR-region selection are forwarded to every milestone.
from tcren import summarize_structure
bundle = summarize_structure("complex.pdb", alanine=True)   # bundle["scores"], ["rank"], ["ddg"], …

# …or the individual steps:
s = parse_structure("complex.pdb.gz")          # also .cif/.cif.gz; import_structure trims the C-gene
classify_chains(s, organism="human")           # TRA/TRB via arda, peptide, MHC
cm = ContactMap.from_structure(s)              # 5 Å contacts + interface partitioning
ranked = score_peptides(cm, ["KQWLVWLFL", "RLLHPHHPL"], tcren())
```

### Batch inputs, gzip, archives

```python
from tcren.structure import iter_structures
for pdb_id, structure in iter_structures("batch.tar.gz"):   # file | directory | .tar.gz
    classify_chains(structure, organism="human")
    ...
```

### Canonical orientation, contacts, docking geometry

```python
from tcren.mhc import annotate_mhc
from tcren.orient import canonicalize_structure, superimpose, docking_angles
from tcren.contacts import multi_contacts, ContactDefinition

annotate_mhc(s)
oriented, info = canonicalize_structure(s)     # frame: z=MHC→TCR, y=peptide, x=thin; chains A–E
oriented, info = superimpose(s)                # orient onto data/Canonical2026 by MHC (class+species ensemble)
layers = multi_contacts(s, ContactDefinition(d1=5, d2=8, d3=12))   # heavy-atom / Cβ / Cα
d = docking_angles(s)                          # crossing (~20–70° αβ) + incident angle
```

### 2D complementarity maps & region-pair contacts

```python
from tcren.project2d import (project_structure, residue_markup_table, contacts_table,
                             region_pair_summary)
from tcren.viz import render_complementarity_map, view_pocket_cdr

proj = project_structure(s)                                   # canonical groove plane
svg  = render_complementarity_map(residue_markup_table(s, proj),
                                  contacts=contacts_table(s, threshold=5.0))
region_pair_summary(s, kind="closest")        # contacts per region pair + bond types (cb/ca too)
view_pocket_cdr(s).show()                      # interactive 3D pocket + CDR overlay (py3Dmol)
```

## Modules

| module | what it does |
|---|---|
| `tcren.structure` | parse/write `.pdb`/`.cif`(`.gz`)/`.tar.gz`; the `Atom`/`Residue`/`Chain`/`Structure` model; `iter_structures` |
| `tcren.annotation` | chain typing — TCR loci/CDRs via `arda`, peptide, MHC; αβ/γδ C-gene call |
| `tcren.mhc` | map MHC chains to allele/class/role; partition the groove (helices/floor); NetMHCpan pseudosequence |
| `tcren.contacts` / `contactmap` | closest-atom 5 Å contacts, Cα distances, multi-layer (5/8/12 Å) contact tables, interface partitioning |
| `tcren.potential` | `Potential` (TCRen/MJ/Keskin); `derive_tcren` (classic/AM/LOO) with non-redundancy filtering |
| `tcren.scoring` / `scoring_rank` | substitution scoring of candidate peptides; percentile rank vs a background |
| `tcren.ddg` | fast virtual-matrix ΔΔG — alanine scan, neoantigen mutants |
| `tcren.binder` | binder/non-binder classifier from AF-orthogonal interface geometry |
| `tcren.orient` | canonical frame, `superimpose` onto the canonical DB, docking angles, reverse-dock detection |
| `tcren.refine` | peptide substitution + refinement (DOPE MC; CCD/OpenMM/ProMod3/FlexPepDock engines); register QC |
| `tcren.clashes` / `mechanics` | steric-clash report; interface spring-network stiffness + rupture model |
| `tcren.project2d` / `viz` | project the interface onto the groove plane; SVG complementarity maps + 3D pocket/CDR views |
| `tcren.pipeline` / `oracle` | one-call end-to-end runs (`run_pipeline`, `summarize_structure`) |
| `tcren.paper` | Nat Comput Sci 2022 reproduction (HF bootstrap, batch annotation, legacy comparison) |

## Data

Structures live in the Hugging Face dataset
[`isalgo/tcren_structures`](https://huggingface.co/datasets/isalgo/tcren_structures), all gzipped:

| folder | contents |
|---|---|
| `Native2022` | the 2022 paper set (oracle) |
| `Native2026` | the comprehensive 2026 TCR:pMHC set the current potential is derived from |
| `Canonical2026` | `Native2026` re-oriented into the canonical frame (`tcren orient`) |

`tcren` reads `.pdb`/`.cif`/`.pdb.gz`/`.cif.gz` and `.tar.gz` batches; an installed library lazily
fetches the canonical reference structures from the Hub when orienting a new complex. The root
`data/` holds `Native2026` (+ `Canonical2026`, gitignored, fetched on demand), `PDB_date.tsv`,
`orient_metadata.json`, and **`TCRen_potential.csv`** — the current potential derived from the
Native2026 set (use it with `tcren score -p data/TCRen_potential.csv`).

## Notebooks

Runnable examples under [`notebooks/`](notebooks/) (rendered in the
[docs](https://docs.isalgo.dev/tcren/)):

- `complementarity_map_2d` — 2D interface maps, multiple structural + map views of 1ao7
- `contact_thresholds_and_bondtypes` — region-pair contact counts (closest/Cβ/Cα) + bond types
- `canonical_frame_figures` — canonical-frame QC across the Native2026 set
- `pymol_canonical_figures` — ray-traced PyMOL panels (overlay, groove, interface) by class/species
- `mhc_pseudosequence_mps` — NetMHCpan MHC pseudosequence (MPS) residues vs. peptide contacts
- `example_gil_a02_rs_motif` — GILGFVFTL/HLA-A*02 and the public CDR3β Arg–Ser motif
- `natcompsci2022/` — full reproduction of the Nat Comput Sci 2022 analyses

## Performance

Per-stage wall time (best of *n*) on a TCR-pMHC complex (1ao7), Apple M-series, single thread
(`RUN_BENCHMARK=1 pytest -k benchmark -s` to reproduce the core stages):

| stage | time | notes |
|---|---|---|
| parse a gzipped structure | ~17 ms | `.pdb.gz` / `.cif.gz` |
| contact map (5 Å, cKDTree) | ~9 ms | per structure |
| score 1000 candidate peptides | ~11 ms | ~10 µs/peptide (vectorised) |
| ΔΔG alanine scan (9-mer) | ~11 ms | virtual-matrix; no atoms move |
| binder P(bind) (features + model) | ~49 ms | native geometry, no external tool |
| peptide refine (2000-step DOPE MC) | ~320 ms | knowledge-based rigid-body refinement |
| annotate (MHC map, 1 structure) | ~670 ms | one mmseqs2 search |
| **annotate (TCR + MHC), batched** | **~0.2 s/structure** | one mmseqs2 call for the whole set; vs ~1.5 s/structure unbatched |
| superimpose onto the canonical DB (per query) | ~2.8 s | aligns to every same-class DB structure |

| peak RSS | value | notes |
|---|---|---|
| single-structure pipeline (no orient) | ~200 MB | parse → annotate → contacts → score → refine |
| + `superimpose` (loads canonical DB) | ~780 MB | holds Canonical2026 in RAM; skip with `--no-superimpose` |

Annotation is the only network/compute-heavy step and is always **batched** (one mmseqs2 search over
all chains; mmseqs2 parallelises internally — never per-structure, never Python-threaded). Threads are
used only for the embarrassingly-parallel, mmseqs-free stages (structural alignment, write, rendering):
`tcren orient -t N`. Screening a peptide/TCR panel is embarrassingly parallel — references are
annotated and oriented **once**, so the hot loop is just refine + contacts + score per complex.

## Tests

```fish
pytest -m "not slow"          # unit + fast regression (the CI gate)
pytest                        # add the arda/mmseqs-backed regression tests
RUN_BENCHMARK=1 pytest -k benchmark -s
```

## Methods appendix

The coordinate-level extensions — backbone-preserving peptide substitution and the potential-guided
Monte-Carlo refinement kernel (energy function, the restraint-necessity argument, sampler, and
citations) — are written up in the technical appendix [`appendix/tcren.tex`](appendix/tcren.tex)
(built with `make -C appendix` → `appendix/tcren.pdf`).

## Citing

**TCRen** is free for academic and non-commercial use. If you use it, please cite our latest 
[Nature Computational Science 2024 paper](https://www.nature.com/articles/s43588-024-00653-0):

```
Karnaukhov VK, Shcherbinin DS, Chugunov AO, Chudakov DM, Efremov RG, Zvyagin IV, Shugay M. Structure-based prediction of T cell receptor recognition of unseen epitopes using TCRen. Nat Comput Sci. 2024 Jul;4(7):510-521. doi: 10.1038/s43588-024-00653-0. Epub 2024 Jul 10. PMID: 38987378.
```
