<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/tcren_dark.svg">
    <img alt="tcren" src="assets/tcren_light.svg" width="340">
  </picture>
</p>

<h1 align="center">tcren — structure-based prediction of TCR–epitope recognition</h1>

<p align="center">
  <a href="https://github.com/antigenomics/tcren/actions/workflows/tests.yml"><img alt="tests" src="https://github.com/antigenomics/tcren/actions/workflows/tests.yml/badge.svg"></a>
  <a href="https://docs.isalgo.dev/tcren/"><img alt="docs" src="https://github.com/antigenomics/tcren/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <img alt="license" src="https://img.shields.io/badge/license-academic%20%2F%20non--commercial-green">
</p>

**TCRen** predicts which epitopes a T-cell receptor recognises from a single TCR–peptide–MHC
structure (experimental or modelled). It extracts the TCR–peptide contact map and scores every
candidate peptide with a **residue-level statistical potential** derived from contact preferences
in TCR:pMHC crystal structures — answering not "what fancy complex can a model draw?" but "is this
binding physically plausible?". Method: Karnaukhov *et al.*, *Nat Comput Sci* 2024
([doi:10.1038/s43588-024-00653-0](https://www.nature.com/articles/s43588-024-00653-0)).

This is a documented, tested, CLI-driven Python library. TCR chains are annotated with the sibling
[`arda`](https://github.com/antigenomics/arda); MHC chains are mapped and the groove partitioned
against a curated reference; structures are oriented into one canonical frame; and the original
contact maps, potential, and scores are reproduced numerically (validated against committed oracles
to floating-point precision).

## Install

```fish
bash setup.sh              # creates the `tcren` conda env, installs arda + tcren
conda activate tcren
```

`setup.sh` expects the sibling `arda` checkout next to this repo (or set `ARDA_DIR`).

## Command line

```fish
# End-to-end candidate-epitope scoring from a structure
tcren score -s complex.pdb -c candidates.txt -o ranked.csv

# Structures: any of .pdb / .cif / .pdb.gz / .cif.gz, a directory, or a .tar.gz batch
tcren contacts -s batch.tar.gz -o contacts.csv --interface tcr_peptide
tcren annotate -s complex.cif.gz -o markup.csv
tcren mhc      -s complex.pdb    -o mhc_calls.csv

# Canonical orientation (writes oriented .pdb.gz; chains A=Vα B=Vβ C=peptide D=MHCα E=MHCβ/β2m).
# Annotation is one batched mmseqs call; -t threads only the structural alignment + write.
tcren orient -s data/Native2026 -o data/Canonical2026 -t 8 \
    --push-to-hub isalgo/tcren_structures --hub-folder Canonical2026

# Fetch recent TCR-pMHC structures from RCSB -> data/pdb_recent (mmCIF .cif.gz, 5-chain validated)
tcren fetch-recent --discover --after 2024-01-01

# Build the MHC reference once (IMGT/HLA + mouse H-2; cached, not committed)
tcren build-mhc-ref

tcren info
```

## Library

```python
from tcren import parse_structure, import_structure, ContactMap, score_peptides
from tcren.annotation import classify_chains
from tcren.potential import tcren

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
from tcren.orient import canonicalize_structure, align_to_canonical, docking_angles
from tcren.contacts import multi_contacts, ContactDefinition

annotate_mhc(s)
oriented, info = canonicalize_structure(s)     # frame: z=MHC→TCR, y=peptide, x=thin; chains A–E
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
`data/` holds `Native2026` (+ `Canonical2026`, gitignored, fetched on demand), `PDB_date.tsv`, and
`orient_metadata.json`.

## Notebooks

Runnable examples under [`notebooks/`](notebooks/) (rendered in the
[docs](https://docs.isalgo.dev/tcren/)):

- `complementarity_map_2d` — 2D interface maps, multiple structural + map views of 1ao7
- `contact_thresholds_and_bondtypes` — region-pair contact counts (closest/Cβ/Cα) + bond types
- `canonical_frame_figures` — canonical-frame QC across the Native2026 set
- `pymol_canonical_figures` — ray-traced PyMOL panels (overlay, groove, interface) by class/species
- `example_gil_a02_rs_motif` — GILGFVFTL/HLA-A*02 and the public CDR3β Arg–Ser motif
- `natcompsci2022/` — full reproduction of the Nat Comput Sci 2022 analyses

## Performance

Per-stage timings on a TCR-pMHC complex (1ao7), Apple M3, single thread (`RUN_BENCHMARK=1 pytest
-k benchmark -s` to reproduce):

| stage | time | notes |
|---|---|---|
| parse a gzipped structure | ~19 ms | `.pdb.gz` / `.cif.gz` |
| contact map (5 Å, cKDTree) | ~9 ms | per structure |
| score 1000 candidate peptides | ~8 ms | ~8 µs/peptide (vectorised) |
| **annotate (TCR + MHC), batched** | **~213 ms/structure** | one mmseqs2 call for the whole set; vs ~1.5 s/structure unbatched |
| peak RSS, single-structure pipeline | ~195 MB | |

Annotation is the only network/compute-heavy step and is always **batched** (one mmseqs2 search over
all chains; mmseqs2 parallelises internally — never per-structure, never Python-threaded). Threads are
used only for the embarrassingly-parallel, mmseqs-free stages (structural alignment, write, rendering):
`tcren orient -t N`.

## Tests

```fish
pytest -m "not slow"          # unit + fast regression (the CI gate)
pytest                        # add the arda/mmseqs-backed regression tests
RUN_BENCHMARK=1 pytest -k benchmark -s
```

## Cite

Karnaukhov *et al.* "Predicting TCR specificity from structure with a residue-level statistical
potential." *Nature Computational Science* (2024).
[doi:10.1038/s43588-024-00653-0](https://www.nature.com/articles/s43588-024-00653-0)

> Free for academic and non-commercial use; for commercial use, contact the study's corresponding author.
