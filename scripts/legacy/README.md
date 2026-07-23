# `scripts/legacy/` — Native2026 → `isalgo/vdjdb_structure_models`

Reformat native TCR:pMHC crystal structures (the `Native2026` set) into the packaged format of
the HuggingFace dataset [`isalgo/vdjdb_structure_models`](https://huggingface.co/datasets/isalgo/vdjdb_structure_models),
and stage every Native2026 structure **not already present** as a patch into a local clone of that
dataset. Structure work uses built-in `tcren`; the coordinate / skeleton-plot / contact *format*
is ported from the legacy `tcr-structures-visualization` pipeline.

## What a native entry is (verified against the shipped dataset)

Keyed by `tcr_pmhc_hash` (= `sha256` of 9 VDJdb identity fields, see `reformat.HASH_KEYS`); the
4-char `<pdbid>` is `meta.structure.id`. Per structure:

| Archive (`data/*.tgz`) | File | Source |
|---|---|---|
| `pdb_files_native.tgz` | `aligned_<pdbid>.pdb` | tcren canonical orientation (chains A/B/C/D/E) |
| `coordinates_aa.tgz` | `<pdbid>_aa_coordinates.tsv` | Cα table, chains mapped `A→TCR_alpha`… |
| `contacts_aa.tgz` | `<hash>_<pdbid>_aa_contacts.tsv` | CDR3α/β–peptide residue pairs ≤5 Å |
| `complementarity_maps[_simplified].tgz` | `<hash>[_simplified].svg` | matplotlib skeleton plot (deduped by hash) |
| `vdjdb_structures_metadata.tsv.gz` | 1 row (45 cols, `is_native=True`) | VDJdb join or tcren annotation |

Model-quality columns (`ranking_confidence, plddt, ptm, iptm, tcr_pmhc_iptm`) are blank for
natives; `num_contacts` (tcren `ContactMap` tcr_peptide count) and, if `stcrpy` is installed,
`scanning_angle`/`pitch_angle` are populated.

## Metadata provenance — joinable vs tcren-annotated

- **Joinable** (pdbid matches a VDJdb `meta.structure.id`): the 34 VDJdb columns and the hash come
  from the VDJdb record. tcren annotation alone mints a *different* hash (arda's `j`-allele and MHC
  nomenclature diverge from VDJdb), so joinable structures must take identity from VDJdb.
- **tcren-annotated** (not in VDJdb — the majority of Native2026's missing set): identity fields
  (CDR3 junction, V/J, MHC allele, epitope, species, class) come from tcren/arda; VDJdb-provenance
  columns (`reference.id`, `method.*`, `meta.*`) are left blank. MHC names follow arda nomenclature
  (a `:UniProt` suffix is stripped), not VDJdb's.

Structures without a complete αβ TCR (pMHC-only, or β-only) are skipped and logged.

## Run order

```zsh
source ../../.venv/bin/activate        # tcren venv (+ pandas matplotlib scikit-learn; optional: stcrpy)
python 00_bootstrap.py                 # diff Native2026 vs dataset -> data_dump/worklist.tsv
python 10_validate.py                  # sanity-check the pipeline vs shipped natives (see VALIDATION.md)
python 20_process.py                   # reformat missing structures -> data_patch/
python 25_backfill_metadata.py         # add rows for present natives missing metadata (VDJdb-recoverable)
python 26_complete_orphans.py          # complete metadata-less natives with no recoverable hash (annotate)
python 21_angles.py                    # (optional) fill scanning/pitch angles via STCRpy — needs ANARCI
python 30_assemble.py                  # DRY-RUN: report deltas into ~/hf/vdjdb_structure_models
python 30_assemble.py --apply          # repack archives + append metadata + copy PCA asset + commit
#   then review and: git -C ~/hf/vdjdb_structure_models push     (or 30_assemble.py --apply --push)
```

`20_process.py` flags: `--limit N`, `--only pdbid,...`, `--no-angles`.

Reproduction evidence vs the shipped dataset: **`VALIDATION.md`** (`python 10_validate.py`).

### ANARCI for angles (optional)

`scanning_angle`/`pitch_angle` come from STCRpy, which needs a built ANARCI
([oxpig/ANARCI](https://github.com/oxpig/ANARCI); STCRpy ships the `anarci-mhc` fork). Set it up in
this venv once:

```zsh
brew install hmmer muscle                 # HMMER3 + muscle (ANARCI build deps)
pip install stcrpy                         # pulls anarci-mhc
python -c "from anarci.build_models import build_models; build_models()"   # germlines -> ~/.anarci
python 21_angles.py                        # then fill the angle columns
```

Without ANARCI the angle columns stay blank (many shipped natives also lack them).

## SOURCES

| What | Where | Notes |
|---|---|---|
| Legacy pipeline (coord/plot/contact port) | `/Users/mikesh/vcs/code/tcr-structures-visualization/produce_plots_pipline/` | `coordinates.py`, `plotting.py` imported at run time |
| PCA projection model | `…/produce_plots_pipline/pca_all_structures.sav` | 3D→2D skeleton projection; copied into the dataset by `30_assemble.py` |
| VDJdb-join table | `…/utils_scripts_and_notebooks/vdjdb_structures_annotated.tsv` | 34 VDJdb cols + `TCR_hash`, keyed by `meta.structure.id` |
| Native2026 source PDBs | `~/hf/tcren_structures/Native2026/` | HF dataset `isalgo/tcren_structures`, folder `Native2026` |
| Target dataset clone | `~/hf/vdjdb_structure_models/` | HF dataset `isalgo/vdjdb_structure_models` (SSH remote) |
| Hash formula | `reformat.tcr_hash` | `sha256(cdr3.alpha+v.alpha+j.alpha+cdr3.beta+v.beta+j.beta+mhc.a+mhc.b+antigen.epitope)` |

Paths are overridable via env vars (`TCREN_LEGACY_REPO`, `TCREN_LEGACY_PCA`, `TCREN_VDJDB_ANNOTATED`,
`TCREN_NATIVE2026`, `TCREN_HF_CLONE`). Everything under `data_dump/` and `data_patch/` is gitignored.

## Dependencies

`tcren` (venv) + `pandas`, `matplotlib`, `scikit-learn`, `biopython`. Optional: `stcrpy` (angles;
needs ANARCI — without it `scanning_angle`/`pitch_angle` stay blank).
