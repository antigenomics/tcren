# tcren notebooks

Runnable examples and analyses on top of the `tcren` library. They share one conda env.

## Environment

The notebooks run in a **separate conda env** (`tcren-nb`) from the lean `tcren` library-dev
env: it is the standard `tcren` environment plus the analysis/viz packages the notebooks need
(Jupyter, matplotlib, rapidfuzz, scikit-learn, logomaker), installed via the `notebooks` extra
in `pyproject.toml`. `arda` (and `tcren` itself) come from that editable install, which requires
`arda-mapper >= 2.5.7` — no separate checkout.

```bash
bash setup.sh        # creates/updates the tcren-nb env, editable installs, Jupyter kernel
```

`setup.sh` reads `environment.yml`, editable-installs `tcren[notebooks,viz]`, and registers the
**"Python (tcren-nb)"** kernel. Select that kernel in VS Code / Jupyter. Override `ENV_NAME` if
your layout differs.

## Notebooks

- `complementarity_map_2d` — 2D interface maps, multiple structural + map views of 1ao7
- `contact_thresholds_and_bondtypes` — region-pair contact counts (closest/Cβ/Cα) + bond types
- `canonical_frame_figures` — canonical-frame QC across the Native2026 set
- `pymol_canonical_figures` — ray-traced PyMOL panels (overlay, groove, interface)
- `mhc_pseudosequence_mps` — NetMHCpan MHC pseudosequence (MPS) residues vs. peptide contacts
- `example_gil_a02_rs_motif` — GILGFVFTL/HLA-A*02 and the public CDR3β Arg–Ser motif
- `pocket_cdr_3d` — 3D peptide-binding pocket with the CDR loops overlaid (py3Dmol)
- `tcren_analysis` — potential heatmaps (TCRen / MJ / Keskin) and contact distributions
- `natcompsci2022/` — full reproduction of the Nat Comput Sci 2022 analyses (see its README)

## The score set, end to end

Two notebooks run the whole path — fetch structures, `tcren features`, score — on the two ranking
tasks the library is built for. Both need the MHC allele reference (`tcren build-mhc-ref`, once) and
fetch their structures from the Hugging Face dataset `isalgo/tcren_structures` into `data/`.

- `score_vdjdb_panel` — **receptor ranking for a fixed epitope**. The balanced VDJdb benchmark: 1,089
  TCRmodel2 complexes over 22 epitope cohorts, 523 real binders against 566 mock mispairings. Scored
  with `tcren.score.score_table` and reported **one cohort at a time**, template-covered (6) and
  template-free (16) apart, with the five channel scores per cohort. The featurisation pass is 1,089
  structures in 169 s on 16 cores; the notebook reuses `data/vdjdb_binder_benchmark/features.tsv`
  when it is already there, and prints the command that regenerates it.
- `rank_peptides_cpl` — **peptide ranking for a fixed receptor**, and nothing is fitted in it. The
  combinatorial-peptide-library set: 7 clones, 2,103 peptide-swap models, each clone's library split
  into the measured best and worst halves. Per-clone ROC, the rank correlation against the assay's
  graded activation score, and a whole position × residue response matrix predicted from one
  template with `tcren cpl`. Featurisation asks for the `energetics` family only: 2,103 structures in
  279 s over 14 calls.

## marimo apps

Four notebooks here are [marimo](https://marimo.io/) apps rather than Jupyter notebooks — plain
Python files, reactive, and runnable as small web apps. They need the `marimo` extra
(`pip install 'tcren[marimo]'`):

```bash
marimo run notebooks/surface_topology.py     # or `marimo edit` to open the cells
marimo run notebooks/pymol_interactive.py
marimo run notebooks/confident_negatives.py
marimo run notebooks/potts_contact_map.py
```

- `surface_topology.py` — pMHC surface topography: elevation / hydropathy / charge over the groove,
  the per-structure scalars, and the featureless-vs-bulged epitope comparison. `surface_topology.ipynb`
  beside it is **generated** from this file by `make -C docs notebooks`; edit the `.py`, never the
  `.ipynb`.
- `pymol_interactive.py` — a PyMOL render explorer over the canonical scenes (overlay, groove,
  interface, residue importance). PyMOL itself is a separate binary the module shells out to, so it
  is not a Python dependency and must be installed separately.
- `confident_negatives.py` — reading a generator's confidence together with the coordinates. Move
  the confidence slider to the top of the range and watch the structural reading stay spread while
  the confidence-only one collapses to a single number. Needs a `tcren features` table with the
  geometry, topology and Potts columns, and falls back to the shipped native-crystal reference so it
  runs with no data of your own.
- `potts_contact_map.py` — the predicted contact-frequency map (CDR loop × peptide position) beside
  the contacts the structure actually made, and its collapse onto a peptide residue-importance
  profile. The same call `tcren potts map` makes; needs only a structure file.
