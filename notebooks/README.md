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

## marimo apps

Two notebooks here are [marimo](https://marimo.io/) apps rather than Jupyter notebooks — plain
Python files, reactive, and runnable as small web apps. They need the `marimo` extra
(`pip install 'tcren[marimo]'`):

```bash
marimo run notebooks/pnative_channels.py     # or `marimo edit` to open the cells
marimo run notebooks/surface_topology.py
marimo run notebooks/pymol_interactive.py
```

- `pnative_channels.py` — the released scoring path end to end: one featurisation pass over a
  directory of structures, then one latent-class fit per channel, then their combination into
  `P_native`. It ends on the correlation between the geometry and energetics posteriors, which is
  the sign that flips when a pose was copied from a template rather than produced by the physics.
  Needs the MHC allele reference (`tcren build-mhc-ref`) and the canonical database
  (`tcren fetch-data`).
- `surface_topology.py` — pMHC surface topography: elevation / hydropathy / charge over the groove,
  the per-structure scalars, and the featureless-vs-bulged epitope comparison. `surface_topology.ipynb`
  beside it is **generated** from this file by `make -C docs notebooks`; edit the `.py`, never the
  `.ipynb`.
- `pymol_interactive.py` — a PyMOL render explorer over the canonical scenes (overlay, groove,
  interface, residue importance). PyMOL itself is a separate binary the module shells out to, so it
  is not a Python dependency and must be installed separately.
