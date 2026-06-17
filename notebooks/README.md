# tcren notebooks

Runnable examples and analyses on top of the `tcren` library. They share one conda env.

## Environment

The notebooks run in a **separate conda env** (`tcren-nb`) from the lean `tcren` library-dev
env: it is the standard `tcren` environment plus the analysis/viz packages the notebooks need
(Jupyter, matplotlib, rapidfuzz, scikit-learn, logomaker), installed via the `notebooks` extra
in `pyproject.toml`. `arda` (and `tcren` itself) come from that editable install, pinned to
`arda@2.0.1` — no separate checkout.

```fish
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
- `natcompsci2022/` — full reproduction of the Nat Comput Sci 2022 analyses (see its README)
