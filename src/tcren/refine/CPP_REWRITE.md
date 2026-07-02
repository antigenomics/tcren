# Native C++ rewrite plan — peptide modelling in tcren

Goal: a **self-contained, dependency-free, license-free** peptide-modelling path inside tcren, fast
enough to run over large repertoires, with **no runtime dependency on OpenMM / ProMod3 / PyRosetta**.
Those three are installed (env `tcren-fold`) as **reference oracles**, not shipped dependencies — we
validate the native C++ against them and retire them from the hot path.

## Principle: oracle, don't wrap

Do **not** reimplement what a mature library already does optimally. OpenMM's force evaluation is
world-class C++; PyRosetta's FlexPepDock is the reference refiner. We keep them as the **accuracy
ceiling** and reimplement only the pieces tcren needs to own: the orchestration, the loop-closure /
backbone-sampling kernels, and a compact restrained refinement — each tuned until it approaches the
oracle on the Native2026 benchmark.

Each engine in `engines/` is swappable behind one interface (`Engine`/`ModelResult`), so a native C++
engine drops in beside the reference ones with no caller change. That separation is the "ready" state.

## Part-by-part disposition

| Part | Reference (conda, oracle) | Native tcren target | Status |
|------|---------------------------|---------------------|--------|
| Anchor prediction | — (heuristic) | `refine/anchors.py` (pure stdlib) | ✅ done, no C++ needed |
| Rigid-body refine | — | `src/_refine/refine.cpp` (DOPE MC) | ✅ done |
| CCD loop closure | ProMod3 `loop` | `src/_fold/fold.cpp` (CCD Cα) | ✅ kernel done; upgrade below |
| Full-atom loop build | ProMod3 loopmodel | `_fold`: N–Cα–C φ/ψ chain + KIC + fragment | ⬜ to write |
| Side-chain repack | ProMod3 sidechain / Rosetta packer | `src/_relax/` rotamer packer (SCWRL/dead-end elim) | ⬜ to write |
| Physics minimisation | **OpenMM** (AMBER) | `src/_relax/`: restrained minimiser over a compact energy (DOPE + soft-sphere + anchor harmonic) — NOT a full MD force field | ⬜ to write; OpenMM stays optional-accuracy |
| Flexible-backbone refine | **PyRosetta FlexPepDock** | `src/_relax/`: native Metropolis MC (backbone small/shear + repack + score) | ⬜ to write; PyRosetta = ceiling |

## What each reference oracle is for

- **OpenMM** — ground-truth *energetics*. Validate that the native restrained minimiser relaxes clashes
  and finds the same local minimum basin (compare final peptide RMSD + relative energy ordering).
- **ProMod3** — ground-truth *loop geometry*. Validate the native full-atom loop build (φ/ψ closure +
  peptide-bond regularity) against ProMod3's loopmodel output on the same anchors.
- **PyRosetta FlexPepDock** — ground-truth *refinement accuracy*. It is the accuracy ceiling in
  `fold_benchmark.py`; the native `_relax` engine must approach its RMSD-to-native from a displaced
  start. FlexPepDock is a *protocol inside PyRosetta* (`pyrosetta.rosetta.protocols.flexpep_docking`),
  driven via the Python API (see `oracle_flexpep.py`), not a standalone binary.

## Validation protocol (already scaffolded)

`scripts/fold_benchmark.py` is the harness: displace the native peptide, re-model with each engine,
measure peptide RMSD to the native crystal pose (MHC-groove superposition), and report the oracle
column. The native C++ engines are "ready to ship" when, on the full Native2026 set (n≈374), they
reach within a target Δ of the PyRosetta/OpenMM oracle RMSD while running ≥10× faster with zero heavy
deps. Guard: `RUN_BENCHMARK=1`.

## Reproducing the reference-oracle env (`tcren-fold`)

The three oracles are installed in a **dedicated** env so the primary `tcren-nb` env stays pristine.
Recipe used (osx-arm64, all three have arm64 builds):

```fish
conda create -y -n tcren-fold -c conda-forge -c bioconda \
    python=3.11 pip cmake cxx-compiler mmseqs2 "numpy>=1.26" "scipy>=1.11" "biopython>=1.84" \
    openmm pdbfixer openstructure promod3           # OpenMM 8.5.2, OST 2.11.1, ProMod3 3.6.0
conda run -n tcren-fold pip install -e .            # tcren editable + arda + rapidfuzz + pytest
conda run -n tcren-fold pip install pyrosetta-installer
conda run -n tcren-fold python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"  # FlexPepDock, academic, ~1.5 GB
conda env config vars set -n tcren-fold ARDA_HOME=/Users/mikesh/vcs/code/arda   # arda VDJ db (wheel install can't self-locate it)
```

Validated: `available_engines()` → `['dope','ccd','openmm','promod3']`; `flexpep_refine` (PyRosetta)
runs. `pyrosetta` / `openmm` / `ost` / `promod3` are **oracle-only** — never added to tcren's
`pyproject.toml` dependencies.

## Build note

New C++ kernels follow the existing stdlib-only pybind11 pattern (`_refine`, `_fold`, `_align`): a
3-line `pybind11_add_module` in `CMakeLists.txt`, no `find_package(OpenMM/OST)` linking (that path is
fragile on osx-arm64 and would re-introduce the dependency we are removing). The reference libraries
are called only through their **Python APIs** in the oracle engines, never linked into tcren's exts.
