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
| Side-chain repack | ProMod3 sidechain / Rosetta packer | `src/_relax/repack` (χ enumeration + DOPE, mean field) | ✅ done 2026-08-17, see below |
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

## Where we actually stand (measured 2026-08-17, Native2026, n = 5, macOS arm64)

Recovery from a 1 Å / 15° rigid displacement, `scripts/fold_benchmark.py --limit 5`:

| engine | backbone RMSD (Å) | ms / structure |
|---|---|---|
| `dope` (native C++) | **0.345** | **159** |
| `ccd` (native C++) | 1.046 | 1 |
| `openmm` | 1.469 | 6284 |
| `promod3` | 1.470 | 2585 |

**Do not read that as a native win.** Two things are wrong with it as an accuracy claim:

1. The task is recovery from a *rigid* displacement, and `dope` is a rigid-body refiner restrained to
   its input. It is being asked to undo exactly the move it is built to undo. OpenMM here is doing
   anchor-restrained local minimisation, which is not that task at all.
2. **`dope` does not build side chains.** Measured on 1ao7: the native crystal peptide has 77 heavy
   atoms, `substitute_peptide` strips it to 44 (N, CA, C, O, CB), and `dope` returns 44 — it never
   rebuilds the other 33. OpenMM and ProMod3 return all 77.

So the all-heavy-atom comparison (`dope` 0.279 Å over 49 atoms vs `openmm` 2.125 Å over 84, mean over
5 structures) is a 44-atom score against a 77-atom score. It says nothing about accuracy, and the
native path cannot approach FlexPepDock on any side-chain-sensitive measure while it declines to
place side chains.

**The gap was exactly the `side-chain repack` row, and `_relax.repack` now closes it** (2026-08-17).

### `_relax.repack` — measured

Like-for-like this time: the **same** wrong-rotamer input (χ1 of every peptide side chain rotated
120°, so a full-atom model whose side chains a predictor placed badly), the **same** 33–42
side-chain atoms, the same crystal reference.

| | peptide side-chain RMSD (Å) | time |
|---|---|---|
| input (wrong χ1) | 4.131 | — |
| **`_relax.repack`** | **2.364** | **6 ms** |
| `openmm` (anchor-restrained minimisation) | 4.133 | 3103 ms |

OpenMM does not move them: 4.133 against 4.131 before. That is not a defect in OpenMM — a
`LocalEnergyMinimizer` is a *local* minimiser, a wrong χ1 sits in a different torsional basin, and
gradient descent cannot cross the barrier between them. Relaxing clashes and re-sampling rotamers
are different operations, and only a discrete packer does the second. So the honest claim is narrow
and strong: the native packer does something the physics minimiser structurally cannot, 500× faster.

Over eight structures with the same perturbation, side-chain RMSD recovers 3.93 Å → 1.66 Å and
all-atom 2.78 Å → 1.20 Å, **8/8 improved**, median 6 ms. 1oga recovers exactly (0.000 Å).

Verified rather than assumed: the kernel reproduces the Python prototype's per-residue energy to
`0.0` (not to a tolerance), and a crystal in gives the crystal back — mean shift 0.06 Å, 2 of 77
atoms moving more than 0.5 Å — because the input conformer is index 0 of every enumeration.

### What is still missing

`repack` rotates the side chains a model **has**. It cannot rebuild ones that are not there, and
`substitute_peptide` strips past Cβ by design, so the 44-vs-77 atom gap on the *substitution* path
needs side-chain **construction** (ideal internal geometry per residue type) — the `Full-atom loop
build` row, not this one. For AlphaFold/TCRmodel output, which is full-atom, `repack` is the whole
answer.

Cost note: the Python prototype was 0.24 s/structure with only 0.06 s in its kernel, so a *one-shot*
repack never needed C++ for speed alone. It needed C++ because a flexible-backbone MC repacks every
cycle, where 0.24 s × 10³ cycles is fatal; at 6 ms that loop is now affordable.

Scale-out: `scripts/fold_benchmark.sbatch` runs the full n ≈ 374 sweep with all oracles on aldan3
(`aldan3 slurm submit scripts/fold_benchmark.sbatch --env tcren-fold --partition medium`).
FlexPepDock is minutes per structure — it burned 21 min of CPU on 6 structures locally without
finishing, which is why the full sweep is a cluster job.

## Reproducing the reference-oracle env (`tcren-fold`)

The three oracles are installed in a **dedicated** env so the primary `tcren-nb` env stays pristine.
Recipe used (osx-arm64, all three have arm64 builds):

```bash
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
