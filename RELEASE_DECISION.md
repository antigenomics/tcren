# tcren release and documentation decision

Evidence base: `/Users/mikesh/vcs/code/tcren` (checkout 3.0.0), `/Users/mikesh/vcs/projects/2026-tcren2-code`, `/Users/mikesh/vcs/manuscripts/2026-tcren2-ms`. Every command below was run in this session unless marked otherwise.

---

## 1. Do we need a release?

**Yes. Cut 3.0.0.** Nine of the eleven load-bearing clauses in the manuscript's Code availability paragraph name software that does not exist in the newest published wheel, so a reader who runs `pip install tcren` today gets 2.23.0 and cannot produce a single score-set number the paper reports.

### The gap

| quantity | value | evidence |
|---|---|---|
| checkout version | 3.0.0 | `/Users/mikesh/vcs/code/tcren/pyproject.toml:7` |
| newest git tag | v2.23.0 | `git tag --sort=v:refname \| tail -1` |
| newest PyPI wheel | 2.23.0 | `/Users/mikesh/vcs/code/tcren/STATUS.md:5` |
| unreleased versions between them | 8 (2.24.0, 2.25.0, 2.26.0, 2.27.0, 2.28.0, 2.29.0, 2.30.0, 3.0.0), dated 2026-08-29 to 2026-09-03 | `CHANGELOG.md` headings at lines 438, 351, 286, 259, 231, 169, 112, 6 |
| descriptor catalogue at 2.23.0 vs now | 123 vs 164 in six families | `CHANGELOG.md:112-230`; `docs/descriptor_table.rst` header |

### What the manuscript promises

**Code availability**, `/Users/mikesh/vcs/manuscripts/2026-tcren2-ms/manuscript-latex/main.tex:722-739`, verbatim in the load-bearing parts:

> "TCRen2 is openly available at https://github.com/antigenomics/tcren and on PyPI as tcren, and every score reported in this paper is produced by that released software rather than by analysis code held privately."

> "The potential is re-derived from the deposited crystals, descriptors are emitted from coordinates in one pass, and the read-outs are read off the frozen model, each by its own documented command."

> "The frozen model itself ships inside the package together with the manifest naming the held-out structures it was fitted on, and refitting it from those deposited structures reproduces the shipped fit, which the package's test suite asserts on every release."

> "Executable notebooks ship with the package that run both receptor and peptide benchmarks end to end, from fetching the deposited structures through featurization to the per-cohort read-outs, and they reproduce the values reported here."

> "A descriptor table also carries a digest of the catalogue that wrote it, and the reader refuses a table written under a different one."

**Data availability**, `main.tex:708-721`: structure sets deposited at `isalgo/tcren_structures` with per-dataset retrieval in the dataset catalogue; the template-split receptor set and the HPVGEADYFEY control as separate collections, each carrying labels, receptor junctions, generator confidences and "the TCRen2 scores"; the presentation potential as AAindex MIYS990106 reproduced in all 400 cells; the 2023 mutational archive public; MD trajectories and unfiltered VDJdb arms on request.

### Which clauses are false against 2.23.0

| clause | first version that makes it true | evidence |
|---|---|---|
| the score set (`peptide_score`, `pose_score`, `confidence_residual`, `binder_score`, `channel_scores`) | 3.0.0 | `CHANGELOG.md:12-44` |
| "read off the frozen model, each by its own documented command" (`tcren assess`, `tcren fit-holdout`) | 3.0.0 | `CHANGELOG.md:12-44`; `src/tcren/cli.py` registers 26 commands including `assess` and `fit-holdout` |
| "the frozen model itself ships inside the package together with the manifest" | 3.0.0 | `src/tcren/data/holdout_model.npz`, `src/tcren/data/holdout_manifest.csv.gz`, both `git ls-files`-tracked and not gitignored |
| "a digest of the catalogue that wrote it, and the reader refuses a table written under a different one" | 2.26.0 | `src/tcren/provenance.py:28-42` (`registry_digest`), `:78-107` (`check` raising `StaleTableError`) |
| "descriptors are emitted from coordinates in one pass", at the 164-descriptor catalogue the reported numbers rest on | 2.30.0 | `CHANGELOG.md:112-168` |
| descriptor column names `Phi_*` / `dPhi_*` that the paper's prose uses | 2.26.0 renamed `F_*` to `Phi_*` with no alias | `CHANGELOG.md:286-350` |
| Data availability's "the TCRen2 scores" in each deposited metadata table | 3.0.0, since those columns are score-set read-outs | `main.tex:713-716` |

Two clauses are **false at 3.0.0 as well** and must be fixed in the release or in the sentence, see the checklist:

- "Executable notebooks ship with the package." No notebook is in either distribution format. `pyproject.toml:66` sets `wheel.packages = ["src/tcren"]` and `:76-80` excludes `/notebooks` from the sdist.
- "which the package's test suite asserts on every release." The test that asserts the shipped 7,584-row fit is `@pytest.mark.slow` and gated on `TCREN_HOLDOUT_FEATURES`; CI runs `.venv/bin/pytest -m "not slow"` (`.github/workflows/tests.yml:44`), so it never runs on a release.

### Version number

**3.0.0**, unchanged. It is already in `pyproject.toml:7`, and `publish.yml` validates the pyproject version against the release tag, so the tag is `v3.0.0`. `CHANGELOG.md:12-44` states the major is for the new public surface and the new package data, not for a removal: every 2.30.0 entry point still works and returns what it returned. Nothing argues for renumbering.

### What must land before the tag can be cut

Two hard blockers, three claim-repairs. All are small. See section 2.

---

## 2. Release checklist

Ordered. Each item carries the command that verifies it.

### Blockers

**1. Commit the five untracked paths in the same commit as the tracked edits that reference them.**

```
docs/beyond-the-contact-sum.rst
notebooks/rank_peptides_cpl.ipynb
notebooks/score_vdjdb_panel.ipynb
tests/assets/score/holdout_slice_features.tsv.gz
tests/assets/score/holdout_slice_model.npz
```

`git commit -a` picks up none of them (`git status --porcelain --untracked-files=all | grep '^??'`). All three CI-relevant consequences are real:

- `docs/index.rst:75, :88-89` (tracked and modified) name `beyond-the-contact-sum`, `notebooks/score_vdjdb_panel` and `notebooks/rank_peptides_cpl` in its toctrees. `docs/conf.py:10-19` copies `../notebooks/*.ipynb` into the gitignored `docs/notebooks/` at build time (`.gitignore:310`), so a fresh checkout that lacks the two notebooks has nothing to copy. The docs job runs `sphinx-build -W --keep-going` (`.github/workflows/docs.yml:36`), which turns three missing-document warnings into a red build.
- `tests/unit/test_score.py` opens `tests/assets/score/holdout_slice_features.tsv.gz` with no existence guard.

Verify: `git status --porcelain --untracked-files=all` returns no `??` lines; then `git stash -u && git stash pop` is not needed because a clean clone check is better: `git clone . /tmp/tcren-clean && ls /tmp/tcren-clean/docs/beyond-the-contact-sum.rst /tmp/tcren-clean/notebooks/score_vdjdb_panel.ipynb /tmp/tcren-clean/tests/assets/score/`.

**2. Fix the live `NameError` in the footprint fallback.** `src/tcren/descriptors/compute.py:307` calls `footprint_topology_features(radii)` and `FOOTPRINT_SIZE_FEATURES` inside the `except` branch of `_footprint_columns`; neither is imported in that module, only `footprint_features` is, function-locally at `:303`.

```
src/tcren/descriptors/compute.py:307:30: F821 Undefined name `footprint_topology_features`
src/tcren/descriptors/compute.py:307:67: F821 Undefined name `FOOTPRINT_SIZE_FEATURES`
```
(`/Users/mikesh/vcs/code/tcren/.venv/bin/ruff check src tests --output-format=concise`)

Consequence: `src/tcren/descriptors/table.py:165-166` wraps the family dispatch in one `except Exception` that returns `{"complex.id": id_, "error": ...}`, so a structure whose footprint cannot be built loses **every** family (placement, interface, energetics, potts, kinetics), where the fallback was written to give it NaN topology columns and keep the rest. This is the same defect class, from the same 2.29.0 module split, that the 3.0.0 changelog records fixing for `contact_table` in `annotation/batch.py`.

Verify: `python -c "import tcren.descriptors.compute as m; m._footprint_columns(object())"` returns a dict of NaNs rather than raising.

### Claim repairs

**3. Decide the notebook-shipping question, then make the manuscript sentence and the README agree with the packaging.** Two options, both one line:

- Package them: drop `"/notebooks"` from `pyproject.toml:77` `sdist.exclude`. This puts them in the sdist only; the wheel still packages `src/tcren` alone (`pyproject.toml:66`), so the sentence would still be false for wheel installs.
- Reword: `main.tex:731` to say the notebooks ship with the repository and are rendered in the documentation, and the same for `README.md:675` and `notebooks/README.md:52-62`, which both present the four marimo apps as shipping with the package.

Recommendation: reword. The notebooks are 0.6 MB and 0.4 MB with outputs (`notebooks/score_vdjdb_panel.ipynb`, `notebooks/rank_peptides_cpl.ipynb`), they need HF-fetched structures to re-run, and the docs already render them.

Verify: `python -m build --sdist && tar -tzf dist/tcren-3.0.0.tar.gz | grep -c notebook` and reconcile the count with what the three prose sites claim.

**4. Decide the "asserted on every release" clause.** Either drop `@pytest.mark.slow` from the full-scale contract in `tests/unit/test_score.py` and add a CI job that fetches the hold-out feature table, or narrow `main.tex:729-731` to say the refit contract is asserted on a committed 362-structure slice of the hold-out (232 binders / 130 non-binders, 10 epitopes, 147 descriptors), which is what CI does assert once item 1 lands.

Verify: `.venv/bin/pytest tests/unit/test_score.py -q` under `pytest -m "not slow"` and read which of the two tests ran.

**5. Clear the ruff gate.** `uvx ruff@0.15.22 check src tests` is the exact CI pin (`.github/workflows/tests.yml:58`) and reports 20 errors today, 16 auto-fixable unused imports, one redefinition at `src/tcren/contacts/definitions.py:16`, and the two F821 from item 2. It was already red at HEAD, so the release would otherwise ship against a red signal that is hiding a real bug.

Verify: `/Users/mikesh/vcs/code/tcren/.venv/bin/ruff check src tests` reports 0 errors.

### Pre-tag verification

**6. Fast suite green, with the repo venv first on PATH.**

```
PATH=/Users/mikesh/vcs/code/tcren/.venv/bin:$PATH .venv/bin/pytest -m "not slow" -q
```
Last measured: 774 passed, 3 skipped, 129 deselected, 0 failed, 144.5 s. The 3 skips are py3Dmol missing and two `RUN_BENCHMARK=1` gates.

**Without the PATH fix the same suite reports 2 failed / 772 passed**, both in `tests/regression/test_shipped_potentials.py::test_recipe_reproduces_the_shipped_matrix`. It is an environment artefact, not a source defect: `_cli()` at `tests/regression/test_shipped_potentials.py:37-39` resolves the CLI with `shutil.which("tcren")`, which finds a stale uv-tool editable install at `/Users/mikesh/.local/bin/tcren` whose module map still points `tcren.mechanics` at `src/tcren/mechanics.py`, a file that became a package directory at 2.29.0. Remove or refresh `/Users/mikesh/.local/share/uv/tools/tcren`, or change `_cli()` to prefer `sys.executable -m tcren`.

**7. Version and provenance stamping, checked three ways.**

```
python -c "import tcren; print(tcren.__version__)"                        # 3.0.0
python -c "from tcren.provenance import registry_digest; print(registry_digest())"
python -c "import numpy as np; print(np.load('src/tcren/data/holdout_model.npz', allow_pickle=True)['meta'])"
```
The digest is currently `194108b01c537472f9f94f6526fcb301d45b31ca7847f67cddd1f4656c6808f1`, and `holdout_model.npz`'s meta carries the identical `catalogue_digest`, so the shipped frozen model and the shipped catalogue agree. `src/tcren/__init__.py::_resolve_version` parses `pyproject.toml` at `parents[2]` from a checkout and falls back to `importlib.metadata` from a wheel; both paths must return 3.0.0 after install.

**8. The burn this project has already taken: the installed wheel must not derive its data root from the source-checkout layout.** At 2.12.1 every on-disk data root came from `parents[2]`, which in a wheel resolves to site-packages' parent, so tcren looked for the MHC allele reference there, failed to annotate every structure, **and exited 0** while writing a table with no energy columns. The fix is `paths.tcren_home()` at `src/tcren/paths.py:26-42`: `$TCREN_HOME`, then the source checkout recognised by its `pyproject.toml`, then `$XDG_CACHE_HOME/tcren`. Re-verify it against the built wheel, not the checkout:

```bash
python -m venv /tmp/whl && /tmp/whl/bin/pip install dist/tcren-3.0.0-*.whl
/tmp/whl/bin/python -c "from tcren.paths import tcren_home, data_dir; print(tcren_home(), data_dir())"   # must be ~/.cache/tcren, never site-packages' parent
/tmp/whl/bin/tcren build-mhc-ref
/tmp/whl/bin/tcren info
/tmp/whl/bin/tcren features -s <one structure> -i placement,interface,topology,energetics -o /tmp/f.tsv
/tmp/whl/bin/tcren assess --features /tmp/f.tsv -o /tmp/s.tsv
```
The `assess` step is the one that proves the availability paragraph: it is the command the paper says the read-outs come from, and it must run from a clean wheel install with no checkout above it. Also confirm `tcren scoring` still exits non-zero when every row carries an error, which is the guard that exists because of the silent-success failure above.

**9. Wheel matrix.** `.github/workflows/publish.yml` runs `build-sdist`, then `build-wheels` with `cibuildwheel@v2.23.3` over `ubuntu-latest`, `windows-latest`, `macos-latest` (Apple Silicon arm64) at `CIBW_BUILD: "cp310-* cp311-* cp312-* cp313-*"` and `CIBW_SKIP: "*-musllinux_* *-manylinux_i686 *-win32"`, then `test-wheels` (a smoke import of the C++ extension across 3 operating systems and 4 Python versions), then `publish`. That is 12 wheels plus 1 sdist, the same 13 distributions 2.12.1 shipped.

**Note the workflow dependency graph**: `publish` needs `[build-sdist, build-wheels, test-wheels]` and does **not** depend on `tests.yml` or the lint job. A red ruff gate would not block the publish, which is exactly why items 2 and 5 have to be done deliberately rather than trusted to CI.

**10. Bump the benchmark pin.** `/Users/mikesh/vcs/projects/2026-tcren2-code/requirements.txt:6` reads `tcren==2.23.0`, committed and unmodified, with a comment saying to bump it deliberately and never float it. 2.23.0 has no `tcren.score`, no `tcren assess`, no `tcren fit-holdout`, no hold-out model or manifest, and 123 descriptors against the current 164, so a recompute under this pin cannot produce the score-set numbers the manuscript reports. Bump to `tcren==3.0.0` after the wheel lands.

Verify: `pip install -r requirements.txt && python -c "import tcren; print(tcren.__version__)"` inside the benchmark env.

**11. Tag and release.** `git tag v3.0.0 && git push --tags`, then publish the GitHub release, which triggers `publish.yml`. Confirm afterwards with `pip index versions tcren` and one clean-room run of item 8 against the PyPI wheel rather than the local one.

---

## 3. README trim specification

The README is **734 lines in the working tree** (`wc -l README.md`) and **864 at HEAD 1ac6544**. The 130-line difference is an in-flight trim that is not yet committed.

**The five "beyond the contact sum" essays are already moved.** In the committed README they sat at roughly lines 397 to 598. In the working tree they are extracted to `docs/beyond-the-contact-sum.rst` (286 lines, six sections, untracked) and `README.md:447-456` is a 10-line pointer at the rendered page. Do not re-extract them; commit the existing move. Starting a second trim from the committed 864-line file would revert finished work and orphan a file `docs/index.rst:75` already names.

### Section-by-section

Line ranges from `grep -n '^#\{2,3\} ' README.md` on the working tree.

| lines | section | verdict | reason |
|---|---|---|---|
| 1-33 | header, badges, what it is | KEEP | Scope statement, the three-interface extension, the arda dependency. 33 lines and every claim resolves. |
| 34-52 | What it is evaluated on | KEEP | The five manuscript benchmark blocks with their entry points. This is the paper's spine and the reason the API has its shape. |
| 53-101 | What it does (30-row task table) | TRIM 49 to 35 | Cut 8 rows for surfaces no benchmark script and no end-user path uses: `tcren refine`, `tcren shuffle`, `tcren substitute-tcr`, `tcren energy`, `tcren fetch-recent`, `peptide_stability`, `repack`, `check_register`. Keep the modules; the row is what goes. |
| 102-137 | Install | KEEP | `pip install tcren` plus the required `tcren build-mhc-ref` step, `setup.sh`, and the five pybind11 extensions (`tcren._align`, `_refine`, `_relax`, `_fold`, `_geom`). Nothing here duplicates a docs page. |
| 138-282 | Command line | TRIM 145 to 40, MOVE the rest | Largest block in the file and the most duplicated. `:223-234` re-documents `features` / `recognize` / `assess`, which `:283-412` then documents properly. `:179-196` documents four `tcren score` flags (`--drop-untyped`, `--position-weights`, `--soft`, `--intra-weight`) that no benchmark script uses. Keep one example each for `features`, `assess`, `recognize`, `scoring`, `ddg`, `cpl`, `surface`, `contacts`, `annotate`, `superimpose`, `orient`, `build-mhc-ref`, `fetch-data`. Move the exhaustive flag reference to `docs/getting-started.rst`, whose "Command line" (`:64`) and "Scoring structures" (`:73`) sections are the landing spot. |
| 283-412 | One table per structure: descriptors, energies and the score set | KEEP whole | The current answer. Carries the family table, the five read-outs and the `tcren assess` subsection at `:322`. Every claim resolves against `src/tcren/score/`. |
| 413-446 | Library, core code block | TRIM 34 to 26 | Keep `run_pipeline`, `parse_structure`, `ContactMap`, `score_peptides`, `classify_chains`, `intra_peptide_energy`. Move the `summarize_structure` oracle-facade block to `docs/getting-started.rst`; see section 4. |
| 447-456 | Beyond the contact sum | KEEP as the pointer it already is | The move is done. Names the six instruments in one sentence with real paths: `tcren.potential`, `tcren.mechanics.dynamics`, `tcren.energetics.rotamers`, `tcren.topology.footprint`, `tcren.topology.surface`, `tcren.stacking`. |
| 457-498 | CPL response matrices from one template | TRIM 42 to 25 | Paper-critical (`tcren.cpl`, `response_matrix`, `mutation_effect`, `position_scan`, `equimolar_effect`, `tcren cpl`). Keep the code block and the two-reference-state table; cut the per-position exposition, which `docs/features.rst` carries. |
| 499-507 | Batch inputs, gzip, archives | FOLD into 413-446 | Three lines of `iter_structures` do not need their own heading. |
| 508-521 | Canonical orientation, contacts, docking geometry | KEEP | `tcren.docking.canonicalize_structure` / `superimpose` / `docking_angles`, `tcren.contacts.multi_contacts` / `ContactDefinition`. Paper-used through `tcren features` and an end-user surface. |
| 522-535 | 2D complementarity maps and region-pair contacts | KEEP whole | Named by the author. `tcren.project2d.project_structure` / `residue_markup_table` / `contacts_table` / `region_pair_summary`, `tcren.viz.render_complementarity_map`, `tcren.viz.view_pocket_cdr`. Two of those are also in-stage paper-used via `bench/eda/permref_markup_full.py`. |
| 536-598 | Publication figures | TRIM 63 to 34 | Named by the author. Keep `tcren.viz.pymol.render` / `overlay_scene` / `groove_scene` / `interface_scene`, the axis-gizmo table, and `residue_importance` / `importance_scene`. Cut the axis-convention essay and the duplicated marimo paragraph to `docs/gallery.rst`, which already carries "Reading the axis gizmo" at `:13`. |
| 599-622 | Modules (20 rows) | TRIM 24 to 20 | Drop the `tcren.binder` row: the package is 91 lines exposing one function (`is_real_interface`), and its own docstring records that the fitted binder score was removed at 2.26.0, so the row's description of a pre-energy dock check names a rule nothing calls. **Correct, do not drop, the `tcren.paper` row**: `tcren.paper.bootstrap.fetch_hf_structures` is what both new end-to-end notebooks call to fetch the deposited structures, so it is now in the reproduction path. Reword it as the Hugging Face structure-set bootstrap plus the 2022 reproduction helpers. Also drop `tcren.pipeline / oracle`'s `summarize_structure` mention if `docs/oracle.rst` is deleted. |
| 623-650 | Data | TRIM 28 to 16 | Keep the three folder rows (`Native2022`, `Native2026`, `Canonical2026`) and the `tcren.paths.tcren_home()` paragraph, which is the 2.12.1 lesson and belongs in the README. Cut the `TCRen_potential.csv` provenance detail to `docs/potentials.rst`, which has "Shipped matrices" at `:17`. |
| 651-683 | Notebooks | TRIM 33 to 16 | Lead with the two score-set notebooks, list the rest as one line pointing at the docs tutorial toctree. **Fix `:675`**: the four marimo apps do not ship with the package. |
| 684-711 | Performance | TRIM 28 to 6 | Duplicates `docs/performance.rst` (86 lines) verbatim. Keep the two headline rows (batched annotation at about 0.2 s per structure, and the score set at about 34 microseconds per structure over 1,089 structures) and point at the page. |
| 712-719 | Tests | KEEP | Eight lines, three commands, all correct. |
| 720-726 | Project state | KEEP | Pointers to CHANGELOG, STATUS, ROADMAP, BENCHMARKS. |
| 727-734 | Citing | KEEP | The Nat Comput Sci 2024 citation. |

### Proposed outline

| section | lines |
|---|---|
| Header, badges, what it is | 33 |
| What it is evaluated on (the five blocks) | 19 |
| What it does (22-row task table) | 35 |
| Install | 36 |
| Command line (13 worked invocations, flags in the docs) | 40 |
| One table per structure: descriptors, energies and the score set | 130 |
| Library: core calls and batch inputs | 26 |
| Library: Beyond the contact sum (pointer) | 10 |
| Library: CPL response matrices | 25 |
| Library: canonical orientation, contacts, docking geometry | 14 |
| Library: 2D complementarity maps and region-pair contacts | 14 |
| Library: publication figures | 34 |
| Modules (19 rows) | 20 |
| Data | 16 |
| Notebooks | 16 |
| Performance (pointer) | 6 |
| Tests | 8 |
| Project state | 7 |
| Citing | 8 |
| **target** | **479** |

Range to hold: **470 to 490 lines**, from 734 in the working tree and 864 at HEAD.

### What is NOT deleted from `src/`

This is a scope trim of the README, not a code deletion. Three modules leave the README and stay in the package because something imports them:

- `tcren.refine` (1,355 lines across 13 files including `engines/{ccd,dope,openmm_engine,promod3_engine}.py`): `refine/interface.py::interface_energy` backs the `tcren energy` command, `refine/anchors.py::native_peptide` is imported by `src/tcren/__init__.py:12`, and `refine/register.py::check_register` is in `__all__`.
- `tcren.oracle` (152 lines): `descriptors/compute.py:143` imports its private `_native_peptide`.
- `tcren.shuffle` (161), `tcren.recent` (145), `tcren.analysis` (149), `tcren.torsions` (136), `tcren.scoring_rank` (183): each backs a registered CLI command or is imported by `src/tcren/__init__.py`.

Nine 15-line deprecation shims sit at the top level (`ddg.py`, `dynamics.py`, `footprint.py`, `interface_graph.py`, `pose.py`, `rotamers.py`, `scoring.py`, `stability.py`, `surface.py`; `grep -l 'Deprecated location; moved to' src/tcren/*.py`). The README correctly does not list them and `docs/tcren.rst` correctly does not automodule them. Do not delete them in this release: 15 benchmark files still import through `tcren.ddg`, 7 of them inside a `recompute.sh` stage (`cpl_predict.py`, `cpl_channel_scaling.py`, `cpl_gate_*.py`, `cpl_potts_diag.py`, `affinity_fits.py`). Point those at `tcren.energetics.mutation` first. `OBSOLETE.md:6` already records that 3.0.0 passed without the scheduled removal.

---

## 4. Docs specification

`docs/` is **4,681 lines across 15 `.rst` files** (`wc -l docs/*.rst`), plus `docs/rederivation.md`. The survey brief's "5,415 lines across 16 .rst" counted the Markdown file.

| file | lines | verdict | reason |
|---|---|---|---|
| `docs/index.rst` | 105 | **PATCH** | Toctree names `beyond-the-contact-sum`, `notebooks/score_vdjdb_panel` and `notebooks/rank_peptides_cpl`, all three of which are untracked. `docs/conf.py:10-19` copies `../notebooks/*.ipynb` into gitignored `docs/notebooks/` at build time, so a clean checkout has nothing to copy and `sphinx-build -W` (`docs.yml:36`) fails. Also remove the `oracle` and `kit` lines once those are deleted. |
| `docs/getting-started.rst` | 359 | **PATCH** | The receiving page for the README's CLI flag reference. Its "Command line" (`:64`) and "Scoring structures" (`:73`) sections already exist. Add the four `tcren score` reweighting flags and the `summarize_structure` sentence from `docs/oracle.rst`. Its `:162-180` "What tcren can answer" table duplicates the README's task table nearly verbatim; keep one and make the other a pointer. |
| `docs/assess.rst` | 123 | **KEEP AS IS** | The 3.0.0 score-set page. Sections: Two commands (`:12`), What each column is for (`:27`), The five channels (`:64`), Reproducing the frozen model (`:87`), The predecessor tier (`:102`), API (`:110`). Its 10 autodoc targets under `.. currentmodule:: tcren.score` (`assess.rst:4`) all resolve on `tcren.score`. |
| `docs/features.rst` | 704 | **PATCH** | Rewrite the "Scores" section (`:667-704`) as a three-line pointer to `assess.rst`. Today it calls `S` "The recommended score" and presents `S` / `Q` / `T` as "the scores the method proposes", where since 3.0.0 the shipped answer is the five-read-out score set. Its fourth row documents `strain_z` while stating at `:695` that "2.26.0 removed the ``s_strain`` output"; delete the row. Everything above `:667` is correct and stays, including "The gap between the two faces" (`:224`) and "Core recognition descriptors (40)" (`:363`), which matches `len(RECOGNITION_FEATURES) == 40`. |
| `docs/descriptor_table.rst` | 1,041 | **KEEP AS IS, NEVER HAND-EDIT** | `:1` says it is generated by `scripts/gen_descriptor_table.py`, and `tests/unit/test_descriptor_table.py` asserts it matches the catalogue. 22 percent of docs by line count. |
| `docs/oracle.rst` | 121 | **DELETE** | A whole page for `tcren.oracle.summarize_structure`, a 152-line facade the page itself scopes to the paper notebooks (`:4-8`, the 2022 reproduction). Zero of the 87 in-stage benchmark scripts import it; only `descriptors/compute.py:143` uses its private `_native_peptide`. Fold "returns five ready-to-tabulate frames" into `getting-started.rst` and drop the toctree line. |
| `docs/potentials.rst` | 261 | **KEEP AS IS** | Shipped matrices (`:17`), Reproducing them (`:77`), the AAindex3 resource (`:158`), Splitting a potential (`:214`). This is where the Data availability paragraph's MIYS990106 claim is documented. |
| `docs/potts.rst` | 374 | **KEEP AS IS** | 13 autodoc directives; covers `tcren.potts` narratively, which is why `docs/tcren.rst` has no `automodule` for it. |
| `docs/beyond-the-contact-sum.rst` | 286 | **COMMIT AS IS** | Untracked. Six sections: what a contact potential can and cannot express, peptide conformational stability, side-chain repack, footprint shape, surface topology, ring stacking. It is where the README's five essays went. |
| `docs/reliability.rst` | 84 | **PATCH, absorbing `kit.rst`** | Keep as the single page for the fit-free predecessor tier. Its `:76` heading, "Nothing here is fitted against a binding label", is the tier's defining property and stays. |
| `docs/kit.rst` | 110 | **DELETE, merge into `reliability.rst`** | Both pages open on the same 26.2 percent top-ipTM-decile non-binder rate (`kit.rst:1-12`, `reliability.rst:6-9`), both carry a "why not a cohort-refit posterior" section (`kit.rst:44-52`, `reliability.rst:28`), and `kit.rst:30-32` defers to `assess.rst` and `reliability.rst` for the read-outs anyway. Three pages for one demoted tier makes it look like the headline result. |
| `docs/gallery.rst` | 181 | **KEEP AS IS** | Worked examples of every PyMOL view with images. Receives the axis-convention essay cut from the README; its "Reading the axis gizmo" section is already at `:13`. |
| `docs/performance.rst` | 86 | **KEEP AS IS** | Becomes the README's target for the 28-line duplicate. |
| `docs/modules.rst` | 7 | **KEEP AS IS** | Seven lines, one toctree entry (`tcren`). |
| `docs/tcren.rst` | 839 | **PATCH** | No broken autodoc targets: all 95 `automodule` targets resolve, and all 114 autodoc targets across `docs/` resolve, including the 17 bare names under `.. currentmodule::` in `assess.rst:4` and `reliability.rst:4`. The gap is coverage, not correctness: 31 source modules have no `automodule` entry. Most are shims or private, but `tcren.score`, `tcren.reliability`, `tcren.potts` and `tcren.potential` are paper read-outs. Either add them or add one sentence saying their API reference lives on the narrative pages (`assess.rst`, `reliability.rst`, `potts.rst`, `potentials.rst`). |
| `docs/rederivation.md` | Markdown | **DELETE or convert** | Inert. `docs/conf.py:35-43` lists `sphinx.ext.autodoc`, `autosummary`, `napoleon`, `viewcode`, `githubpages` and `nbsphinx`, with no `myst_parser`, so Sphinx never reads it, and no toctree includes it. |

### Toctree state

- **Currently resolving**: `docs/index.rst:67-99` lists 13 documents plus 11 tutorial notebooks; all 13 `.rst` files exist in the working tree and all 11 `.ipynb` sources exist under `notebooks/`.
- **Will break on a clean checkout**: `beyond-the-contact-sum`, `notebooks/score_vdjdb_panel`, `notebooks/rank_peptides_cpl`, because all three are untracked. This is the docs half of blocker 1.
- **Included by no toctree**: `docs/rederivation.md` only. No `.rst` in `docs/` is orphaned.

---

## 5. Stale-reference kill list

### Already fixed by tonight's uncommitted edits

`git diff --stat` in `/Users/mikesh/vcs/code/tcren` shows 27 tracked files, +1,464 / -800.

| site | what was wrong | fix |
|---|---|---|
| `src/tcren/annotation/batch.py` | `contact_table` left behind by the 2.29.0 module split | import restored |
| `src/tcren/cli.py`, four help strings | `ΔΔF` where the quantity is `ΔΔG`; `tcren score --model` where the command is `tcren assess --model`; `F` / `dF` where 2.26.0 renamed to `Phi` / `dPhi`; the `footprint --score` description | corrected |
| `src/tcren/descriptors/catalogue.py:25` | "35-descriptor recognition vector the frozen recognizers consume" | now "40 descriptors", which matches `len(RECOGNITION_FEATURES) == 40` against `len(descriptors()) == 164` |
| `README.md` five essays at roughly `:397-598` | a library README carrying five method essays | extracted to `docs/beyond-the-contact-sum.rst`, pointer at `README.md:447-456` |
| `STATUS.md`, `ROADMAP.md`, `BENCHMARKS.md`, `OBSOLETE.md`, `DATASETS.md`, `CLAUDE.md`, `skills/tcren/SKILL.md` | version and surface drift | updated in the same change set |

### Remaining: a removed entry point named as live

| site | text | truth |
|---|---|---|
| `src/tcren/descriptors/catalogue.py:131-133` | "Fitted and cohort-relative composites (``p_real``, ``p_real_bn``, ``p_forced``, ``p_bind``, ``q_bind``, ``s_strain``) are listed under ``score``. ... :func:`descriptors` excludes them by default." | 2.26.0 removed the `score` family and the `score` invariance class outright (`CHANGELOG.md:305`). `DESCRIPTORS` has six families totalling 164 (placement 31, interface 26, topology 70, energetics 15, potts 5, kinetics 17), none of them `score`, so there is nothing for `descriptors()` to exclude. This is a docstring a user reads to learn what `descriptors()` returns. |
| `data/Native2026_recognize.README.md:7` | reproduction command `tcren recognize --full --scores -s "$TCREN_DATA_DIR/Native2026" -o Native2026_recognize.tsv` | `--scores` was removed from `tcren recognize` at 2.26.0 (`CHANGELOG.md:300`) with the composites it switched on. The documented command fails on any tcren at or after 2.26.0. The file is git-tracked and is the stated provenance of a committed table; it is excluded from the sdist by the `/data` rule (`pyproject.toml:77`), so it does not ship. |
| `docs/features.rst:667-704`, the `strain_z` row at `:695` | documents `strain_z` as a score-table entry while stating in the same cell that 2.26.0 removed the `s_strain` output | Delete the row with the section rewrite in section 4. |

### Remaining: a claim about packaging that is false

| site | text | truth |
|---|---|---|
| `README.md:675` | "Four **marimo** apps ship alongside them (`pip install 'tcren[marimo]'`, then `marimo run <file>`)" | The `[marimo]` extra installs marimo, matplotlib, py3Dmol and pillow. It delivers no file under `notebooks/`. `pyproject.toml:66` packages `src/tcren` into the wheel and `:77` excludes `/notebooks` from the sdist, so `marimo run notebooks/surface_topology.py` works only from a git clone. |
| `notebooks/README.md:52-62` | the same four apps presented as installable through the extra | Same. |
| `manuscript-latex/main.tex:731` | "Executable notebooks ship with the package that run both receptor and peptide benchmarks end to end" | Same, and this one is in a published availability statement. |

### Remaining: a number that does not check out against its own artefact

| site | text | truth |
|---|---|---|
| `CHANGELOG.md:42`, `BENCHMARKS.md:187`, `ROADMAP.md:56`, `skills/tcren/SKILL.md:696` | "the 8,292 structures the fit used" | `holdout_manifest()` is 8,292 rows (7,029 binders / 1,263 non-binders, 31 epitopes, 832 with null ipTM). `src/tcren/score/fit.py:58` inner-joins manifest to features and `:64-66` drops rows with any non-finite descriptor, and the shipped model's meta records `n_pos = 6,429` and `n_neg = 1,155`, so **7,584 complete-case rows entered the fit** and 708 manifest rows (600 binders, 108 non-binders) never reach it. The accurate form is "the 8,292-row manifest, of which 7,584 complete-case rows entered the fit". The manuscript's "the manifest naming the held-out structures it was fitted on" (`main.tex:729`) inherits the imprecision. |

### Remaining: a stated reason that is no longer the code's behaviour

| site | text | truth |
|---|---|---|
| `pyproject.toml:29` | rapidfuzz is "imported at package-import time (potential/__init__), so core, not optional." | The only rapidfuzz import in `src/` is function-local, at `src/tcren/potential/redundancy.py:36-40`, with a comment saying it is lazy precisely so `import tcren` does not hard-require it. Keeping the dependency core may still be right; the stated reason is wrong. |

### Remaining: cross-repo

| site | text | truth |
|---|---|---|
| `/Users/mikesh/vcs/projects/2026-tcren2-code/requirements.txt:6` | `tcren==2.23.0` | Cannot produce the score set. Item 10 of the checklist. |
| `/Users/mikesh/vcs/manuscripts/2026-tcren2-ms/CLAUDE.md`, "Three repositories, one boundary" | "the checkout is ahead at **2.30.0**, 164 descriptors" | `pyproject.toml:7` reads 3.0.0. The 164-descriptor half is correct. |
| `/Users/mikesh/vcs/manuscripts/2026-tcren2-ms/CLAUDE.md`, "Feature pipeline" | "Two commands per structure set, and they are the only supported path: `tcren features` ... `tcren recognize`", and "`recognize` turns that table into ... `S` with its calibrated $p_{\mathrm{binder}}$" | `p_binder` and the frozen Platt links were removed at 2.28.0 (`docs/reliability.rst:78-81`). `recognize` now emits `Q`, `T` and `S` only. The pair that produces the score set is `tcren features` then `tcren assess`, which `src/tcren/cli.py:1372` documents as the one command to run on a folder of AlphaFold models. A stale pipeline recipe in the manuscript repo is how a wrong command reaches an availability statement. |

### The one that is a bug, not a reference

`src/tcren/descriptors/compute.py:307` (`F821` twice). Covered as blocker 2 in section 2, listed here because it is the second instance of the same 2.29.0 module-split defect the 3.0.0 changelog records fixing for `contact_table`.

---

## 6. Notebooks

### What exists

`ls notebooks/` gives 11 Jupyter notebooks and 4 marimo apps:

| kind | files |
|---|---|
| Jupyter, tracked | `canonical_frame_figures`, `complementarity_map_2d`, `contact_thresholds_and_bondtypes`, `example_gil_a02_rs_motif`, `mhc_pseudosequence_mps`, `pocket_cdr_3d`, `pymol_canonical_figures`, `surface_topology`, `tcren_analysis`, plus `natcompsci2022/` |
| Jupyter, **untracked** | `score_vdjdb_panel.ipynb`, `rank_peptides_cpl.ipynb` |
| marimo apps | `surface_topology.py`, `pymol_interactive.py`, `confident_negatives.py`, `potts_contact_map.py` |

`docs/index.rst:88-99` lists all 11 Jupyter notebooks in its tutorial toctree, and `docs/conf.py:10-19` copies them into gitignored `docs/notebooks/` at build time.

### What the manuscript implies must exist

`main.tex:731-734`: "Executable notebooks ship with the package that run both receptor and peptide benchmarks end to end, from fetching the deposited structures through featurization to the per-cohort read-outs, and they reproduce the values reported here."

That sentence asserts four things. Three are now true and one is not:

| assertion | status |
|---|---|
| notebooks exist that run the receptor benchmark end to end | **true**, `notebooks/score_vdjdb_panel.ipynb` |
| notebooks exist that run the peptide benchmark end to end | **true**, `notebooks/rank_peptides_cpl.ipynb` |
| they go from fetching the deposited structures through featurization to the per-cohort read-outs | **true**, see below |
| they **ship with the package** | **false** for both the wheel and the sdist, at 3.0.0 and every earlier version (`pyproject.toml:66`, `:76-80`) |

### The two notebooks the author asked for

Both exist, both are fully executed, and both do what was asked. They are **untracked in git**, which is the only thing standing between them and being real.

**Scoring the 22 individual VDJdb epitopes: `notebooks/score_vdjdb_panel.ipynb`** (23 cells, 11 code). Its own recorded outputs:

- Cell 3 fetches with `tcren.paper.bootstrap.fetch_hf_structures(DATA, folders=("vdjdb_binder_benchmark",))` and untars, then prints: 1,089 structures on disk, 1,089 metadata rows, 523 binders against 566 mock negatives, **22 epitope cohorts**, 297 of 1,089 rows template-covered.
- Cell 5 shells out to `python -m tcren features -s '<set>/*/*.pdb' -i placement,interface,topology,energetics,potts,kinetics -t 0 -o features.tsv`: 1,089 structures featurised in 169 s, 1,089 rows by 164 descriptors.
- Cell 7 calls `tcren.score.holdout_model()` and `score_table(...)`, printing "model: tcren 3.0.0, 6429 binders / 1155 non-binders over 31 hold-out epitopes" and the ten scored columns (`pose_score`, `binder_score`, five `channel_*`, `peptide_score`, `confidence_residual`, `binder_iptm`).
- Cell 9 reads `holdout_manifest()` and states the two overlaps apart: **0 of 1,089 structures shared with the fit**, 11 of 22 epitopes shared.
- Cell 11 prints one row per epitope for all 22 cohorts, with `n`, `n_pos`, the template flag and every arm.
- Cell 13 gives the per-stratum summary: 6 template-covered cohorts and 16 template-free, with the median per-cohort ROC-AUC and the cleared count for each of ipTM, pLDDT, `pose_score`, `binder_score`, `binder_iptm`, `confidence_residual` and the five channels.

This is exactly the per-cohort, template-split reporting the benchmark manifest specifies, and it is generated by the released command line plus the frozen model, with no label anywhere in the scoring path.

**CPL peptide ranking: `notebooks/rank_peptides_cpl.ipynb`** (21 cells, 10 code). Its own recorded outputs:

- Cell 3 fetches `folders=("cpl",)`: 7 clones (1e6, 4c6, 868, ila1, mel5, mel8, sb27), 14 clone-half directories, **2,103 structures**, 2,102 assay rows carrying the graded activation read-out.
- Cell 5 runs `tcren features -i energetics` per clone-half: 2,103 structures featurised in 279 s, 2,103 rows by 15 descriptors, with the per-clone-half counts (161 / 164 for 1e6, down to 162 / 64 for sb27).
- Cell 7 calls `tcren.score.peptide_score(feats)`, showing `dPhi_tcr_pep` and `dPhi_pep_mhc` beside the read-out.
- Cell 9 makes the point the poly-alanine reference exists for: `Phi_tcr_mhc` spreads by 3.86 to 8.33 energy units within a single clone, because every model carries its own generated pose.
- Cell 11 prints the per-clone ROC-AUC and PR-AUC table, best half against worst half, and the summary: **median per-clone ROC-AUC 0.999 over 7 clones, 7 of 7 above chance, 6 of 7 above 0.90**.

### What has to happen

1. **Commit both** (blocker 1). Until then `docs/index.rst` names two documents a clean checkout does not have, and `sphinx-build -W` fails.
2. **Decide the shipping question** (checklist item 3). Either drop `/notebooks` from `sdist.exclude`, or reword `main.tex:731`, `README.md:675` and `notebooks/README.md:52-62` to say the notebooks ship with the repository and are rendered in the documentation.
3. **One optional tightening.** Both notebooks featurise through the CLI (`python -m tcren features`) but read the score set through the Python API (`tcren.score.score_table`, `tcren.score.peptide_score`) rather than through `tcren assess`. The availability sentence says the read-outs are "read off the frozen model, each by its own documented command". Either add one `tcren assess --features` cell beside the API call, or drop "command" for "documented entry point". The API path is the more readable one in a notebook, so the sentence is the cheaper thing to change.