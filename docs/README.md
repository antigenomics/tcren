# TCRen — structural evidence for T-cell recognition

[![PyPI](https://img.shields.io/pypi/v/tcren)](https://pypi.org/project/tcren/)
[![Tests](https://github.com/antigenomics/tcren/actions/workflows/tests.yml/badge.svg)](https://github.com/antigenomics/tcren/actions/workflows/tests.yml)
[![Documentation](https://github.com/antigenomics/tcren/actions/workflows/docs.yml/badge.svg)](https://docs.isalgo.dev/tcren/)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
[![GPLv3](https://img.shields.io/badge/license-GPLv3-green)](https://github.com/antigenomics/tcren/blob/master/LICENSE)

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/antigenomics/tcren/master/docs/_static/tcren-splash-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/antigenomics/tcren/master/docs/_static/tcren-splash-light.svg">
  <img alt="TCRen contact graphs branch into a statistical-potential route for peptide ranking and measured affinity, response and kinetic associations, and an ipTM–Betti–enrichment route for binder/non-binder and binder/mock ranking. Geometry, visualization and statistics support both routes." src="https://raw.githubusercontent.com/antigenomics/tcren/master/docs/_static/tcren-splash.png" width="100%">
</picture>

TCRen turns TCR–peptide–MHC structures into **contact-based scores, receptor
rankings and interpretable structural readouts**. It brings two approaches into
one Python library and command-line toolkit:

- **Residue statistical potentials** rank candidate epitopes for a fixed TCR and
  support studies of measured binding affinity, response magnitude and kinetics.
  Contact energies, structural alanine references and mechanical readouts retain
  their separate definitions.
- **ipTM + Betti fraction + contact enrichment** rank candidate receptors against
  one pMHC, including binder/non-binder and binder/mock comparisons. The default
  combines generator confidence with global contact organisation and local
  evidence calibrated against frozen mock structures.

Geometry, surface chemistry, contact attribution, molecular visualization and
matched statistics help explain both routes.

## Start with your question

| Question | Route | Worked guide |
|---|---|---|
| Which peptides rank highest for this TCR? | `score` for fixed-template threading; `screen --task epitope` for separately modelled poses | [Epitope ranking](https://docs.isalgo.dev/tcren/ranking.html) |
| Which modelled receptors rank highest for this pMHC? | `screen --task tcr`: ipTM, Betti fraction and enrichment | [Receptor screening](https://docs.isalgo.dev/tcren/screening.html) |
| How do structural scores relate to affinity or response magnitude? | Join measured endpoints, then `evaluate --kind continuous` | [Matched evaluation](https://docs.isalgo.dev/tcren/evaluation.html) |
| How do contacts and mechanics relate to kinetics? | `md`, `mechanics`, then matched endpoint evaluation | [Dynamics and kinetics](https://docs.isalgo.dev/tcren/md-validation.html) |
| Why does a model receive this score? | `features`, `contacts`, `export`, `surface`, `superimpose` | [Inspect a model](https://docs.isalgo.dev/tcren/kit.html) |
| How can I render or compare the evidence? | PyMOL scenes, contact maps, `associate`, `evaluate` | [Visualization gallery](https://docs.isalgo.dev/tcren/gallery.html) |

## Install and rank a receptor cohort

```sh
pip install tcren
tcren build-mhc-ref

tcren screen -s models.tar.gz --task tcr \
  --confidence confidence.tsv -o scores.tsv
```

Supply separately modelled candidate receptors for **one peptide and MHC
background**. `confidence.tsv` contains unique `complex.id`, `iptm` on [0,1] and
`plddt` on [0,100]; IDs match the structure stems. TCRen reads generator
confidence from metadata or this table. Annotation is batched internally, and
`-t 0` selects the allocated CPU budget by default.

Wheels include five native extensions. [ARDA](https://github.com/antigenomics/arda)
provides TCR annotation and its first-use reference/search setup; the MHC
reference is built on demand. The [installation tutorial](https://docs.isalgo.dev/tcren/getting-started.html)
covers raw structures, input formats and first-run setup.

For peptides threaded onto one fixed template:

```sh
tcren score -s template.pdb -c peptides.txt --potential tcren2 -o energies.csv
```

For separately modelled peptides, use `screen --task epitope`. Its structural
alanine reference truncates peptide side chains in three dimensions and rebuilds
contacts before combining the energy difference with ipTM. Fixed-template
threading and separately modelled poses are different protocols.

## Read the receptor score

```text
S = (1 - 0.41*q - 0.22*(1-q))*ipTM + 0.41*q*B + 0.22*(1-q)*M
```

| Component | Structural meaning |
|---|---|
| `ipTM` | Generator interface confidence |
| `B` | Recognition graph cycle fraction, β₁/E, over CDR:peptide, peptide:MHC and CDR:MHC contacts |
| `M` | All-six-CDR/peptide local contact-cycle evidence, calibrated against 2,056 whole mock graphs |
| `q` | Inclusive percentile of unlabelled cohort-mean pLDDT against five frozen reference means |

The weights are nonnegative and sum to one. Mock histograms and compressed
whole-control calibration ship with the package; scoring fits no new parameters.
Missing active values use flagged frozen training imputation. Read `screen_score`,
component values, weights, missingness and availability reasons together.

**Cohort mean remains the default.** Models without TCR:peptide contacts retain a
reason and no combined score; they do not count towards the default 30-model
floor. For one structure, explicitly select a representative confidence context:

```sh
tcren screen -s model.pdb --task tcr --confidence confidence.tsv \
  --confidence-context nearest-plddt --min-cohort 1 -o score.tsv
```

`nearest-mean-plddt` uses a small same-pMHC set's observed mean. Both alternatives
report support and distance to the selected frozen representative. See
[confidence contexts](https://docs.isalgo.dev/tcren/motif-confidence.html).
Supplied MHC names retain their declared resolution, separately from sequence
representatives; exact allele identity is not a term in the default score.

**Interpretation.** These are structural rankings and measured-endpoint
associations. Statistical contact energy is not physical binding free energy;
the combined receptor score is development-selected and is not a calibrated
binding probability. Representative confidence bins approximate cohort context.
Keep matched populations, component readouts and every cohort gain or loss visible.
[Score definitions](https://docs.isalgo.dev/tcren/screening.html) and
[reference dependencies](https://docs.isalgo.dev/tcren/fit-free.html) give the full contracts.

## Explore an application or reproduce a comparison

| Example | What you can reproduce |
|---|---|
| [Confidence-motif score notebook](https://github.com/antigenomics/tcren/blob/master/notebooks/confidence_motif_score.ipynb) | Self-contained frozen-score example, components and confidence contexts |
| [TCRvdb frozen-score reproduction](https://github.com/antigenomics/tcren/blob/master/notebooks/tcrvdb_frozen_score_reproduction.ipynb) | Accepted-score replay and matched comparisons from a supplied retained input bundle |
| [Contacts and interaction types](https://github.com/antigenomics/tcren/blob/master/notebooks/contact_thresholds_and_bondtypes.ipynb) | Atom/residue contacts and chemical interaction labels |
| [Canonical molecular figures](https://github.com/antigenomics/tcren/blob/master/notebooks/pymol_canonical_figures.ipynb) | Reusable PyMOL scenes and canonical views |
| [Surface topology](https://github.com/antigenomics/tcren/blob/master/notebooks/surface_topology.ipynb) | Geometry and receptor-free pMHC surface readouts |
| [CLI reproduction guide](https://docs.isalgo.dev/tcren/reproduction.html) | Raw coordinates → exported components → matched binary or continuous evaluation |

The TCRvdb input bundle is supplied separately and is not distributed with the
package. Frozen-feature replay and raw-coordinate reproduction are distinct
checks. Public structure notebooks document their reference-data prerequisites;
PyMOL examples additionally need the PyMOL executable. For scientific comparisons,
pin the installed artifact and reference identities, preserve original labels and
report exact membership alongside the metrics. [SOURCES](https://github.com/antigenomics/tcren/blob/master/SOURCES)
records the bundled parameter provenance.

## Runtime and batching

One public 1AO7 structure (4,848 coordinate records), exporting 143 feature
columns on an Apple M3 Max with one compute CPU and existing annotation
references, took **7.70 seconds**. The matched pre-optimization observation took
8.65 seconds; outputs were byte-identical. Each timing is one unprofiled run.
Sampled process-tree RSS was 1.12 GB after the change and can double-count shared
pages. [Recorded workload and reference hashes](https://github.com/antigenomics/tcren/blob/master/docs/_static/runtime-1ao7.json)
make the scope explicit.

Pass a whole cohort in one command. Annotation searches and coarse descriptor
workers use separate bounded phases; avoid an outer pool of all-core commands.
Reference acquisition/index construction is separate setup work. See
[performance](https://docs.isalgo.dev/tcren/performance.html) for CPU budgets,
measurement and memory guidance; this single-structure example does not estimate
cohort throughput.

## Documentation and development

- **Learn:** [installation and first run](https://docs.isalgo.dev/tcren/getting-started.html),
  [worked score notebook](https://docs.isalgo.dev/tcren/notebooks/confidence_motif_score.html).
- **Do a task:** [ranking](https://docs.isalgo.dev/tcren/ranking.html),
  [inspection](https://docs.isalgo.dev/tcren/kit.html), [reproduction](https://docs.isalgo.dev/tcren/reproduction.html).
- **Look up a definition:** [score contract](https://docs.isalgo.dev/tcren/screening.html),
  [features](https://docs.isalgo.dev/tcren/features.html), [API](https://docs.isalgo.dev/tcren/modules.html).
- **Understand the evidence:** [reference dependencies](https://docs.isalgo.dev/tcren/fit-free.html),
  [confidence contexts](https://docs.isalgo.dev/tcren/motif-confidence.html).

```sh
./setup.sh --tests
.venv/bin/ruff check src tests
source .venv/bin/activate
make -C docs html SPHINXOPTS=-W
```

[Changelog](https://github.com/antigenomics/tcren/blob/master/docs/CHANGELOG.md) ·
[API skill](https://github.com/antigenomics/tcren/blob/master/skills/tcren/SKILL.md) ·
[Editable splash source](https://github.com/antigenomics/tcren/blob/master/docs/_static/tcren-splash.tex)
