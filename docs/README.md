<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/antigenomics/tcren/master/assets/tcren_dark.png">
    <img alt="tcren" src="https://raw.githubusercontent.com/antigenomics/tcren/master/assets/tcren_light.png" width="340">
  </picture>
</p>

# tcren — rank and inspect TCR–pMHC models

[![PyPI](https://img.shields.io/pypi/v/tcren)](https://pypi.org/project/tcren/)
[![tests](https://github.com/antigenomics/tcren/actions/workflows/tests.yml/badge.svg)](https://github.com/antigenomics/tcren/actions/workflows/tests.yml)
[![docs](https://github.com/antigenomics/tcren/actions/workflows/docs.yml/badge.svg)](https://docs.isalgo.dev/tcren/)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
[![GPLv3](https://img.shields.io/badge/license-GPLv3-green)](https://github.com/antigenomics/tcren/blob/master/LICENSE)

tcren ranks candidate αβ TCR:peptide–MHC structures and exposes the evidence behind
the ranking: generator confidence, structural contact energy, interface geometry
and contact topology. Start with a cohort of separately modelled candidates, or
use the toolkit to inspect one structure or score peptides on a fixed template.

The scores are dimensionless structural readouts. They are not binding free
energies, calibrated binding probabilities or affinity estimates. For receptor
ranking, keep ipTM primary and inspect the structural components alongside it;
the combination is not an established universal improvement over ipTM.

## Install

```sh
pip install tcren
tcren build-mhc-ref
```

Binary wheels contain five native extensions. The runtime dependency
[arda-mapper](https://github.com/antigenomics/arda) annotates TCR chains and fetches
its reference/search backend on first use. The MHC allele reference is built from
IMGT on demand. See the [installation guide](https://docs.isalgo.dev/tcren/getting-started.html)
for reference data, supported environments and first-run behavior.

For development, use the repository's isolated environment:

```sh
bash setup.sh
source .venv/bin/activate
```

## Start with one ranking task

```sh
tcren screen -s models.tar.gz -o scores.tsv
# Optional table with unique complex.id and iptm columns:
tcren screen -s models.tar.gz --confidence confidence.tsv -o scores.tsv
```

A cohort holds one side fixed:

| Task | Fixed side | Varying side | Structural alanine reference |
|---|---|---|---|
| Epitope ranking | receptor | peptide | epitope positions |
| TCR ranking | pMHC | receptor | contacted CDR residues |

Automatic detection reads the annotated peptide and CDR3 pair. Keep one MHC
background and verify experimental identities before joining labels. A mixed set
or a set where nothing varies is refused. Missing ipTM is reported; tcren reads it
from metadata or the confidence table and does not calculate generator confidence.

![Task-matched screening workflow](https://raw.githubusercontent.com/antigenomics/tcren/master/docs/_static/screening-workflow.svg)

The structural reference truncates side chains to alanine in three dimensions,
retains backbone and Cβ atoms, and rebuilds contacts. It is not a virtual identity
substitution on a fixed contact map. The selected side follows the ranking task.

## Read the evidence and the score

`screen` reports raw components, standardized contributions, missingness, cohort
flags and a provenance sidecar containing the frozen-reference identity.

- **Confidence:** the generator's `iptm`.
- **Energy:** the structural alanine difference on the TCR:peptide interface.
- **Geometry:** height-field Pearson complementarity, `sc_shape`.
- **Topology:** CDR:pMHC graph cycle fraction, `g_cyclo_frac`; carried for TCR ranking.

Components use frozen raw-column means and scales, then a declared direction.
Weighted Stouffer combines the observed standardized components. TCR ranking gives
confidence weight 0.8 and shares 0.2 across the structural components; epitope
ranking gives energy and confidence equal weights. These empirical reference
constants are estimated dependencies, not query-derived parameters. Nothing is
fitted or normalized against the query cohort.

Models with zero TCR:peptide contacts remain in output with a discard reason and
no score. They do not count toward the default floor of 30 models with an
interface. Smaller cohorts report components alone. Flags report by default.
`--on-flag` explicitly selects `report`, `structure`, `geometry`, `withhold` or
`reverse`; responses cannot restore unavailable scores. Reversal is not recommended
as a correction for modelling failures.

The [screening guide](https://docs.isalgo.dev/tcren/screening.html) defines every
component, its units, the combination equations, missingness and output fields.
The [inspection guide](https://docs.isalgo.dev/tcren/kit.html) explains how to read
those outputs. Keep MHC class, epitope and template coverage explicit when making
comparisons. High peptide-ranking accuracy alone does not establish receptor
specificity across epitopes.

## Use the toolkit for a particular quantity

| Need | Entry point |
|---|---|
| Rank from an existing structural epitope-energy/ipTM table | `tcren rank-epitopes` |
| Thread candidate peptides on one fixed structure | `tcren score` |
| Structural alanine difference for candidate poses | `tcren scoring --delta --reference-mode structural` |
| Structural contacted-CDR reference | `tcren export --quantity receptor-reference --reference-mode structural` |
| Individual geometry, topology and energy descriptors | `tcren features` |
| Chain/region annotation and residue contacts | `tcren annotate`, `tcren contacts` |
| Canonical MHC frame | `tcren superimpose` |
| pMHC surface relief, charge and hydropathy | `tcren surface` |
| Pose/confidence diagnostics and matched comparisons | `tcren inspect-models`, `tcren associate` |
| Trajectory and mechanical readouts | `tcren md`, `tcren mechanics` |

Fixed-template threading and separately modelled candidates are distinct protocols.
`features` retains historical fixed-map `dPhi` columns; use the explicit structural
scoring/export routes above for the current alanine reference. The full
[API reference](https://docs.isalgo.dev/tcren/modules.html) documents library calls,
experimental contact layers and supported historical reproduction APIs.

Historical `recognize`, `assess` and fitted Potts readouts retain their dependencies
and reproduction contracts; they are not the current screen. The integral/GLM
commands were withdrawn. See [reference dependencies](https://docs.isalgo.dev/tcren/fit-free.html)
and [SOURCES](https://github.com/antigenomics/tcren/blob/master/SOURCES) for parameter provenance and resource identities.

## Development and verification

```sh
./setup.sh --tests
.venv/bin/ruff check src tests
source .venv/bin/activate
make -C docs html SPHINXOPTS=-W
```

Use `.venv/bin/python -m pytest -m 'not slow'` for the non-slow suite; scientific
benchmarks are separate from software checks. Rebuild the editable install after
native changes. Scientific campaigns, old results and publication figures are
maintained with the benchmark and manuscript repositories. [CHANGELOG.md](https://github.com/antigenomics/tcren/blob/master/docs/CHANGELOG.md)
records software changes. The [API skill](https://github.com/antigenomics/tcren/blob/master/skills/tcren/SKILL.md) owns API guidance.
