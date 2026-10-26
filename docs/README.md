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
ranking, retain ipTM and the structural components alongside the combined score;
performance depends on the task and evaluation population.

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
# Optional table with unique complex.id, iptm and plddt columns:
tcren screen -s models.tar.gz --confidence confidence.tsv -o scores.tsv
```

A cohort holds one side fixed:

| Task | Fixed side | Varying side | Structural alanine reference |
|---|---|---|---|
| Epitope ranking | receptor | peptide | epitope positions |
| TCR ranking | pMHC | receptor | contacted CDR residues |

Automatic detection reads the annotated peptide and CDR3 pair. Keep one MHC
background and verify experimental identities before joining labels. Automatic mode refuses mixed sets or sets where nothing varies; the explicit
single-model route below declares its task. Missing ipTM is reported; tcren reads it
from metadata or the confidence table and does not calculate generator confidence.

![Task-matched screening workflow](https://raw.githubusercontent.com/antigenomics/tcren/master/docs/_static/screening-workflow.svg)

The structural reference truncates side chains to alanine in three dimensions,
retains backbone and Cβ atoms, and rebuilds contacts. It is not a virtual identity
substitution on a fixed contact map. The selected side follows the ranking task.

## Read the evidence and the score

`screen` reports raw components, standardized contributions, missingness, cohort
flags and a provenance sidecar containing the frozen-reference identity.

- **Confidence:** the generator's `iptm`.
- **Enrichment:** all-six-CDR/peptide contact-cycle evidence against frozen mocks.
- **Energy:** the structural alanine difference on the TCR:peptide interface.
- **Geometry:** height-field Pearson complementarity, `sc_shape`.
- **Topology:** three-interface CDR:peptide + peptide:MHC + CDR:MHC graph cycle fraction,
  `g_recognition_cyclo_frac`; carried for TCR ranking.

TCR ranking uses the frozen confidence-dependent motif score:

```text
S = (1 - 0.41*q - 0.22*(1-q))*ipTM + 0.41*q*B + 0.22*(1-q)*M
```

`B` is the recognition graph cycle fraction; `M` is whole-mock-calibrated
contact enrichment from all six CDR loops and the peptide. `q` is the inclusive
percentile of the unlabelled cohort mean pLDDT against five frozen reference
means. Supply generator `plddt` on the 0–100 scale alongside `iptm`. Missing
values use flagged frozen training imputation. All weights and the score are
bounded in [0,1]; the score is not a binding probability. The packaged compressed
mock reference requires no download or fitting at scoring time.

Cohort-mean confidence remains the default: changing the supplied cohort can
change its weights. Template status is recorded without choosing a different
coefficient. Energy and geometry remain diagnostics; epitope ranking and legacy
custom references retain their existing rules. The preceding two-weight TCR
reference is retained as `screen_reference_tcr_v44.json` for reproduction.

Models with zero TCR:peptide contacts remain in output with a discard reason and
no score. They do not count toward the default floor of 30 models with an
interface. Smaller cohorts report components alone. Flags report by default.
`--on-flag` explicitly selects `report`, `structure`, `geometry`, `withhold` or
`reverse`; responses cannot restore unavailable scores. Reversal is not recommended
as a correction for modelling failures.

For one model, explicitly select a representative pLDDT context:

```sh
tcren screen -s model.pdb --task tcr --confidence confidence.tsv \
  --confidence-context nearest-plddt --min-cohort 1 -o score.tsv
```

Use `nearest-mean-plddt` for a declared small same-pMHC set. These opt-in modes
select the nearest frozen pLDDT representative and report the observed count and
distance. They do not estimate a binding probability or replace cohort scoring.
Missing pLDDT withholds the proxy score. See [confidence contexts](https://docs.isalgo.dev/tcren/motif-confidence.html).

The [worked notebook](https://docs.isalgo.dev/tcren/notebooks/confidence_motif_score.html)
introduces the score on a small reproducible example. The [installation how-to](https://docs.isalgo.dev/tcren/getting-started.html)
walks through a raw-structure run.

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
and reproduction contracts; they are not the current screen. Withdrawn integral commands remain unavailable; explicitly named experimental
contact GLMs have separate reproduction contracts. See [reference dependencies](https://docs.isalgo.dev/tcren/fit-free.html)
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
