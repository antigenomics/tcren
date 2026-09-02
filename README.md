<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/antigenomics/tcren/master/assets/tcren_dark.png">
    <img alt="tcren" src="https://raw.githubusercontent.com/antigenomics/tcren/master/assets/tcren_light.png" width="340">
  </picture>
</p>

<h1 align="center">tcren — structure-based prediction of TCR–epitope recognition</h1>

<p align="center">
  <a href="https://pypi.org/project/tcren/"><img alt="PyPI" src="https://img.shields.io/pypi/v/tcren"></a>
  <a href="https://github.com/antigenomics/tcren/actions/workflows/tests.yml"><img alt="tests" src="https://github.com/antigenomics/tcren/actions/workflows/tests.yml/badge.svg"></a>
  <a href="https://docs.isalgo.dev/tcren/"><img alt="docs" src="https://github.com/antigenomics/tcren/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-GPLv3-green"></a>
</p>

**TCRen** predicts which epitopes a T-cell receptor recognises from a single TCR–peptide–MHC
structure (experimental or modelled). It extracts the TCR–peptide contact map and scores every
candidate peptide with a **residue-level statistical potential** derived from contact preferences
in TCR:pMHC crystal structures — answering not "what fancy complex can a model draw?" but "is this
binding physically plausible?".

This is a documented, tested, CLI-driven Python library. TCR chains are annotated with the sibling
[`arda`](https://github.com/antigenomics/arda); MHC chains are mapped and the groove partitioned
against a curated reference; structures are oriented into one canonical frame; and the original
contact maps, potential, and scores are reproduced numerically (validated against committed oracles
to floating-point precision).

Where the original tcren scored TCR:peptide contacts alone, this version also scores the TCR:MHC
and peptide:MHC interfaces, which a full picture of TCR:pMHC binding mechanics and any ΔΔG
estimate both need.

## What it is evaluated on

Five benchmark blocks, chosen to cover the problems users actually bring to a modelled TCR:pMHC
complex. They are the backbone of the accompanying manuscript, and they are why the API looks the
way it does.

| block | the question | entry points |
|---|---|---|
| **Combinatorial peptide libraries** | which peptides does this receptor read? | `tcren pipeline --delta`, `tcren.ddg.neoantigen_ddg` |
| **A functionally validated repertoire screen** | which receptors read this epitope? | `tcren features`, `tcren recognize` |
| **A balanced epitope panel, template-stratified** | the same, where no related complex has been solved | as above, with template availability reported rather than inferred |
| **Molecular dynamics with measured kinetics** | may a single static structure be scored at all, and what does its energy reach? | the three interface energies, `tcren.potts` contact marginals |
| **Model-confidence diagnostics** | which confidently modelled complexes are not real? | `tcren assess`, `tcren.score.confidence_residual`, `tcren.reliability.af_band` |

The last is the one most users reach for first: you already have an AlphaFold model and want to know
whether to trust it. Everything in the score set is frozen on a hold-out that ships inside the wheel
and refits from a manifest that ships with it, so nothing is estimated from the rows you score and no
number depends on what was scored beside it.

## What it does

> **Scope.** `tcren` is for **αβ TCR : peptide–MHC** complexes and nothing else: a TRA and a
> TRB chain, a peptide of standard amino acids, and a class I or class II groove. γδ
> receptors, single-chain constructs, pMHC with no TCR and non-peptidic ligands are out of
> scope — `derive-potential` drops them rather than deriving from them, and there is no flag
> to widen it.

From one TCR–peptide–MHC structure (crystal or model), each task is one command or one call:

| task | command | library |
|---|---|---|
| Score candidate epitopes for a TCR | `tcren score` | `score_peptides` |
| **Rank peptides for a fixed receptor** — the poly-alanine-referenced energy | `tcren assess --peptide` | `score.peptide_score` |
| Percentile-rank a peptide vs background | `tcren rank` | `percentile_rank` |
| ΔΔG of mutations (alanine scan / neoantigen) | `tcren ddg` | `alanine_scan`, `neoantigen_ddg` |
| **Predict a CPL response matrix from a template** | `tcren cpl` | `response_matrix`, `mutation_effect`, `position_scan`, `equimolar_effect` |
| **Binder vs non-binder for a TCR model** | `tcren features` + `tcren assess` | `score.binder_score`, `score.score_table` |
| **Which part of the structure says so** — the five channels | `tcren assess` | `score.channel_scores` |
| **Is the reported confidence warranted?** | `tcren assess` | `score.confidence_residual` |
| **Every interface descriptor, 164 in six families (four by default)** | `tcren features` | `recognition_table(include=...)`, `descriptors` |
| **All interface descriptors, one row per structure** | `tcren recognize` | `recognition_features` |
| **`S`** — geometry, footprint shape and energy in native-sd units | `tcren recognize --features` | `reliability.s_score` |
| **Is *this* model worth believing?** — the score set, `S` and the generator diagnostic | `tcren assess` | `score.score_table`, `reliability.s_score`, `af_band` |
| Refit the frozen model from its shipped manifest | `tcren fit-holdout` | `score.holdout_model`, `score.holdout_manifest` |
| **The contact map as a probability model** — energy, partition function, per-pair contact probability | `tcren potts fit` / `score` / `contacts` | `potts.fit_potts`, `score_sites`, `contact_probabilities` |
| Three-interface energy Φ, poly-Ala ΔΦ, interface geometry | `tcren scoring` | `run_pipeline` |
| Annotate chains + region markup | `tcren annotate` | `classify_chains`, `annotate_mhc` |
| Interface contact table (5/8/12 Å) | `tcren contacts` | `ContactMap`, `multi_contacts` |
| Orient into the canonical MHC frame | `tcren superimpose` / `orient` | `superimpose`, `canonicalize_structure` |
| Graft a TCR onto another pMHC (chimera) | `tcren substitute-tcr` | `substitute_tcr` |
| Wrong-TCR decoy set (recognition negatives) | `tcren shuffle` | `make_decoys`, `graft_tcr` |
| Substitute a peptide + refine its pose | `tcren refine` | `substitute_peptide`, `refine_peptide` |
| **Surface topology of the pMHC face — is this epitope featureless?** | `tcren surface` | `surface_map`, `surface_stats`, `surface_distance` |
| **Backbone dynamics — does the peptide hold its TCR-facing conformation?** | — | `peptide_stability`, `stability_table` |
| Repack side chains into their preferred rotamers | `tcren refine --repack` | `repack` |
| DOPE interface energy (ΔΔG `e_native`) | `tcren energy` | `interface_energy` |
| Interface mechanics — koff proxies (stiffness / rupture) | `tcren recognize --mechanics`, or `tcren mechanics` alone | `interface_mechanics` |
| Re-derive the statistical potential | `tcren derive-potential` | `derive_tcren` |
| Steric-clash / wrong-register QC | — | `interface_clashes`, `check_register` |
| 2D complementarity map + 3D pocket/CDR view | — | `render_complementarity_map`, `view_pocket_cdr` |
| **Publication PyMOL figures, with a labelled axis gizmo** | — | `viz.pymol.render`, `overlay_scene`, `groove_scene`, `interface_scene` |

**Scope — ranking, not affinity.** TCRen ranks peptide/TCR *specificity* for a given receptor (and the
`ddg` matrix is a fast triage, not a free energy). It is **not** an affinity model: on the ATLAS SPR
benchmark neither the raw contact energy nor its poly-alanine difference predicts Kd/ΔG/koff/kon
(|ρ|≤0.3). The one affinity-adjacent quantity a structure predicts is the off-rate koff, via interface
mechanics (`tcren mechanics`) — not the contact sum.

## Install

```bash
pip install tcren          # from PyPI — binary wheels ship the C++ extension; pulls in arda-mapper
tcren build-mhc-ref        # once: builds the MHC allele reference from IMGT (not bundled in the wheel)
```

**`tcren build-mhc-ref` is a required one-time step after a `pip install`.** The curated MHC
allele reference is built from IMGT on demand rather than shipped in the wheel, and every command
that annotates a structure needs it.

For development (a repo-local `.venv` via [`uv`](https://docs.astral.sh/uv/), an editable
install, and the reference data fetched into `data/`):

```bash
bash setup.sh                    # uv venv + editable install + arda + fetch data/ (no conda)
source .venv/bin/activate
```

`setup.sh` needs only `uv` and a C++ compiler (macOS: `xcode-select --install`); it never
touches conda. Pass `--tests` to run the fast suite after install.

tcren ships five small **pybind11/C++ extensions**, built on install by `scikit-build-core`
(which fetches `cmake`+`ninja` automatically): `tcren._align` (MHC-pseudosequence fitting
alignment; a Biopython fallback runs if unbuilt), `tcren._refine` (DOPE atom-level Monte-Carlo
peptide refinement), `tcren._relax` (DOPE interface energy for `tcren energy` / ΔΔG),
`tcren._fold` (CCD loop closure) and `tcren._geom` (interface geometry, clash detection and
contact stability). TCR
annotation is provided by [`arda`](https://github.com/antigenomics/arda), a runtime dependency
published to PyPI as [`arda-mapper`](https://pypi.org/project/arda-mapper/) (it imports as
`arda`); `uv`/`setup.sh` pull it automatically, and from `arda-mapper >= 2.5.7` it auto-fetches
both its own reference **and a static `mmseqs2` binary** on first use — so no conda/bioconda and
no `ARDA_HOME` to set (override the binary with `$ARDA_MMSEQS`). `setup.sh` also runs `tcren
fetch-data` to populate `data/` with the reference structure sets (`Native2026`, `Canonical2026`)
used by `orient`/`superimpose` (set `TCREN_NO_FETCH=1` to skip).

## Command line

```bash
# Score structures: the three interface contact energies (TCRen for TCR↔peptide, MJ for
# TCR↔MHC and peptide↔MHC) and their total Φ. One row per structure.
tcren scoring -s complex.pdb.gz -o scores.csv

# Inputs: a file, a directory, a .tar.gz, a quoted glob, a .txt manifest (one path per line),
# a comma-separated list, or a repeated -s. Mix freely.
tcren scoring -s a.pdb.gz -s b.pdb.gz -o scores.csv
tcren scoring -s 'models/*.pdb.gz' -o scores.csv
tcren scoring -s models/ --delta --geometry -t 8 -o scores.csv   # a directory, 8 workers
tcren scoring -s models.txt -o scores.csv

# --delta adds the poly-alanine reference ΔΦ per interface (ΔΦ_TCR:MHC is identically 0).
# Use ΔΦ, not Φ, when each candidate carries its OWN generated pose: raw Φ then partly reads
# the pose the predictor chose rather than the peptide.
tcren scoring -s 'models/*.pdb.gz' --delta -o scores.csv

# --geometry adds the interface descriptors and Q, the directional decorrelated
# interface-quality score (native-crystal calibrated, so it is defined for a single structure).
tcren scoring -s complex.pdb.gz --delta --geometry -o scores.csv

# Configurable per-interface potential: swap a bundled name (tcren2|karnaukhov2022|mj|keskin),
# a CSV, or None for any interface; default reproduces the built-in per-interface families exactly.
tcren scoring -s complex.pdb -o scores.csv --tcr-mhc-potential keskin

# Opt-in TCR framework regions: --regions {all,cdr,cdr+fr} chooses which TCR regions
# contribute on the TCR side (cdr = CDR1-3 only; cdr+fr adds FR1-3; all = unfiltered, default).
tcren score -s complex.pdb -c candidates.txt -o ranked.csv --regions cdr+fr

# Surface topology: the pMHC face a TCR meets, BEFORE any TCR is there. A height field over
# the groove with hydropathy and charge painted on, per structure, plus the scalars that make
# "featureless" a number: relief, peak_to_valley, frac_above_ridge.
tcren surface -s complex.pdb -o surface.csv

# --compare writes the pairwise map distance (SURFMAP's Manhattan metric) so epitopes cluster;
# --svg writes one figure per structure; --cells writes the long per-cell table.
tcren surface -s models/ -o surface.csv --compare dist.csv --svg figs/ --channel h
tcren surface -s models/ -o surface.csv --channel phobic --scale kd   # or --scale mj

# Three opt-in reweightings of the SAME energy sum, all off by default so nothing moves
# unless asked:
#   --drop-untyped      ignore contacts that are only proximity (no h-bond / salt bridge /
#                       stacking / hydrophobic / polar chemistry)
#   --position-weights  weight by where a contact sits on the peptide (central | tcr_facing);
#                       a clash at an anchor the TCR never touches is not a clash at P5
#   --soft              replace the hard 5 A cutoff with a contact PROBABILITY averaged over
#                       side-chain rotamers, Boltzmann-weighted under DOPE
tcren score -s complex.pdb -c candidates.txt -o ranked.csv --drop-untyped
tcren score -s complex.pdb -c candidates.txt -o ranked.csv --position-weights central
tcren score -s complex.pdb -c candidates.txt -o ranked.csv --soft

# Opt-in intra-peptide term: every interface energy sums over contacts between two DIFFERENT
# chains, so a candidate held in the template's conformation by its own side chains costs the
# same as one that is not. --intra-weight w adds score = Φ + w·E_intra (5 Å, |i-j| >= 3, MJ).
# Sparse by design: an extended class-I 9-mer makes zero to two internal contacts. w=0 = off.
tcren score -s complex.pdb -c candidates.txt -o ranked.csv --intra-weight 0.5
tcren scoring -s complex.pdb -o scores.csv --intra-weight 0.5   # reports Phi_pep_int separately

# Percentile-rank the native (or candidate) peptide's TCRen energy against a random pMHC
# background — small rank_pct = the peptide scores among the best binders.
tcren rank -s complex.pdb -o rank.csv

# Fast ΔΔG of peptide point mutations (virtual-matrix path: no atoms move, no re-docking).
# Requires --native (the peptide) and exactly one mode: --alanine-scan or --mutant.
# ddG = E(native) - E(mutant), and lower energy binds better, so POSITIVE = stabilising.
tcren ddg -s complex.pdb --native EPITOPE --alanine-scan -o ddg.csv

# Predict a combinatorial-peptide-library (CPL) response matrix from ONE template TCR:pMHC
# structure: every peptide position x all 20 residues, threaded on the template's own contact map.
# Every cell sums BOTH peptide-bearing interfaces (TCRen over TCR:peptide + Miyazawa-Jernigan over
# peptide:MHC), because the assay reads activation, which needs presentation as well as engagement.
tcren cpl -s complex.pdb -o cpl_matrix.csv
# Two reference states, both emitted, and a cell means nothing except against one of them:
#   effect_equimolar  vs the 1/20 mixture  -> the CPL background; compare against a measured matrix
#   effect_wild_type  vs the template residue -> the mutation-scan / neoantigen question
# Positive is favourable on both. Three narrower questions off the same matrix:
tcren cpl -s complex.pdb --position 5                  # every substitution at position 5, best first
tcren cpl -s complex.pdb --position 5 --mutation W     # just that one cell
tcren cpl -s complex.pdb --position 5 --to-mixture     # cost of giving position 5 up to the mixture

# Rank candidate receptors against a fixed pMHC. `assess` is the score set: is the pose real, is it
# a binder, and which part of the structure says so. Every read-out is defined for ONE structure --
# the transform, the class means and the covariance are frozen on a hold-out that ships in the wheel.
tcren features -s candidates/ -o feats.tsv
tcren assess   --features feats.tsv -o scores.tsv
tcren assess   --features cpl_feats.tsv --peptide -o cpl.tsv   # when the PEPTIDE is what varies

# The fit-free predecessor tier on the same table: Q (interface geometry), T (footprint shape) and
# their composition with the contact energy, S. Still shipped, still reported, and it COMPOSES with
# the score set rather than being replaced by it.
tcren recognize --features feats.tsv -o qts.tsv

# One TSV per structure: every interface descriptor (geometry + energies).
tcren recognize -s my_pdbs/ -o recognize.tsv          # descriptors, one row per PDB

# End-to-end candidate-epitope scoring from a structure
tcren score -s complex.pdb -c candidates.txt -o ranked.csv

# Wrong-TCR decoys: keep each ORIENTED complex's pMHC, graft on 10 other complexes' TCRs (within
# MHC class, no real pairing). Real-vs-decoy trains a label-free TCR-recognition classifier.
tcren orient -s natives/ -o oriented/          # inputs must share the canonical MHC frame
tcren shuffle -s oriented/ -o shuffled/ --n 10

# Substitute a peptide and refine its pose (knowledge-based MC scored by the DOPE atom-level
# statistical potential — independent of the TCRen/MJ scoring potentials, restrained to the input).
# Not physics relaxation — use Rosetta FlexPepDock for that.
tcren refine -s complex.pdb -o refined/ --substitute KQWLVWLFL

# Structures: any of .pdb / .cif / .pdb.gz / .cif.gz, a directory, or a .tar.gz batch
tcren contacts -s batch.tar.gz -o contacts.csv --interface tcr_peptide

# Per-residue markup: TCR (CDR/FR) + MHC groove (helix/floor) + peptide in one table.
# --regions all|tcr|mhc|peptide filters; --pseudo also marks NetMHCpan groove residues (MPS).
tcren annotate -s complex.cif.gz -o markup.csv --regions mhc --pseudo

# Superimpose structure(s) onto the canonical frame, by MHC, against the canonical database
# (data/Canonical2026, fetched at install). Detects MHC class + species and averages the
# superposition over every database structure of that class/species. Chains -> A=Vα B=Vβ
# C=peptide D=MHCα E=MHCβ/β2m. -s takes a file / directory / .tar.gz / glob; -o is a directory,
# or a single structure file (one input) whose extension must match --mmCIF/--compress; -t threads.
tcren superimpose -s complex.pdb -o oriented.pdb           # single file
tcren superimpose -s 'data/*.pdb' -o oriented/ -t 8        # glob -> directory, threaded

# Build a canonical database from native complexes (how Canonical2026 is produced). Annotation
# is one batched mmseqs call; -t threads only the structural alignment + write.
tcren orient -s data/Native2026 -o data/Canonical2026 -t 8

# Structure outputs are plain .pdb by default; add --mmCIF for .cif and --compress for .gz.
tcren superimpose -s complex.pdb -o oriented/ --mmCIF --compress   # -> oriented/<id>.cif.gz

# Fetch recent TCR-pMHC structures from RCSB -> data/pdb_recent (mmCIF .cif.gz, 5-chain validated)
tcren fetch-recent --discover --after 2024-01-01

# Build the MHC reference once (IMGT/HLA + mouse H-2; cached, not committed)
tcren build-mhc-ref

tcren info
tcren --install-completion        # shell tab-completion (bash/zsh)
```

`tcren orient` and `tcren superimpose` need the reference sets in `data/` (`Native2026`,
`Canonical2026`); `setup.sh` fetches them at install via `tcren fetch-data` (re-run it any time).

## One table per structure: descriptors, energies and the score set

Two commands, two jobs. **`tcren features` reads structures and writes descriptors**;
**`tcren recognize` turns descriptors into scores.** The feature pass is the expensive half, so it
runs once and the scoring pass can be repeated for nothing.

```bash
tcren features  -s my_pdbs/ -o feats.tsv                   # the four default families (--all adds kinetics)
tcren features  -s my_pdbs/ -o shape.tsv -i topology       # one family -- and only it is computed
tcren recognize --features feats.tsv -o scores.tsv         # Q, T, S
```

The 164 descriptors are catalogued in six **families**, four of them computed by default (`potts`
and `kinetics` are opt-in), split by what each is invariant under — which is also the axis along
which they carry independent evidence:

| family | what it is | invariance |
|---|---|---|
| `placement` | where the receptor sits in the groove frame — angles, TCRdock parameters, ride height, shift, offset, the CDR3 loop frames | frame-**dependent** |
| `interface` | how much contact and of what chemical kind — buried area, contact counts and types, hydrogen bonds, clashes | SE(3)-invariant |
| `topology` | the **shape** of the contact set, free of its size — coverage entropy, Hill numbers, Betti numbers, persistence entropy, canonical preference | SE(3)-invariant |
| `energetics` | statistical-potential interface energies `Phi` and their references `dPhi` — poly-alanine, and the smoothed background form | SE(3)-invariant |
| `potts` | the contact map's own energy against a Boltzmann distribution over it — `neg_energy`, `log_z`, `log_lik` (off unless asked) | SE(3)-invariant |
| `kinetics` | the interface as a spring network — stiffness, rupture, coupling residues (off unless asked) | — |

`tcren recognize -s my_pdbs/` reads the structures itself, skipping the feature file:

```bash
tcren recognize -s my_pdbs/ -o recognize.tsv               # descriptors only
tcren recognize -s my_pdbs/ -o scored.tsv --mechanics      # + the spring-network kinetics terms
```

| what you want | columns in `recognize.tsv` |
|---|---|
| **(a) energy** — `Phi` per interface (TCRen on TCR:peptide, MJ on presentation) + references `dPhi` + loop parts | `Phi_tcr_pep`, `Phi_tcr_mhc`, `Phi_pep_mhc`, `dPhi_tcr_pep`, `dPhi_pep_mhc`, `Phi_cdr12`, `Phi_cdr3a`, `Phi_cdr3b`, `dPhi_{pep,tcr,tra,trb}_soft`, `varPhi_{pep,tcr}_soft` |
| **(a′) intra-peptide** (`--full`) — the peptide's contacts with *itself*, which every interface sum omits | `Phi_pep_int`, `n_pep_int` |
| **(b) geometry** — every docking + interface descriptor | `pitch`, `crossing`, `crossing_signed`, `dock_d`, `dock_torsion`, `dock_{tcr,mhc}_u{y,z}`, `extent`, `chain_balance`, `burial`, `n_contacts_{tp,tm}`, `n_pep_contacted`, `ct_{tp,tm}_*` |
| **(c) scores** — written by `tcren recognize --features`, not by `-s`. No training set and no binding label enters any of them | `Q` — interface quality; `T` — footprint shape; `S` — the blocks combined with the contact energy. See [`tcren.cohort`](src/tcren/cohort.py), [`tcren.reliability`](src/tcren/reliability.py). The two-class read-outs are `tcren assess`, below |

### Is *this* model worth believing? — `tcren assess`

A co-folding model will seat any TCR against any peptide, binding or not. `assess` reads the
coordinates it produced and answers four separate questions about them, and **every answer is defined
for a single structure**: the transform, the class means and the covariance are all frozen on a
hold-out that ships with the package, so nothing is estimated from the rows you pass and a score does
not move depending on what was scored beside it.

```bash
tcren features -s models/ -o feats.tsv                 # the expensive pass, once
tcren assess   --features feats.tsv -o scores.tsv      # arithmetic over that table
```

**The score set — five read-outs of one frozen object.** Higher is better throughout.

| read-out | tier | what is estimated | what it answers |
|---|---|---|---|
| `peptide_score` | 0 | nothing; the direction is fixed by the potential | which peptide does this receptor read? |
| `pose_score` | 1 | a covariance over hold-out binders — **no negative, no label** | is this the kind of interface real complexes make? |
| `confidence_residual` | 1 | the same covariance, read as a conditional mean | is the reported confidence warranted? |
| `binder_score` | 2 | class means and covariances, from hold-out binder labels | binder or not? |
| `channel_scores` | 2 | the same object, marginalized to one descriptor family | which part of the structure says so? |

`binder_iptm` is `binder_score + logit(ipTM)`: two log-odds added, no coefficient to fit, and still
defined for one structure. It is the recommended read when a confidence is available.

**The five channels** are named in physics and geometry terms — `placement` (where the receptor sits
in the groove frame), `interface` (how much interface it makes, of what chemistry), `shape` (the
footprint free of its size), `energetics` (the contact chemistry in kT), `mechanics` (the interface as
a network of breakable springs). A marginal of a Gaussian is a sub-block of its covariance — exact,
closed form, no re-fit — so attributing a score to a part of the structure costs an index and nothing
else. They do not sum to `binder_score` and should not: the whole model also reads the correlations
*between* channels. Sometimes a channel is the better instrument: on template-free cohorts of the
22-cohort VDJdb panel `channel_shape` reads 0.637 median ROC-AUC against the full posterior's 0.615.

Pass `--peptide` when the **peptide** is what varies across the structures being compared, as in a
combinatorial library or a mutational scan. Otherwise the five descriptors computed without the
receptor are marginalized out, because they are constant across every structure of one epitope on one
allele and a model reading them reaches the cohort's name without reading an interface.

`assess` also emits, on the same rows: the fit-free predecessor tier `S`; the rank and percentile
within the set with the expected mean score at a recall budget; and — when the table carries ipTM —
the generator diagnostic, which band the model falls in, how often models in that band turned out to
be non-binders, and what `S` still separates *inside* it. On the balanced 22-cohort VDJdb panel the
top ipTM decile is 26.2% [18.7, 35.5] non-binders.

**The coefficients are frozen, and the inputs they were frozen against are named.** That is the
contract the withdrawn cohort-refit posterior could not offer:

```bash
tcren fetch-data                                       # the structure sets the manifest names
tcren features -s <those structures> -o hold.tsv
tcren fit-holdout --features hold.tsv -o refit.npz     # matches the shipped model bit for bit
```

`tcren.score.holdout_manifest()` returns the 8,292 structures with their dataset, epitope, label and
ipTM. From Python the whole set is one call:

```python
import polars as pl
from tcren import score_table

scores = score_table(pl.read_csv("feats.tsv", separator="\t"))
```

**(c) physics of the interaction.** The koff proxies fold into the same table with `--mechanics`;
only the mutation scan, which is per-residue rather than per-structure, needs its own command:

```bash
tcren recognize -s models/ --mechanics -t 0 -o out.tsv    # every per-structure descriptor, one table
tcren ddg       -s complex.pdb -o ddg.csv     # per-residue alanine / neoantigen ΔΔG (fast virtual matrix)
```

`--mechanics` is how to ask for the stiffness tensor, steered rupture and coupling residues on a
cohort. `tcren mechanics` still exists and gives the same numbers, but as a second command it
repeats the parse and both mmseqs searches to return a second table — CSV, keyed `pdb.id` rather
than `complex.id` — that then has to be joined. Inside `recognize` the structures are already
annotated, so the flag costs only the mechanics arithmetic (12 crystals: 19.0 s → 19.5 s, against
22.5 s for the two commands).

(Per the affinity scope caveat above, structures predict the **off-rate koff** via the mechanics
columns, not Kd/ΔG/kon.) From Python:

```python
from tcren.recognition import recognition_features
from tcren.reliability import s_score

feats = recognition_features("complex.pdb")    # dict of the 40 descriptors (RECOGNITION_FEATURES)
score = s_score({k: [v] for k, v in feats.items()})[0]     # one structure is enough
```

## Library

```python
from tcren import run_pipeline, parse_structure, import_structure, ContactMap, score_peptides
from tcren.annotation import classify_chains
from tcren.potential import tcren

# One call: annotate -> superimpose -> contacts -> per-interface energies + total
res = run_pipeline("complex.pdb")              # res.scores, res.markup, res.contacts, res.oriented
res = run_pipeline("complex.pdb", reference_aa="A")  # + delta_* : the poly-alanine ΔΦ per interface

# Oracle facade: one structure -> a bundle of ready-to-tabulate frames for the paper
# notebooks (scores, percentile rank, ΔΔG alanine scan, markup, contacts). Configurable
# per-interface potentials and TCR-region selection are forwarded to every milestone.
from tcren import summarize_structure
bundle = summarize_structure("complex.pdb", alanine=True)   # bundle["scores"], ["rank"], ["ddg"], …

# …or the individual steps:
s = parse_structure("complex.pdb.gz")          # also .cif/.cif.gz; import_structure trims the C-gene
classify_chains(s, organism="human")           # TRA/TRB via arda, peptide, MHC
cm = ContactMap.from_structure(s)              # 5 Å contacts + interface partitioning
ranked = score_peptides(cm, ["KQWLVWLFL", "RLLHPHHPL"], tcren())

# Opt-in intra-peptide term: the contacts the peptide makes with ITSELF, which every
# interface energy omits. Off by default (intra_weight=0 leaves every score untouched).
from tcren import intra_peptide_energy
from tcren.potential import mj
cm = ContactMap.from_structure(s, peptide_internal=True)
intra_peptide_energy(cm, mj())                                  # the native peptide's own energy
intra_peptide_energy(cm, mj(), peptide="KQWLVWLFL")             # a candidate on the same pose
score_peptides(cm, cands, tcren(), intra_weight=0.5, intra_potential=mj())
res = run_pipeline("complex.pdb", intra_weight=0.5)             # + scores["peptide_internal"]
```

### Beyond the contact sum

Five things a TCR:pMHC interface does are invisible to a sum over a contact list, and each has its own
instrument here: the one-body / pair split of a potential (`tcren.potential`), peptide backbone
stability under Monte Carlo (`tcren.mechanics.dynamics`), discrete side-chain repacking
(`tcren.energetics.rotamers`), footprint shape (`tcren.topology.footprint`), the pMHC surface a TCR
meets before it binds (`tcren.topology.surface`) and ring-stacking geometry (`tcren.stacking`). What
each measures, and what it was measured on, is in
[**Beyond the contact sum**](https://docs.isalgo.dev/tcren/beyond-the-contact-sum.html).

### CPL response matrices from one template structure

A positional-scanning combinatorial peptide library fixes position *i* to residue *a* and leaves
every other position an **equimolar 1/20 mixture**, so a measured cell is an ensemble mean,
`R[i,a] = E[response | x_i = a]`. `tcren.cpl` predicts that matrix from a single template complex —
each of the twenty residues threaded through the template's own contact map, nothing re-docked,
nothing fitted to any assay.

```python
from tcren import (ContactMap, parse_structure, response_matrix,
                   mutation_effect, position_scan, equimolar_effect)
from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc

s = parse_structure("3HG1.pdb", pdb_id="3HG1")
classify_chains(s, organism="human")
annotate_mhc(s)                       # REQUIRED: without it peptide:MHC is empty and anchors zero out
rm = response_matrix(ContactMap.from_structure(s, cutoff=5.0))

rm.to_frame()                         # the whole matrix, one row per (position, amino acid) cell
position_scan(rm, 5)                  # every substitution at position 5, best first
mutation_effect(rm, 5, "W")           # one cell
equimolar_effect(rm, 5)               # cost of giving position 5 up to the 1/20 mixture
```

**Every cell sums both peptide-bearing interfaces** — TCRen over TCR:peptide plus Miyazawa–Jernigan
over peptide:MHC — because the assay reads *activation*, which needs the peptide presented as well as
the receptor engaged. A position the receptor never touches is an anchor; its TCR term is constant
along the row, so the sum degrades to presentation alone rather than to a special case.

**Two reference states, and a cell is meaningless except against one of them.** A raw Φ carries a
large per-position offset that says only how many contacts the position makes:

| `reference` | cell value | use it for |
|---|---|---|
| `"equimolar"` (default) | `mean_b Φ(x_{i→b}) − Φ(x_{i→a})` | comparing against a **measured** CPL matrix — the mixture is the assay's own background |
| `"wild_type"` | `Φ(x_{i→wt}) − Φ(x_{i→a})` | **mutation scan** / neoantigen ranking off the residue the template carries |

They differ by a per-position constant — how far the template's residue sits above its column mean.
Positive is favourable on both, since lower energy is the better binder. Under `"wild_type"` the
template's own cell is identically zero; under `"equimolar"` it is an ordinary measurement.

### Batch inputs, gzip, archives

```python
from tcren.structure import iter_structures
for pdb_id, structure in iter_structures("batch.tar.gz"):   # file | directory | .tar.gz
    classify_chains(structure, organism="human")
    ...
```

### Canonical orientation, contacts, docking geometry

```python
from tcren.mhc import annotate_mhc
from tcren.docking import canonicalize_structure, superimpose, docking_angles
from tcren.contacts import multi_contacts, ContactDefinition

annotate_mhc(s)
oriented, info = canonicalize_structure(s)     # frame: z=MHC→TCR, y=peptide, x=thin; chains A–E
oriented, info = superimpose(s)                # orient onto data/Canonical2026 by MHC (class+species ensemble)
layers = multi_contacts(s, ContactDefinition(d1=5, d2=8, d3=12))   # heavy-atom / Cβ / Cα
d = docking_angles(s)                          # crossing (~20–70° αβ) + incident angle
```

### 2D complementarity maps & region-pair contacts

```python
from tcren.project2d import (project_structure, residue_markup_table, contacts_table,
                             region_pair_summary)
from tcren.viz import render_complementarity_map, view_pocket_cdr

proj = project_structure(s)                                   # canonical groove plane
svg  = render_complementarity_map(residue_markup_table(s, proj),
                                  contacts=contacts_table(s, threshold=5.0))
region_pair_summary(s, kind="closest")        # contacts per region pair + bond types (cb/ca too)
view_pocket_cdr(s).show()                      # interactive 3D pocket + CDR overlay (py3Dmol)
```

### Publication figures

`tcren.viz.pymol` drives a headless PyMOL to ray-trace figure panels of oriented complexes. Three
scenes cover the usual views, and every panel carries a **labelled axis gizmo** in its corner:

Figures need the `viz` extra (`pip install "tcren[viz]"`) for Pillow, plus a `pymol` binary on
PATH — PyMOL is a separate install, not a Python dependency.

```python
from tcren.viz.pymol import render, overlay_scene, groove_scene, interface_scene
render(groove_scene("1ao7", "data/Canonical2026"), "groove.png")            # peptide in the cleft
render(groove_scene("1ao7", "data/Canonical2026", surface=True), "s.png")   # + molecular surface
render(overlay_scene(ids, "data/Canonical2026"), "overlay.png")             # ensemble, side-on
render(interface_scene("1ao7", "data/Canonical2026", cdr), "iface.png")     # peptide + CDR loops
```

A canonically-oriented structure is only interpretable if the reader can tell which way the frame
points, and `x/y/z` does not tell them — so the arrows are named for what they mean:

| axis | label | direction |
|---|---|---|
| x | `width` | groove width, across the cleft (α1↔α2) |
| y | `N→C` | groove axis, toward the peptide C-terminus |
| z | `TCR` | docking normal, MHC floor → TCR |

The triad is thin, arrow-headed, and turns with the camera. An axis pointing at the viewer
foreshortens to a dot and its label drops to the lower left of it, the usual convention for an axis
normal to the page. These are the three directions the docking-geometry literature uses (SwiftTCR,
TCR3d); only the principal-component ranking differs, because `tcren.docking.frame` fits the whole
complex where those fit the MHC groove alone.

**Colour by which residues carry the score.** Φ is a sum over residue–residue contacts, so it
decomposes exactly: a residue's share is the sum of `φ(a_i, a_j)` over the contacts it makes. The
total says how large the score is; this says what it is made of.

```python
from tcren.viz.pymol import residue_importance, importance_scene
imp = residue_importance(structure)                 # phi + n_contacts, per residue
render(importance_scene("1ao7", CANON, imp), "importance.png")                     # energy share
render(importance_scene("1ao7", CANON, imp, by="n_contacts",
                        spectrum="white_red"), "contacts.png")                     # geometric share
```

CDR3 and peptide residues become sticks on a ramp, everything else stays pale. Blue is favourable
and red unfavourable — the ramp is centred on zero rather than fitted to the range, so those words
keep their meaning even when every contact in an interface is stabilising. Each contact is
attributed to *both* residues it joins, so the per-residue values sum to twice Φ: an attribution,
not a partition.

`render()` is deliberately not a `tcren` subcommand: a figure is a handful of styling choices that
want editing, not a fixed flag set. Pass any PyMOL script body as the scene.

**Explore it interactively** with the [marimo](https://marimo.io) app — pick a structure and scene,
swing the camera and watch the gizmo follow, restyle it, colour by importance with the numbers
beside the render, and rotate a live 3Dmol.js view with the mouse:

```bash
pip install "tcren[marimo]"
marimo run notebooks/pymol_interactive.py       # or `marimo edit` to change the code
```

Worked examples of every view, with images: **[Figure gallery](https://docs.isalgo.dev/tcren/gallery.html)**.

## Modules

| module | what it does |
|---|---|
| `tcren.structure` | parse/write `.pdb`/`.cif`(`.gz`)/`.tar.gz`; the `Atom`/`Residue`/`Chain`/`Structure` model; `iter_structures` |
| `tcren.annotation` | chain typing — TCR loci/CDRs via `arda`, peptide, MHC; αβ/γδ C-gene call |
| `tcren.mhc` | map MHC chains to allele/class/role; partition the groove (helices/floor); NetMHCpan pseudosequence |
| `tcren.contacts` / `contactmap` | closest-atom 5 Å contacts, Cα distances, multi-layer (5/8/12 Å) contact tables, interface partitioning |
| `tcren.potential` | `Potential` (TCRen/MJ/Keskin/MJ1996 + MJ partition energies); `decompose` / `hydrophobicity_fit` — the one-body vs pair split; `derive_tcren` (classic/AM/LOO) with non-redundancy filtering |
| `tcren.stacking` | ring-stacking **geometry** (centroid distance, interplanar angle, vertical/lateral offset) — the directional signal a contact potential cannot see |
| `tcren.energetics` | the interface energy sum (`scoring`), ΔΔG on mutation (`mutation`), rotamer-averaged contacts (`rotamers`); `scoring_rank` percentile-ranks a peptide against a background |
| `tcren.cpl` | CPL response-matrix prediction from one template complex; equimolar and wild-type references; per-position and per-cell queries |
| `tcren.binder` | the pre-energy check that an interface is a plausible dock at all — a rule over contact count and docking geometry, not a model |
| `tcren.recognition` / `descriptors` | the descriptor catalogue: 164 columns in six families, what each means, its units and its known defects (`DESCRIPTORS`, `STATUS`), plus the 40-column interface block this layer computes itself |
| `tcren.score` | the score set — one frozen object, five read-outs (`peptide_score`, `pose_score`, `confidence_residual`, `binder_score`, `channel_scores`), each defined for a single structure |
| `tcren.cohort` / `reliability` | the fit-free predecessor tier: `Q`, `T`, `S`, the AlphaFold band table and the screening cut |
| `tcren.docking` | canonical frame, `superimpose` onto the canonical DB, docking angles, reverse-dock detection |
| `tcren.refine` | peptide substitution + refinement (DOPE MC; CCD/OpenMM/ProMod3/FlexPepDock engines); register QC |
| `tcren.clashes` / `mechanics` | steric-clash report; interface spring-network stiffness + rupture model; peptide backbone dynamics |
| `tcren.topology` | footprint **shape**: coverage entropy / Hill numbers over the CDR-loop × target partition, canonical germline-MHC vs CDR3-peptide preference, α/β contact imbalance, and the footprint's topology (patches, holes, H₀ persistence); the pMHC surface as a height field and the gap between the two faces — no potential, no reference, orientation-free |
| `tcren.project2d` / `viz` | project the interface onto the groove plane; SVG complementarity maps + 3D pocket/CDR views |
| `tcren.pipeline` / `oracle` | one-call structure scoring (`run_pipeline` → Φ, ΔΦ per interface; `summarize_structure`) |
| `tcren.paper` | Nat Comput Sci 2022 reproduction (HF bootstrap, batch annotation, legacy comparison) |

## Data

Structures live in the Hugging Face dataset
[`isalgo/tcren_structures`](https://huggingface.co/datasets/isalgo/tcren_structures), all gzipped:

| folder | contents |
|---|---|
| `Native2022` | the 2022 paper set (oracle) |
| `Native2026` | the comprehensive 2026 TCR:pMHC set the current potential is derived from |
| `Canonical2026` | `Native2026` re-oriented into the canonical frame (`tcren orient`) |

`tcren` reads `.pdb`/`.cif`/`.pdb.gz`/`.cif.gz` and `.tar.gz` batches; an installed library lazily
fetches the canonical reference structures from the Hub when orienting a new complex.

**Where it all lives: `tcren.paths.tcren_home()`.** That one root holds the MHC allele reference
(`database/mhc/`, written by `tcren build-mhc-ref`), its mmseqs index (`data/mhc_cache/`) and the
structure sets (`data/`). It resolves to `$TCREN_HOME` when set; otherwise to the source checkout,
recognised by its `pyproject.toml`, so a development install uses the repo's own `data/`;
otherwise to `$XDG_CACHE_HOME/tcren` (in practice `~/.cache/tcren`), which an installed wheel can
write and which survives an upgrade. `$TCREN_DATA_DIR` overrides the `data/` subdirectory alone.

That `data/` holds `Native2026` (+ `Canonical2026`, gitignored, fetched on demand), `PDB_date.tsv`,
and **`TCRen_potential.csv`** — the 2022 (`karnaukhov2022`) matrix, kept for reproducing
published results; the current default is the bundled `tcren2` (use `-p karnaukhov2022`
for the old one). `Canonical2026`'s
`orient_metadata.json` ships inside the package (`src/tcren/data/`), because the fetch brings down
structures only and an installed library has no repo `data/`.

## Notebooks

Runnable examples under [`notebooks/`](notebooks/) (rendered in the
[docs](https://docs.isalgo.dev/tcren/)):

- `complementarity_map_2d` — 2D interface maps, multiple structural + map views of 1ao7
- `contact_thresholds_and_bondtypes` — region-pair contact counts (closest/Cβ/Cα) + bond types
- `canonical_frame_figures` — canonical-frame QC across the Native2026 set
- `pymol_canonical_figures` — ray-traced PyMOL panels (overlay, groove, interface) by class/species
- `mhc_pseudosequence_mps` — NetMHCpan MHC pseudosequence (MPS) residues vs. peptide contacts
- `example_gil_a02_rs_motif` — GILGFVFTL/HLA-A*02 and the public CDR3β Arg–Ser motif
- `pocket_cdr_3d` — 3D peptide-binding pocket with the CDR loops overlaid (py3Dmol)
- `tcren_analysis` — potential heatmaps (TCRen / MJ / Keskin) and contact distributions
- `natcompsci2022/` — full reproduction of the Nat Comput Sci 2022 analyses

Two of them run the score set end to end, from fetching structures to ranking them:

- `score_vdjdb_panel` — receptor ranking for a fixed epitope on the balanced 22-cohort VDJdb panel
  (1,089 TCRmodel2 models), reported one cohort at a time with template-covered and template-free
  apart, and the five channel scores per cohort
- `rank_peptides_cpl` — peptide ranking for a fixed receptor on the combinatorial-peptide-library set
  (7 clones, 2,103 models): per-clone ROC, the graded activation read-out, and a whole response
  matrix predicted from one template

Four **marimo** apps ship alongside them (`pip install 'tcren[marimo]'`, then `marimo run <file>`):

- `surface_topology.py` — elevation / charge / hydropathy maps over the groove, and the
  featureless-vs-bulged epitope comparison against the structures the literature names
- `pymol_interactive.py` — a PyMOL render explorer over the canonical scenes (overlay, groove,
  interface, residue importance)
- `confident_negatives.py` — a generator's confidence read together with the coordinates
- `potts_contact_map.py` — the predicted contact-frequency map beside the contacts a structure made

## Performance

Per-stage wall time (best of *n*) on a TCR-pMHC complex (1ao7), Apple M-series, single thread
(`RUN_BENCHMARK=1 pytest -k benchmark -s` to reproduce the core stages):

| stage | time | notes |
|---|---|---|
| parse a gzipped structure | ~17 ms | `.pdb.gz` / `.cif.gz` |
| contact map (5 Å, cKDTree) | ~9 ms | per structure |
| score 1000 candidate peptides | ~11 ms | ~10 µs/peptide (vectorised) |
| ΔΔG alanine scan (9-mer) | ~11 ms | virtual-matrix; no atoms move |
| the score set, all six read-outs | ~5 ms for one structure, ~34 µs/structure over 1,089 | frozen transform + two Gaussians; no structure re-read |
| peptide refine (2000-step DOPE MC) | ~320 ms | knowledge-based rigid-body refinement |
| annotate (MHC map, 1 structure) | ~670 ms | one mmseqs2 search |
| **annotate (TCR + MHC), batched** | **~0.2 s/structure** | one mmseqs2 call for the whole set; vs ~1.5 s/structure unbatched |
| superimpose onto the canonical DB (per query) | ~2.8 s | aligns to every same-class DB structure |

| peak RSS | value | notes |
|---|---|---|
| single-structure pipeline (no orient) | ~200 MB | parse → annotate → contacts → score → refine |
| + `superimpose` (loads canonical DB) | ~780 MB | holds Canonical2026 in RAM; skip with `--no-superimpose` |

Annotation is the only network/compute-heavy step and is always **batched** (one mmseqs2 search over
all chains; mmseqs2 parallelises internally — never per-structure, never Python-threaded). Threads are
used only for the embarrassingly-parallel, mmseqs-free stages (structural alignment, write, rendering):
`tcren orient -t N`. Screening a peptide/TCR panel is embarrassingly parallel — references are
annotated and oriented **once**, so the hot loop is just refine + contacts + score per complex.

## Tests

```bash
pytest -m "not slow"          # unit + fast regression (the CI gate)
pytest                        # add the arda/mmseqs-backed regression tests
RUN_BENCHMARK=1 pytest -k benchmark -s
```

## Project state

- [CHANGELOG.md](CHANGELOG.md) — what has landed, per release, with the measurement for each.
- [STATUS.md](STATUS.md) — where the modules stand, and the known caveats.
- [ROADMAP.md](ROADMAP.md) — where it is going, and what each direction is waiting on.
- [BENCHMARKS.md](BENCHMARKS.md) — achieved accuracy.

## Citing

**TCRen** is free for academic and non-commercial use. If you use it, please cite our latest 
[Nature Computational Science 2024 paper](https://www.nature.com/articles/s43588-024-00653-0):

```
Karnaukhov VK, Shcherbinin DS, Chugunov AO, Chudakov DM, Efremov RG, Zvyagin IV, Shugay M. Structure-based prediction of T cell receptor recognition of unseen epitopes using TCRen. Nat Comput Sci. 2024 Jul;4(7):510-521. doi: 10.1038/s43588-024-00653-0. Epub 2024 Jul 10. PMID: 38987378.
```
