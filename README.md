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
| Percentile-rank a peptide vs background | `tcren rank` | `percentile_rank` |
| ΔΔG of mutations (alanine scan / neoantigen) | `tcren ddg` | `alanine_scan`, `neoantigen_ddg` |
| **Predict a CPL response matrix from a template** | `tcren cpl` | `response_matrix`, `mutation_effect`, `position_scan`, `equimolar_effect` |
| Binder vs non-binder for a TCR model | `tcren features` + `tcren recognize --features` | `cohort.p_native`, `cohort.q_score` |
| **Every interface descriptor, in five families (four by default)** | `tcren features` | `recognition_table(include=...)`, `descriptors` |
| **All interface descriptors + joint P(real)** | `tcren recognize` | `recognition_features`, `real_probability` |
| **P(native)** — the channels combined by a latent-class Bayes network | `tcren recognize --features` | `cohort.p_native` |
| **Is *this* model worth believing?** — `S_free`, a calibrated `p_binder`, and the generator diagnostic | `tcren assess` | `reliability.s_free`, `p_binder`, `af_band` |
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
`tcren._fold` (CCD loop closure) and `tcren._geom` (interface geometry for `tcren binder`). TCR
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
tcren scoring -s complex.pdb -o scores.csv --intra-weight 0.5   # reports F_pep_int separately

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

# Rank candidate receptors against a fixed pMHC. S_free is the score to ship -- fit-free, and
# defined for one structure. P_native comes out of the same call: a latent class over geometry,
# footprint topology and contact energetics, fitted by EM with no binding label (TCRvdb macro
# ROC 0.832 / PR 0.849, against AlphaFold ipTM 0.795 / 0.783), but cohort-relative, so score the
# whole candidate set together rather than one structure at a time.
tcren features  -s candidates/ -o feats.tsv
tcren recognize --features feats.tsv -o scores.tsv

# One TSV per structure: every interface descriptor (geometry + energies) + joint P(real).
tcren recognize -s my_pdbs/ -o recognize.tsv          # descriptors + p_real + p_real_bn, one row/PDB

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

## One table per structure: descriptors, energies & the joint recognizer

Two commands, two jobs. **`tcren features` reads structures and writes descriptors**;
**`tcren recognize` turns descriptors into scores.** The feature pass is the expensive half, so it
runs once and the scoring pass can be repeated for nothing.

```bash
tcren features  -s my_pdbs/ -o feats.tsv                   # the four default families (--all adds kinetics)
tcren features  -s my_pdbs/ -o shape.tsv -i topology       # one family -- and only it is computed
tcren recognize --features feats.tsv -o scores.tsv         # Q, P_native + the G/T/E channels
```

Descriptors are catalogued in five **families**, four of them computed by default (`kinetics` is
opt-in), split by what each is invariant under — which is also the axis along which they carry
independent evidence:

| family | what it is | invariance |
|---|---|---|
| `placement` | where the receptor sits in the groove frame — angles, TCRdock parameters, ride height, shift, offset, the CDR3 loop frames | frame-**dependent** |
| `interface` | how much contact and of what chemical kind — buried area, contact counts and types, hydrogen bonds, clashes | SE(3)-invariant |
| `topology` | the **shape** of the contact set, free of its size — coverage entropy, Hill numbers, Betti numbers, persistence entropy, canonical preference | SE(3)-invariant |
| `energetics` | statistical-potential interface energies `F` and their poly-alanine references `dF` | SE(3)-invariant |
| `kinetics` | the interface as a spring network — stiffness, rupture, coupling residues (off unless asked) | — |

`tcren recognize -s my_pdbs/` reads the structures itself, skipping the feature file:

```bash
tcren recognize -s my_pdbs/ -o recognize.tsv               # descriptors + p_real
tcren recognize -s my_pdbs/ -o scored.tsv --mechanics      # + the spring-network kinetics terms
```

| what you want | columns in `recognize.tsv` |
|---|---|
| **(a) energy** — `F` per interface (TCRen on TCR:peptide, MJ on presentation) + poly-alanine `dF` + loop parts | `F_tcr_pep`, `F_tcr_mhc`, `F_pep_mhc`, `dF_tcr_pep`, `dF_pep_mhc`, `F_cdr12`, `F_cdr3a`, `F_cdr3b` |
| **(a′) intra-peptide** (`--full`) — the peptide's contacts with *itself*, which every interface sum omits | `F_pep_int`, `n_pep_int` |
| **(b) geometry** — every docking + interface descriptor | `pitch`, `crossing`, `crossing_signed`, `dock_d`, `dock_torsion`, `dock_{tcr,mhc}_u{y,z}`, `extent`, `chain_balance`, `burial`, `n_contacts_{tp,tm}`, `n_pep_contacted`, `ct_{tp,tm}_*` |
| **(c) cohort scores** — no training set, no binding label; written by `tcren recognize --features`, not by `-s` | `Q` — interface quality; `G` / `T` / `E` — the geometry, topology and energetics channels on their own; `P_native` — the three combined; `S_free` and its calibrated `p_binder`. See [`tcren.cohort`](src/tcren/cohort.py), [`tcren.reliability`](src/tcren/reliability.py) |
| **(d) joint P(real)** ~ Bayesian model over energy + geometry | `p_real` — distribution-aware Bayesian **logistic** (5-fold CV AUC 0.885); `p_real_bn` — the Gaussian **BN** variant |

### Is *this* model worth believing? — `tcren assess`

A co-folding model will seat any TCR against any peptide. `assess` answers the three questions a
caller actually has about one structure, from coordinates alone:

```bash
tcren features -s models/ -i placement,interface,topology,energetics -o feats.tsv
tcren potts   score -s models/ -o potts.tsv            # writes neg_energy, the energy block
# join potts.tsv's neg_energy onto feats.tsv on the structure id, then:
tcren assess --features joined.tsv -o assessed.tsv
```

- **Reliability** — `S_free` = `Q/sd_Q + T/sd_T + (Pi - mu)/sd_Pi`, three fit-free directional
  blocks each divided by its own native spread, so they carry equal weight in native-sd units.
  Nothing is fitted at score time, so **it is defined for a single structure**. `p_binder` turns it
  into a probability through a frozen out-of-fold Platt link.
- **Ranking** — the structure's rank and percentile inside the set, and the expected precision at a
  recall budget.
- **The generator diagnostic** — which AlphaFold confidence band the model falls in, how often
  models in that band turned out to be non-binders, and what `S_free` still separates *inside* it.
  On a balanced 22-cohort VDJdb panel the top ipTM decile is 26.2% [18.7, 35.5] non-binders, and is
  also the band where `S_free` reads highest.

Without the joined `neg_energy` column `assess` emits the two-block `Q + T` form and says so in its
report rather than imputing the missing block. `P_native` is still emitted by `recognize`, now
documented as cohort-refit and not the recommended score: it refits per call, raises when a cohort
has fewer rows than features, and its value depends on what else was scored alongside it.

**Where the joint model lives.** `p_real` is the frozen recognizer we derive from real crystals vs
wrong-TCR *shuffled* decoys: code in [`tcren.recognition`](src/tcren/recognition.py)
(`recognition_features` → `real_probability`), coefficients shipped in
`src/tcren/data/shuffle_logistic.json.gz`, and the full derivation (PyMC fit, encoding, ROC/PR,
posterior forest) in the technical appendix, which lives with the manuscript rather than here —
`logistic_stan/`, with the Gaussian-BN companion in `shuffle_bn/`. Decoys come from
`tcren shuffle`.

**(c) physics of the interaction.** The koff proxies fold into the same table with `--mechanics`;
only the mutation scan, which is per-residue rather than per-structure, needs its own command:

```bash
tcren recognize -s models/ --mechanics -t 0 -o out.tsv    # every per-structure descriptor, one table
tcren ddg       -s complex.pdb -o ddg.csv     # per-residue alanine / neoantigen ΔΔF (fast virtual matrix)
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
from tcren.recognition import recognition_features, real_probability
feats = recognition_features("complex.pdb")    # dict of the 34 descriptors (RECOGNITION_FEATURES)
p = real_probability(feats)                     # {"logistic": P(real), "bn": P(real)}
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

### What a contact potential can and cannot express

A contact energy is not purely an interaction: burying a residue against *any* partner costs
something that depends on that residue alone. `decompose()` separates the two exactly, and only the
pair part `J` is beyond what a per-position model can already write down.

```python
from tcren.potential import mj, mj1996, mj_partition_energy

d = mj1996().decompose()          # e(a,b) = mean + H(a) + H(b) + J(a,b), J double-centred
d.h("F"), d.j("F", "W")           # one-body term; the genuinely pairwise remainder
d.energy("F", "W")                # reassembles the original value

f = mj1996().hydrophobicity_fit()  # C0 + C1(q_a + q_b) + C2 q_a q_b
f.r2, f.eigenvalue_share           # 0.98 on MJ1996, 0.85 on the bundled mj
mj_partition_energy()["F"]         # 4.37 — MJ's own one-body scale (larger = more hydrophobic)
```

Where a potential has that shape the interaction term is only `C2·q_a·q_b`, so it **cannot prefer
one pair of side chains over another of equal hydrophobicity**. Both calls refuse a directed
potential — TCRen is TCR→peptide and must not be split this way.

### Peptide conformational stability: what a contact model cannot see

A contact potential scores whichever conformation it is handed. It cannot tell a peptide that its
own side chains **hold** in the TCR-facing conformation from one that merely happens to have been
modelled there — both present the same contact list. `tcren.dynamics` puts the backbone in motion:
it samples peptide φ/ψ by Metropolis Monte Carlo against DOPE and reports how far the peptide
wanders, not a better pose.

```python
from tcren import peptide_stability, stability_table
peptide_stability(structure).rmsf                       # ensemble spread, A -- larger = floppier
stability_table([s1, s2]) ["delta_rmsf"]                # intra-peptide term ON vs OFF, paired
```

**The hypothesis it was built to test** (Sewell, 2026-08): intra-peptide interactions stabilise the
productive bulge a TCR reads, and *"poor binders could perhaps still make many contacts but fail to
stabilise the productive peptide conformation"* — which would explain why an additive contact model
describes some systems well and others badly.

Tested on the CPL set: ~160 best-binder and ~160 worst-binder modelled complexes for each of seven
clones, 2102 structures (`scripts/sewell_stability.py`). AUC is best-vs-worst discrimination.

| clone | contact energy | **stability** |
|---|---|---|
| ila1 | 0.348 | **0.862** |
| 868 | 0.537 | **0.677** |
| sb27 | 0.570 | **0.934** |
| mel8 | 0.690 | **0.876** |
| 4c6 | **0.955** | 0.519 |
| 1e6 | **0.973** | 0.707 |
| mel5 | **0.974** | 0.859 |

**Stability beats the contact energy in 4/4 clones where the contact model fails, and 0/3 where it
works.** Mean AUC over the failing clones goes 0.536 → 0.837; over the working ones the contact
energy stays ahead (0.967 vs 0.695). Combining the two (within-clone z-sum) lifts the mean AUC from
0.721 to 0.826, improved in 5/7 — though with seven clones that paired test is underpowered
(Wilcoxon p = 0.22).

**The intra-peptide term is a switch, and flipping it does what the hypothesis says.** Removing the
peptide's contacts with itself lets the *best* binders' backbones wander further (Δrmsf = +0.021 Å,
SE 0.005, i.e. 4.4σ) and leaves the *worst* binders unchanged (+0.002 Å, SE 0.007); best vs worst
p = 0.042. The same term also sharpens the stability discrimination itself, by +0.024 AUC on average
and in 5/7 clones (Wilcoxon p = 0.078).

So the **mechanism** Sewell proposed is supported, while the **system** he guessed is not: 4c6 is one
of the clones the contact model handles well here (0.955), and the ones it fails on are ila1, 868,
sb27 and mel8. Caveats worth carrying: these are modelled structures, the MC is knowledge-based
rather than MD (no solvent, no force field, no time), Δrmsf is a mechanistic signal and not a useful
classifier on its own (AUC 0.526), and every clone-level test has n = 7.

### Side-chain repack: what a local minimiser cannot do

`tcren.repack` (native `_relax.repack`) places every side chain in the χ rotamer the DOPE potential
prefers. The rigid-body refiner moves the peptide and leaves every χ where it found it, so a
full-atom model whose side chains a predictor placed keeps them — which is most of why a pairwise
contact energy stops discriminating on AlphaFold poses.

```python
from tcren import repack
fixed, report = repack(structure)        # report: n_conformers, energy, p_best per residue
```

Like-for-like — same wrong-rotamer input (χ1 rotated 120°), same 33–42 side-chain atoms, same
crystal reference:

| | peptide side-chain RMSD (Å) | time |
|---|---|---|
| input (wrong χ1) | 4.131 | — |
| **`repack`** | **2.364** | **6 ms** |
| OpenMM (anchor-restrained minimisation) | 4.133 | 3103 ms |

OpenMM leaves them where they are. That is not a defect in OpenMM: a local minimiser cannot cross
the torsional barrier between two rotamer basins, so relaxing clashes and re-sampling rotamers are
different operations and only a discrete packer does the second. Over eight structures the packer
recovers side-chain RMSD 3.93 → 1.66 Å, 8/8 improved, median 6 ms.

It rotates the side chains a model **has** — it cannot rebuild ones `substitute_peptide` stripped;
that is side-chain *construction*, still open (`refine/CPP_REWRITE.md`).

### Footprint shape: what the contacts say before they are scored

Every other scorer here sums over contacts. `tcren.footprint` reads the same contact map as a
**shape** — which of the six CDR loops touched what, and whether the resulting footprint is one
connected patch. No potential, no reference structure, no fitted parameter, and **no canonical
orientation**: every feature is invariant under rigid motion, so unaligned inputs are fine. Only
chain typing and CDR markup are needed, which the CLI does in one batched annotation pass.

Coverage is the composition over cells — the 6 CDR loops × {peptide, MHC}, optionally splitting the
peptide into thirds — summarised by the normalised Shannon entropy and by the Hill numbers
([Hill 1973](https://doi.org/10.2307/1934352), [Jost 2006](https://doi.org/10.1111/j.2006.0030-1299.14714.x)),
where `D2` is the *effective number of engaged cells*. Topology joins the contacted pMHC residues at
a Cα threshold and builds the flag complex: `fp_b0_*` counts footprint patches and `fp_b1_*` its
holes. Coverage and topology are only weakly related, which is why they belong in one channel and
why that channel is read as `T` — the shape posterior of `p_native`, not a hand-written z-sum.

```bash
tcren features -s structures/ -i topology -o shape.tsv
```

```python
from tcren.cohort import p_native
from tcren.footprint import footprint_batch, footprint_features

row = footprint_features(structure)          # one dict, 29 features at the default two radii
row["D2_pep24"], row["fp_b0_r7"], row["L_canon"]

table = footprint_batch("structures/")       # polars frame, one row per structure
T = p_native(table, channels=("topology",))  # the shape score, cohort-relative
```

Note the cyclomatic number of the bipartite contact graph (`E − V + C`) is deliberately not offered:
with of order thirty contacts among of order thirty residues it is dominated by `E` and simply
tracks interface size. The patch count is scale-free instead.

### Surface topology: what a TCR meets before it binds

A contact potential scores an interface that already exists. `tcren.surface` describes the pMHC
*beforehand*: the peptide sits in a groove between two helices, and a TCR coming down meets one
surface, so the descriptor is a height field `h(x, y)` over that groove with hydropathy and charge
painted on. Method follows [SURFMAP](https://doi.org/10.1021/acs.jcim.1c01269) (surface shell,
per-cell feature, 8-neighbour smoothing, Manhattan map distance, hierarchical tree) and
[Protein Surface Topography](https://doi.org/10.1074/jbc.RA119.010494) (centre the chart on the
functional site). A flat raster rather than SURFMAP's equal-area spherical chart, because the
TCR-facing surface is an open, near-planar patch that a plane does not distort.

```python
from tcren import surface_map, surface_stats, surface_distance, surface_tree
smap = surface_map(structure)              # channels: h, phobic, charge; source: peptide/helix/floor
surface_stats(smap)["frac_above_ridge"]    # how much peptide surface clears the MHC helix crests
ids, d = surface_distance([m1, m2, m3])    # pairwise map distance -> epitopes cluster
```

Two things worth knowing, because both were defects first:

* **The frame is refit from every structure** — z from the groove-floor plane normal, **y from the
  peptide**, origin on the peptide centroid. The floor's own principal axis is *not* the groove axis
  (its β-strands run across the groove), which put the two helices diagonally across the map. Because
  the frame is intrinsic, maps compare without prealigning the inputs — SURFMAP's standing caveat.
* **Heights come from ray casting in the groove frame**, not from Shrake-Rupley surface points.
  Sphere sampling is fixed in global axes, so the same structure rotated gave a different map (median
  cell moved 1.35 Å, `relief` by 19%). Ray casting is exactly equivariant and needs no probe test —
  the highest surface in a column is by definition the one nothing is above.

**"Featureless" becomes a number.** Over the 374 Canonical2026 complexes (230 distinct epitopes), the
epitopes the literature *names* as featureless and as bulged separate completely:

| epitope | source | rank by `frac_above_ridge` | `frac_above_ridge` | `relief` (Å) |
|---|---|---|---|---|
| LPEPLPQGQLTAY | EBV BZLF1 13-mer, HLA-B\*35 — **bulged** | **2 / 230** | 0.749 | 2.81 |
| HPVGEADYFEY | HCMV pp65 11-mer, HLA-B\*35:08 — **bulged** | 5 / 230 | 0.562 | 3.59 |
| EPLPQGQLTAY | EBV BZLF1 11-mer, HLA-B\*35 — **bulged** | 8 / 230 | 0.416 | 2.54 |
| LLFGYPVYV | HTLV-1 Tax, HLA-A\*02:01 — prominent P5-Tyr | 46 / 230 | 0.145 | 2.15 |
| GILGFVFTL | influenza M1, HLA-A\*02:01 — **featureless** | 139 / 230 | **0.000** | 1.16 |
| TAFTIPSI | HIV RT 8-mer, HLA-B\*51:01 — **featureless** | 205 / 230 | **0.000** | 0.95 |

Five of the eight most-protruding epitopes are literature-named bulged HLA-B\*35 epitopes; both named
featureless ones have *no* peptide surface clearing the helix crest at all. Structure-level AUC is
1.000 on `relief`, `peak_to_valley` and `frac_above_ridge` (p ≤ 0.001, 9 featureless vs 5 bulged
structures) — though with two distinct epitopes per group that is a 2-vs-2 comparison, so the
properly-powered evidence is the trend over all 279 class-I structures: `frac_above_ridge` rises
0.054 (8-mers) → 0.569 (13-mers), Spearman on `relief` +0.414, p = 5.5e-13.

`notebooks/pnative_channels.py` (marimo) runs the released scoring path over a directory of
structures — one featurisation pass, one latent-class fit per channel, then `P_native` — and ends on
the correlation whose sign says whether a pose was copied from a template.
`notebooks/surface_topology.py` (marimo) draws the elevation / charge / hydropathy maps and
reproduces this comparison.

### Ring stacking: the geometry identity cannot carry

A contact potential scores a pair by identity, so two rings face-to-face at 3.5 Å score exactly like
the same two residues brushing past edge-on. This measures the difference and returns **no energy**:

```python
from tcren import ring_stacking
ring_stacking(structure, cutoff=7.5)   # centroid_distance, interplanar_angle, vertical, lateral
```

`interplanar_angle` near 0 is face-to-face, near 90 edge-to-face; a parallel-displaced stack shows a
small `vertical` with a few Å of `lateral`. Proline is included — its pyrrolidine ring packs face-on
against aromatics through CH–π contacts.

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
from tcren.orient import canonicalize_structure, superimpose, docking_angles
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
TCR3d); only the principal-component ranking differs, because `tcren.orient.frame` fits the whole
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
| `tcren.scoring` / `scoring_rank` | substitution scoring of candidate peptides; percentile rank vs a background |
| `tcren.ddg` | fast virtual-matrix ΔΔG — alanine scan, neoantigen mutants |
| `tcren.cpl` | CPL response-matrix prediction from one template complex; equimolar and wild-type references; per-position and per-cell queries |
| `tcren.binder` | binder/non-binder classifier from AF-orthogonal interface geometry |
| `tcren.recognition` | 34-descriptor extractor (`recognition_features`) + frozen real-vs-shuffled recognizers — distribution-aware Bayesian logistic + Gaussian BN — for joint `P(real)` |
| `tcren.orient` | canonical frame, `superimpose` onto the canonical DB, docking angles, reverse-dock detection |
| `tcren.refine` | peptide substitution + refinement (DOPE MC; CCD/OpenMM/ProMod3/FlexPepDock engines); register QC |
| `tcren.clashes` / `mechanics` | steric-clash report; interface spring-network stiffness + rupture model |
| `tcren.footprint` | footprint **shape**: coverage entropy / Hill numbers over the CDR-loop × target partition, canonical germline-MHC vs CDR3-peptide preference, α/β contact imbalance, and the footprint's topology (patches, holes, H₀ persistence) — no potential, no reference, orientation-free |
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

Two **marimo** apps ship alongside them (`pip install 'tcren[marimo]'`, then `marimo run <file>`):

- `surface_topology.py` — elevation / charge / hydropathy maps over the groove, and the
  featureless-vs-bulged epitope comparison against the structures the literature names
- `pymol_interactive.py` — a PyMOL render explorer over the canonical scenes (overlay, groove,
  interface, residue importance)

## Performance

Per-stage wall time (best of *n*) on a TCR-pMHC complex (1ao7), Apple M-series, single thread
(`RUN_BENCHMARK=1 pytest -k benchmark -s` to reproduce the core stages):

| stage | time | notes |
|---|---|---|
| parse a gzipped structure | ~17 ms | `.pdb.gz` / `.cif.gz` |
| contact map (5 Å, cKDTree) | ~9 ms | per structure |
| score 1000 candidate peptides | ~11 ms | ~10 µs/peptide (vectorised) |
| ΔΔG alanine scan (9-mer) | ~11 ms | virtual-matrix; no atoms move |
| binder P(bind) (features + model) | ~49 ms | native geometry, no external tool |
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
