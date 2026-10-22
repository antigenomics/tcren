# Quantitative contact-prior dossier — 18 October 2026

Purpose: resolve the source plan's literature placeholders into inspectable
quantitative inputs before implementing a new contact model. This is source
research, not a native-cohort rerun or a fitted resource. Metadata returned by
PubMed is preserved in `contact-prior-metadata.json`; no bibliographic fields
were reconstructed from memory. No external software was executed or installed.

## Resolved antibody reference

The source plan's unnamed large antibody survey is identifiable as
[PMID 38161735](https://pubmed.ncbi.nlm.nih.gov/38161735/),
[DOI 10.1016/j.csbj.2023.11.056](https://doi.org/10.1016/j.csbj.2023.11.056).
The rights-eligible PMC full text was retrieved completely through the PubMed
connector (PMCID PMC10755492; 45003 body-text characters). The paper reports
1833 nonredundant complexes and more than 850000 atom contacts, using an inclusive
5 Å heavy-atom cutoff. It identifies the deposited tabular output as
[DTU.22555672](https://doi.org/10.11583/DTU.22555672), with
[producer code](https://github.com/andreasvisbech/Ab_interface_mapping).

The pinned producer commit is `39ac0b3db00f925493b13d6b38d0998bce86bda3`.
`S3_extract_protein_multi.py`, Git blob
`8ecfa29d4b9e6dd4008ff1e07fdf5ea033fe7a3f`, records both contacting atom
coordinates, identities, residue/chain context and distance, and writes
`Output_contact.parquet` (or CSV). This makes atom identities and relative
geometry potentially recoverable from the published output without rerunning
the full antibody corpus. Its endpoint coordinates alone do not provide every
antecedent atom needed for hydrogen-bond angles. Verify the deposited schema,
license, raw identities, count/duplicate conventions and TCR overlap before use.

Retrieval status: the public Figshare metadata API returned HTTP 403; the exact
DTU landing-page route returned HTTP 403 through web retrieval. A separate curl
landing-page request returned success with zero bytes and is not a successful
data retrieval. No deposited contact file or license has been verified. These
responses do not establish dataset absence. The raw data remain a retrieval
dependency, rather than permission to invent frequencies from plotted results.

Subsequent direct in-app browser navigation **did reach the public deposit**.
It lists `Output_contact.parquet` (38.66 MB; file ID 42341739),
`Output_ref.parquet` (62.71 MB; file ID 42341748), their larger CSV counterparts,
and `README.txt` (5.4 kB; file ID 43561332). Thus file availability is confirmed
in the visible deposit, superseding the earlier endpoint-only uncertainty.
The README preview was unavailable; its browser download attempt stalled and
the tool terminated with no verified file. Direct institutional curl downloads
yielded zero bytes. The later successful retrieval and validation below supersede
these access limits.

### Downloaded deposit and integrity checks

The visible deposit specifies **CC BY 4.0** and publishes contact-file MD5
`9bdbfbab026ebe4727748af608dd5c5f`. Ordinary GET requests through Figshare's
standard public download endpoint succeeded on 18 October 2026. HEAD requests
failed at the storage endpoint and are not a valid download-availability test.
The author-supplied downloaded contact file independently matches the same MD5
and SHA256. No author file was modified.

Raw inputs are retained under `data/contact-priors/DTU-22555672/`, excluded from
Git. Exact download commands, run from the repository root:

```sh
curl -fLsS --max-time 60 https://ndownloader.figshare.com/files/42341739 -o data/contact-priors/DTU-22555672/Output_contact.parquet
curl -fLsS --max-time 90 https://ndownloader.figshare.com/files/42341748 -o data/contact-priors/DTU-22555672/Output_ref.parquet
curl -fLsS --max-time 30 https://ndownloader.figshare.com/files/43561332 -o data/contact-priors/DTU-22555672/README.txt
```

| File | Bytes | SHA256 |
|---|---:|---|
| Output_contact.parquet | 40532886 | bc746a214988f16bb376c8d56bf879739b7704b6dafa2eac032893509d7a7d4b |
| Output_ref.parquet | 65760941 | f1b736fa32b05c2b4fc3efbf335ae066bdea7b6823cb740b215f2c206741b4a1 |
| README.txt | 5525 | 3a0af3e6c60b0022fdb1d82f8dbefd414a522ac5b72b94e57722bc79aa410f3a |

Integrity/schema inspection used installed Polars and NumPy with one compute
thread, projecting the contact identity, coordinate, context and distance columns.
It did not regenerate the reference corpus or fit a prior/scorer. Observed counts:

- **856751 deposited contact rows from 1833 distinct PDB IDs**, representing
  856222 distinct physical atom-pair keys and 101974 physical residue-pair keys.
  The physical atom key is `(PDB, Ab chain ID, Ab resi ID, Ab atom,
  Ag chain ID, Ag resi ID, Ag atom)`; the residue key omits atom names.
- All **529 excess rows** repeat keys in **4k94** exactly twice. Repeated keys
  agree in both endpoint coordinates, distance, antigen type and antibody type;
  no conflicting geometry/context was found. Preserve the raw file; count each
  physical key once when constructing the prior. These are not 529 new contacts.
- None of the 15 projected identity/context/coordinate/distance columns contains
  nulls. All 856751 stored distances are finite and within 5 Å; their minimum
  is 1.3316469192504883 Å and maximum is 4.999998092651367 Å. Recomputing distances
  from the stored float32 endpoint coordinates gives maximum absolute difference
  4.76837158203125e-7 Å. A short distance is not by itself an interaction label.
- **5154707 reference atom rows from 1833 PDB IDs**, with no repeated
  `(PDB, Ab chain ID, Ab resi ID, Ab atom)` key. The reference table contains
  antibody atoms only, as the README explains; it is not a complete antigen
  coordinate table. Thus it cannot guarantee both endpoints' Cβ orientations
  or all hydrogen-bond antecedents without additional raw structures.
- The downloaded contact corpus shares **two deposited PDB IDs** with the
  374-file Native2026 raw inventory: **3gjf** (582 contact rows) and **3hae**
  (693 contact rows). PDB-ID overlap does not establish identical interfaces;
  do not claim an independent external corpus or silently count these as
  independent support. Exclude overlapping structures from the external prior
  before native updating, with exact retained counts recorded by its producer.

| Antigen type | Antibody type | Deposited atom-contact rows | Distinct PDB IDs |
|---|---|---:|---:|
| peptide | Fv | 151799 | 367 |
| peptide | VH sdAb | 5377 | 14 |
| peptide | VL sdAb | 377 | 2 |
| protein | Fv | 515651 | 1058 |
| protein | VH sdAb | 182243 | 388 |
| protein | VL sdAb | 1304 | 4 |

These are deposited-row counts before physical-key deduplication and overlap
exclusion. The six PDB counts sum to 1833. No quantitative plausibility prior,
effective prior strength, transformed-angle model or classification improvement
has yet been derived from the download.

## PPI geometry and atomic-property source

[PICCOLO, PMID 21801404](https://pubmed.ncbi.nlm.nih.gov/21801404/),
[DOI 10.1186/1471-2105-12-313](https://doi.org/10.1186/1471-2105-12-313),
Table 1 defines van der Waals contact as distance less than summed atom radii
plus 0.5 Å; a van der Waals clash is below summed radii. Its hydrogen-bond rule
uses donor–acceptor distance <3.9 Å, hydrogen–acceptor distance <2.5 Å and
three specified angles >90 degrees. Those rules require hydrogen placement
and acceptor antecedents; heavy-atom distance alone is insufficient.

The paper links additional file 1 (`atomic_properties.xls`, residue-specific
radii and property flags) and additional file 2 (observed/expected residue
contact matrices). The first supplement request encountered a browser check;
the second failed retrieval. Neither binary file was inspected. The paper's
ionic/aromatic 6 Å criteria exceed the requested 5 Å observation window;
conditioning/truncation must be explicit. These are classification definitions,
not Gaussian means/variances or native-validity probabilities. The available
paper portion supplies a concrete route for PPI prior counts, not a validated
prior table. Source license: CC BY 2.0 according to PMC rights metadata.

## Inspectable Arpeggio criteria

The maintained [PDBe Arpeggio configuration](https://github.com/PDBeurope/arpeggio/blob/258855b8ba13447f2776b232ca32884d637c6a9c/arpeggio/core/config.py)
was retrieved at commit `258855b8ba13447f2776b232ca32884d637c6a9c`, Git blob
`1316e9e1c45a44bbf6d97890c32ee18549d375b3`. Its constants include:

| Configuration quantity | Value | Units |
|---|---:|---|
| Hydrogen-bond distance | 3.9 | Å |
| Hydrogen-bond angle | 90 | degrees |
| Angle-free polar distance | 3.5 | Å |
| Ionic distance | 4.0 | Å |
| Hydrophobic distance | 4.5 | Å |
| Aromatic centroid distance | 6.0 | Å |

These are source constants, not a complete executable rule specification;
comparison operators, angle vertices, atom typing and aromatic subtypes must
be traced through `interactions.py` before implementation. The configuration's
`CONTACT_TYPES_DIST_MAX = 4.5` is not the plan's 5 Å contact definition.
Do not silently replace the latter. The
[paper record](https://pubmed.ncbi.nlm.nih.gov/27964945/)
([DOI](https://doi.org/10.1016/j.jmb.2016.12.004)) distinguishes polar labels
without angular criteria from hydrogen-bond assignments. Its feature labels
are useful chemistry annotations, not a model of statistical rarity.

## Other verified source routes

- [SIRIUS, PMID 2308168](https://pubmed.ncbi.nlm.nih.gov/2308168/),
  [DOI 10.1016/0022-2836(90)90268-Q](https://doi.org/10.1016/0022-2836(90)90268-Q):
  verified metadata/abstract identifies the underlying geometry-analysis method
  and a 400-combination side-chain database. No numeric atlas distribution has
  been retrieved; book-review descriptions are not numerical prior evidence.
- [Antibody interaction survey, PMID 24488795](https://pubmed.ncbi.nlm.nih.gov/24488795/),
  [DOI 10.1002/prot.24527](https://doi.org/10.1002/prot.24527): metadata and
  abstract verified; quantitative tables/criteria remain uninspected.
- [Interface hydrogen-bond survey, PMID 27153634](https://pubmed.ncbi.nlm.nih.gov/27153634/),
  [DOI 10.1093/bioinformatics/btw197](https://doi.org/10.1093/bioinformatics/btw197):
  metadata verified; full-text connector failed twice and web retrieval met a
  browser check. No numerical distribution has been extracted from its abstract.
- [Antibody/single-domain comparison, PMID 37533864](https://pubmed.ncbi.nlm.nih.gov/37533864/),
  [DOI 10.3389/fimmu.2023.1231623](https://doi.org/10.3389/fimmu.2023.1231623):
  complete rights-eligible body retrieved. Its interaction analysis uses a
  4.5 Å cutoff, so its counts cannot be pooled into 5 Å statistics unchanged.

## Consequences for the model design

1. The antibody reference is now identified and its producer schema inspected.
   Its exact deposit and README are now verified. Collapse
   unordered endpoint tuples only after preserving which endpoint is antibody
   or antigen; do not mix atom-pair frequency with residue-edge frequency.
2. Separate native plausibility from interaction assignment. Failure to form an
   H-bond does not make a chemically valid proximal contact unnatural, and a
   geometrically valid H-bond does not certify a binding model. A low-frequency
   category does not by itself deserve a positive energy penalty.
3. Available published cutoffs define support/typing, not a normal distribution.
   Do not translate a cutoff interval into a guessed mean, variance or prior
   effective sample count. Normal-transformed priors still require quantitative
   reference observations or a source that reports their distribution.
4. Existing `tcren.clashes` uses element-specific Bondi radii and an overlap
   tolerance. PICCOLO uses residue-specific radii and a different clash rule.
   Preserve these definitions separately; importing the latter is a deliberate
   model change, not a correction to historical output.
5. The raw antibody table may support marginal geometry statistics, but the prior
   strength and atom-to-residue aggregation remain explicit modelling choices.
   The positive-energy candidate menu does not resolve those dependencies.

## Recheck recipes and current scope

Use PubMed `get_articles` for PMIDs `38161735`, `21801404`, `27153634`,
`2308168`, `24488795`; `get_full_text_articles` for `PMC10755492` and
`PMC10392943`. PICCOLO's 50000-character retrieval is truncated near the end;
its Table 1 and supplement descriptions were read, not its unreturned tail.
GitHub `fetch_file` at the pinned commits above reproduces the two inspected
producer/configuration files. Resolve DTU DOI `10.11583/DTU.22555672` and
inspect its public file metadata. The successful public download recipes and
hashes above now resolve the original contact-table retrieval dependency.

A selective filename search in the local `~/hf` mirror trees and benchmark
checkout found no matching `Ab_interface`, `22555672`, `PICCOLO`, `Arpeggio`
or antibody-contact filename. This is a limited search, not a global absence
claim. Existing archived docking priors are not atom-contact prior substitutes.
No library score, frozen table, benchmark label, scientific output or MD job
was changed by this dossier.
