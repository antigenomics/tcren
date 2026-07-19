# Reproduction validation — re-implementation vs shipped `vdjdb_structure_models`

Evidence that the `scripts/legacy/` re-implementation (tcren structure work + ported legacy
format) reproduces the existing native entries in `isalgo/vdjdb_structure_models`. Regenerate with
`python 10_validate.py`.

## 1. Hash algorithm — exact

`tcr_pmhc_hash = sha256(cdr3.alpha + v.alpha + j.alpha + cdr3.beta + v.beta + j.beta + mhc.a +
mhc.b + antigen.epitope)` (UTF-8, no separators). Fed the *shipped* VDJdb fields, it reproduces
the published hash byte-for-byte:

| record | reproduced hash | matches published |
|---|---|---|
| idx 1159 (predicted, HLA-B*08:01/RPIIRPATL) | `e65469…b14319` | ✅ |
| 1d9k (native, H-2Aa/GNSHRGAIEWEGIESG) | `82945838…37d0d5` | ✅ |

## 2. CDR3 junction — exact

tcren stores IMGT CDR3 (anchors excluded); the VDJdb junction is reconstructed as
`Cys104 + CDR3 + Phe/Trp118`. On 1d9k this reproduces VDJdb exactly: CDR3α `CAATGSFNKLTF`,
CDR3β `CASGGQGRAEQFF`.

## 3. Contacts + `num_contacts` — high overlap on present natives

Per present native: our pipeline output vs the shipped files. `num_contacts` = tcren
`ContactMap.interface("tcr_peptide", tcr_regions="all")` (cutoff 5 Å). Contacts overlap is the
Jaccard of **undirected** residue pairs (the shipped `chain_from`/`chain_to` direction is not
consistent across the dataset — some files list peptide first, others CDR3 first — so direction
is ignored).

| pdbid | class | num_contacts (ours vs shipped) | contacts overlap |
|---|---|---|---|
| 1ao7 | I  | 29 vs 28 | 30/34 = 0.88 |
| 1d9k | II | 18 = 18  | 22/26 = 0.85 |
| 2bnq | I  | 29 vs 25 | 29/31 = 0.94 |
| 2gj6 | II | 21 = 21  | 26/28 = 0.93 |
| 1mi5 | I  | 15 vs 14 | 20/20 = 1.00 |
| 2ckb | I  | 17 = 17  | 11/11 = 1.00 |
| 5c0b | I  | 20 vs 19 | 17/20 = 0.85 |
| 1oga | I  | 19 = 19  | 20/24 = 0.83 |
| 2p5e | I  | 30 vs 26 | 31/31 = 1.00 |
| 3qib | I  | 32 = 32  | 32/32 = 1.00 |

Residual differences are 5 Å-cutoff borderline pairs and minor PDB-parsing differences (tcren vs
the legacy PyMOL-processed files). `num_contacts` matches exactly or differs by a few contacts;
where it differs, ours is the tcren re-count on the primary complex.

## 4. A6 TCR (3D3V)

Already present in the dataset (`aligned_aligned_3d3v.pdb` + contacts + coords), so **not** in the
add set. Our re-implementation recovers the A6 identity correctly: CDR3α `CAVTTDSWGKLQF`,
CDR3β `CASRPGLAGGRPEQYF`. Two caveats specific to this entry (not our pipeline): the shipped 3D3V
has **no metadata row** and a reduced contact set (13 pairs), and the Tax peptide extracts as
`LLFGPVYV` because residue Tyr5 is unresolved in the crystal — so contact overlap here is lower
(13/24 = 0.54) than the well-resolved structures above.

## 5. Orientation caveat (by design)

Coordinates are **not** byte-identical to the shipped files: tcren's canonical frame differs from
the legacy PyMOL frame by a rigid transform. Contacts (distance-based) and `num_contacts` are
frame-invariant and reproduce as shown; each structure's coords/skeleton map are self-consistent.

## 6. Full-run coverage

100 Native2026 structures missing from the dataset → **95 reformatted**, 5 skipped (logged):
`3gjf`, `3hae` (pMHC-only, no TCR); `8yiv`, `8yj2` (β-chain only); `3tf7` (α-chain only). All 5
lack a complete αβ TCR, so they cannot form a valid `tcr_pmhc_hash`.

## 7. Metadata completeness audit

The shipped dataset had **28 native structures with files but no metadata row** (e.g. the A6 TCR
3D3V). Our update resolves them:

- **21 backfilled** (`25_backfill_metadata.py`) — the hash embedded in their existing files is
  present in the metadata (10) or the VDJdb-annotated table (11), so the row is rebuilt from the
  VDJdb record and linked by that exact hash. No new files.
- **2 completed** (`26_complete_orphans.py`, `6uln`/`6ulr`) — no recoverable hash, but they are
  complete αβ TCR:pMHC (KRAS-G12V neoantigen TCRs on HLA-C*08), so identity + hash come from
  annotating the shipped aligned PDB. Metadata-only; `num_contacts` left blank because the complex
  forms under crystal symmetry (the asymmetric unit separates TCR and pMHC).
- **5 excluded** — not TCR:pMHC complexes, so no valid row is possible: `6dxf` (MHC only);
  `6mt3`, `6uli`, `8rnh`, `8rop` (pMHC only, no TCR). These are mis-included in the native archive.

Total metadata rows the update *adds*: **115** (95 new + 19 backfilled + 1 completed, after the
Native2026 rule in §8 drops 3 bogus additions).

**Reverse check — a structure for every entry.** Every metadata row references an existing
structure (predicted → `pdb_files.tgz` by hash; native → `pdb_files_native.tgz` by pdbid): 0
predicted and 0 native rows lack one, both in the shipped metadata and in the final state after the
patch (15,018 rows over 10,884 predicted + 369 native structures). The metadata is therefore
consistent in both directions — every structure has a row (§7) and every row has a structure.

## 8. Native2026 enforcement (`is_native` ⟺ Native2026)

A structure can be native only if it is in the Native2026 source set. `30_assemble.py` enforces
this on the final metadata and native archives, removing **11 bogus natives** (present as native
files/rows but not in Native2026):

- metadata rows removed: `2xn9`, `7l1d`, `7rrg` (3);
- native files removed: 11 pdb + 11 coords + 7 contacts (natives carry no maps; predicted files,
  matched only by native naming, are untouched);
- bogus patch additions never added: `6ulk`, `6uln`, `8rng` (`6uln` is the KRAS structure whose
  partner `6ulr` *is* in Native2026 and is kept).

After enforcement: 0 native rows reference a non-Native2026 structure; no Native2026 native is
dropped. Native rows 294 → 406 (+115 patch − 3 bogus).

## 9. Angles

`scanning_angle`/`pitch_angle` via STCRpy (ANARCI backup models) for **80** of the 118 rows; the
rest are structures STCRpy could not process or metadata-only rows with no local aligned PDB.
