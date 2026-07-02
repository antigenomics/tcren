"""Rescore the Bobisse & Bigot neoantigen benchmarks with the v1 (2022-paper) TCRen potential.

Closes the last legacy-reproduction gap: the v2 numbers already exist
(benchmark_bobisse.csv: cognate rank 1/79; benchmark_bigot.csv: median rank 13/43). This
script reruns the *same* scoring pipelines as 04_benchmark_neoepitopes.ipynb (Bobisse) and
06_benchmark_models.ipynb cell 4 (Bigot Part B) but swaps ONLY the potential matrix for the
v1 potential, so the geometry / candidate sets / thresholds are byte-for-byte identical to v2
and the only moving part is v1-vs-v2 of the energy table.

v1 potential choice
-------------------
The v1/2022 potential is results_new/TCRen_2022.csv (column `value`). This is the table both
02_benchmark_cognate_unrelated.ipynb (cell: 'TCRen 2022': mat_from(read_csv('TCRen_2022.csv')))
and results_new/benchmark_cognate_ranks.csv (row 'TCRen 2022', median 8% / AUC 0.81) treat as
"TCRen 2022" throughout the v2 reproduction. (The separately-shipped legacy table
src/tcren/data/TCRen_potential.csv is a *different* v1 candidate — corr 0.85, max |Δ| 1.6 vs
TCRen_2022.csv — and is not what the notebooks label "TCRen 2022"; per the task's disambiguation
rule we use the notebook/benchmark_cognate_ranks.csv definition.)

v2 potentials (for reference, NOT used here): Bobisse used the 6p64 leave-one-out potential
(TCRen_2026_LOO.csv, col TCRen.LOO); Bigot used the TCRen 2026 full potential (TCRen_2026.csv,
col value). This script substitutes the single v1 matrix into both.

Outputs (in results_new/)
-------------------------
  raw_bobisse_v1.csv   rank, peptide, HLA, pep.len, TCRen.score, real   (matches benchmark_bobisse.csv)
  raw_bigot_v1.csv     pdb.id, cognate, rank, n                          (matches benchmark_bigot.csv)

Run from notebooks/natcompsci2022/ (conda env with tcren + arda importable).
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, polars as pl
from tcren.structure import parse_structure
from tcren.annotation import classify_chains
from tcren.annotation.arda_adapter import _import_arda
from tcren.paper import contact_table, annotate_structure_set
from tcren.paper.helpers import _batch_annotate

R = 'results_new'
AA = list('LFIMVWYCHAGPTSQNDERK'); AIDX = {a: i for i, a in enumerate(AA)}


def mat_from(df, vcol):
    m = np.zeros((20, 20))
    for r in df.iter_rows(named=True):
        i, j = AIDX.get(r['residue.aa.from']), AIDX.get(r['residue.aa.to'])
        if i is not None and j is not None:
            m[i, j] = r[vcol]
    return m


# The single v1 (2022-paper) potential matrix, used for BOTH benchmarks.
M_V1 = mat_from(pl.read_csv(f'{R}/TCRen_2022.csv'), 'value')


# ============================================================================================
# Bobisse — 302TIL neoepitope ranking (mirrors 04_benchmark_neoepitopes.ipynb cells 2-3)
# ============================================================================================
COGNATE = 'KQWLVWLFL'

files = sorted(Path('../data/Bobisse').glob('*.pdb.gz'))
structs = [parse_structure(p, pdb_id=p.stem) for p in files]
records = _batch_annotate(structs, _import_arda())

geom = {}                                    # (HLA, length) -> (contact positions, TCR aa indices)
for idx, (s, p) in enumerate(zip(structs, files)):
    classify_chains(s, organism='human', autodetect_species=True, precomputed_records=records[idx])
    ab = contact_table(s).filter(pl.col('chain.type.from').is_in(['TRA', 'TRB']))
    pos = np.array(ab['pos.to'].to_list())
    tcr = np.array([AIDX.get(a, -1) for a in ab['residue.aa.from'].to_list()])
    keep = (pos >= 0) & (tcr >= 0)
    hla, length = p.stem.split('_')[1], int(p.stem.split('_')[2])
    geom[(hla, length)] = (pos[keep], tcr[keep])

cand = (pl.read_csv('data_legacy/Bobisse/Bobisse_peptides.tsv.gz', separator='\t')
        .select(peptide='Peptide', HLA='HLA')
        .with_columns(HLAn=pl.col('HLA').str.replace_all(':', ''),
                      L=pl.col('peptide').str.len_chars())
        .unique())

rows = []
for r in cand.iter_rows(named=True):
    pep, key = r['peptide'], (r['HLAn'], r['L'])
    if key not in geom or any(a not in AIDX for a in pep):
        continue
    pos, tcr = geom[key]
    aa = np.array([AIDX[a] for a in pep])
    score = float(M_V1[tcr, aa[pos]].sum())
    rows.append({'peptide': pep, 'HLA': r['HLA'], 'pep.len': r['L'],
                 'TCRen.score': score, 'real': pep == COGNATE})
bob = pl.DataFrame(rows).sort('TCRen.score').with_row_index('rank', offset=1)
bob.write_csv(f'{R}/raw_bobisse_v1.csv')

n_bob = bob.height
bob_rank_v1 = bob.filter(pl.col('real'))['rank'][0]


# ============================================================================================
# Bigot — 14 patient TCRs, cognate rank among 44 candidate 9-mers
# (mirrors 06_benchmark_models.ipynb cell 4, Part B)
# ============================================================================================
def geom_ab(ab, pid):
    sub = ab.filter(pl.col('pdb.id') == pid)
    pos = np.array(sub['pos.to'].to_list())
    tcr = np.array([AIDX.get(a, -1) for a in sub['residue.aa.from'].to_list()])
    keep = (pos >= 0) & (tcr >= 0)
    return pos[keep], tcr[keep]


def score_pep(pos, tcr, pep, M):
    keep = pos < len(pep); p, t = pos[keep], tcr[keep]
    aa = np.array([AIDX[a] for a in pep])
    return float(M[t, aa[p]].sum())


cbig, _ = annotate_structure_set('../data/Bigot')
ab_big = cbig.filter(pl.col('chain.type.from').is_in(['TRA', 'TRB']))

candidates = pl.read_csv('data_legacy/Bigot/Bigot_candidate_epitopes.txt.gz')['peptide'].to_list()
cognate = {str(r['pdb.id']): r['cognate.peptide']
           for r in pl.read_csv('data_legacy/Bigot/Bigot_cognate_epitopes.csv.gz').iter_rows(named=True)}

bigot_rows = []
for pid in sorted(ab_big['pdb.id'].unique().to_list()):
    cog = cognate.get(str(pid))
    if not cog:
        continue
    peps = sorted(set(candidates) | {cog})
    pos, tcr = geom_ab(ab_big, pid)
    sc = {p: score_pep(pos, tcr, p, M_V1) for p in peps if all(a in AIDX for a in p)}
    order = sorted(sc, key=sc.get)                          # ascending energy
    rank = order.index(cog) + 1
    bigot_rows.append({'pdb.id': pid, 'cognate': cog, 'rank': rank, 'n': len(sc)})
big = pl.DataFrame(bigot_rows).sort('rank')
big.write_csv(f'{R}/raw_bigot_v1.csv')

big_median_v1 = float(big['rank'].median())
big_n = int(big['n'][0]) if big.height else 0


# ============================================================================================
# v1-vs-v2 comparison
# ============================================================================================
bob_v2 = pl.read_csv(f'{R}/benchmark_bobisse.csv')
bob_rank_v2 = bob_v2.filter(pl.col('real'))['rank'][0]
big_v2 = pl.read_csv(f'{R}/benchmark_bigot.csv')
big_median_v2 = float(big_v2['rank'].median())

print('=== v1 (TCRen 2022) rescore — neoepitope benchmarks ===')
print(f'v1 potential: results_new/TCRen_2022.csv (col "value") — the table the notebooks label "TCRen 2022"')
print()
print('Bobisse (302TIL, cognate KQWLVWLFL):')
print(f'  candidates scored: {n_bob}')
print(f'  cognate rank  v1: {bob_rank_v1}/{n_bob}   |   v2: {bob_rank_v2}/{bob_v2.height}')
print()
print('Bigot (14 patient TCRs):')
print(f'  patients scored: {big.height}   candidates per patient: {big_n}')
print(f'  median cognate rank  v1: {big_median_v1:.0f}/{big_n}   |   v2: {big_median_v2:.0f}/{int(big_v2["n"][0])}')
print()
print('per-patient Bigot ranks (v1):')
print(big)
print()
print(f'wrote {R}/raw_bobisse_v1.csv  ({bob.height} rows)')
print(f'wrote {R}/raw_bigot_v1.csv    ({big.height} rows)')
