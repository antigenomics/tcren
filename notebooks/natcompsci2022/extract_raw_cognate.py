"""Extract per-structure RAW cognate-benchmark scores for the TCRen 2026 (static) potential.

Reproduces the logic of 02_benchmark_cognate_unrelated.ipynb (cells 1-4 for the random-decoy
variant, cells 7-8 for the IEDB-decoy variant) exactly, but saves per-item outputs instead of
only the summary AUC.

Outputs (in results_new/):
  raw_cognate.csv              structure_id, potential, cognate_rank_pct
  raw_cognate_scores.csv       structure_id, is_cognate(1/0), score   (long-form full score vectors)
  raw_cognate_iedb.csv         structure_id, potential, cognate_rank_pct   (IEDB-decoy variant)
  raw_cognate_iedb_scores.csv  structure_id, is_cognate(1/0), score

Run from notebooks/natcompsci2022/.
"""
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, polars as pl

R = 'results_new'
AA = list('LFIMVWYCHAGPTSQNDERK'); AIDX = {a: i for i, a in enumerate(AA)}
POTENTIAL = 'TCRen 2026'  # the static (non-LOO) 2026 potential, per task

rng = np.random.default_rng(0)  # matches notebook cell 3 seed

contacts = pl.read_csv(f'{R}/contacts_2026.csv')
markup = pl.read_csv(f'{R}/markup_2026.csv')
loo = pl.read_csv(f'{R}/TCRen_2026_LOO.csv')
nonred = sorted(loo['pdb.id'].unique().to_list())
peptide_of = {r['pdb.id']: r['peptide'] for r in markup.iter_rows(named=True) if r['peptide']}


def mat_from(df, vcol='value'):
    m = np.zeros((20, 20))
    for r in df.iter_rows(named=True):
        i, j = AIDX.get(r['residue.aa.from']), AIDX.get(r['residue.aa.to'])
        if i is not None and j is not None:
            m[i, j] = r[vcol]
    return m


M = mat_from(pl.read_csv(f'{R}/TCRen_2026.csv'), 'value')
ab = contacts.filter(pl.col('chain.type.from').is_in(['TRA', 'TRB']))

# ---------------------------------------------------------------------------
# Variant 1: random anchor-preserving decoys (notebook cells 3-4)
# ---------------------------------------------------------------------------
N_DECOY = 1000
rank_rows, score_rows = [], []
for pid in nonred:
    cog = peptide_of.get(pid)
    if not cog or len(cog) < 4 or any(a not in AIDX for a in cog):
        continue
    sub = ab.filter(pl.col('pdb.id') == pid)
    pos = np.array(sub['pos.to'].to_list())
    tcr = np.array([AIDX.get(a, -1) for a in sub['residue.aa.from'].to_list()])
    keep = (pos >= 0) & (pos < len(cog)) & (tcr >= 0)
    pos, tcr = pos[keep], tcr[keep]
    if len(pos) == 0:
        continue
    L = len(cog); cogv = np.array([AIDX[a] for a in cog])
    dec = rng.integers(0, 20, size=(N_DECOY, L)); dec[:, 1] = cogv[1]; dec[:, L - 1] = cogv[L - 1]
    allp = np.vstack([cogv[None, :], dec])       # (N+1, L), row 0 = cognate
    contact_aa = allp[:, pos]                     # (N+1, K)
    sc = M[tcr[None, :], contact_aa].sum(axis=1)  # (N+1,)  energy, lower = better binder
    rank_pct = float((sc[1:] < sc[0]).mean() * 100)
    rank_rows.append({'structure_id': pid, 'potential': POTENTIAL, 'cognate_rank_pct': rank_pct})
    # long-form full score vector: row 0 cognate, rest decoys
    score_rows.append({'structure_id': pid, 'is_cognate': 1, 'score': float(sc[0])})
    for s in sc[1:]:
        score_rows.append({'structure_id': pid, 'is_cognate': 0, 'score': float(s)})

pl.DataFrame(rank_rows).write_csv(f'{R}/raw_cognate.csv')
pl.DataFrame(score_rows).write_csv(f'{R}/raw_cognate_scores.csv')

ranks = np.array([r['cognate_rank_pct'] for r in rank_rows])
print('=== raw_cognate.csv (random decoys) ===')
print(f'  n structures         : {len(ranks)}          (target 218)')
print(f'  median rank_pct      : {np.median(ranks):.1f}%      (target 5%)')
print(f'  mean rank_pct        : {np.mean(ranks):.1f}%      (target 12%)')
print(f'  rank-AUC = 1-mean/100: {1 - np.mean(ranks) / 100:.3f}    (target 0.88)')
print(f'  raw_cognate_scores.csv rows: {len(score_rows)}')

# ---------------------------------------------------------------------------
# Variant 2: real IEDB decoys of same MHC group + length (notebook cells 7-8)
# ---------------------------------------------------------------------------
mhc_path = Path(f'{R}/mhc_2026.csv')
if not mhc_path.exists():
    print('\n[skip] IEDB variant: results_new/mhc_2026.csv not found')
    raise SystemExit(0)

mhc = pl.read_csv(mhc_path)


def mhc_norm(s):
    if not s:
        return None
    m = s.split(':')[0]
    if m.startswith('H2-') or m.startswith('H-2'):
        return 'H-2'
    if m.startswith('HLA-D'):
        return m[:6]
    return m


mhc_of = {r['pdb.id']: mhc_norm(r['mhc.allele']) for r in mhc.iter_rows(named=True)}

AA_RE = r'^[LFIMVWYCHAGPTSQNDERK]+$'
iedb = (pl.read_csv('data_legacy/iedb_slim.csv.gz')
        .filter((pl.col('Assay.1') == 'cellular MHC/mass spectrometry') &
                (pl.col('Assay.5') == 'Positive') & (pl.col('Epitope.1') == 'Linear peptide'))
        .with_columns(grp=pl.col('MHC').str.split(':').list.first())
        .with_columns(
            mhcnorm=pl.when(pl.col('grp').str.starts_with('H2-') | pl.col('grp').str.starts_with('H-2'))
                      .then(pl.lit('H-2'))
                      .when(pl.col('grp').str.starts_with('HLA-D')).then(pl.col('grp').str.slice(0, 6))
                      .otherwise(pl.col('grp')),
            pep=pl.col('Epitope.2'))
        .filter(pl.col('pep').str.contains(AA_RE))
        .with_columns(L=pl.col('pep').str.len_chars()))
pools = {(r['mhcnorm'], r['L']): r['peps']
         for r in iedb.group_by(['mhcnorm', 'L']).agg(pl.col('pep').unique().alias('peps')).iter_rows(named=True)}

rng2 = np.random.default_rng(1)  # matches notebook cell 8 seed
rank_rows_iedb, score_rows_iedb = [], []
for pid in nonred:
    cog = peptide_of.get(pid); mn = mhc_of.get(pid)
    if not cog or mn is None or len(cog) < 4 or any(a not in AIDX for a in cog):
        continue
    decoys = [p for p in pools.get((mn, len(cog)), []) if p != cog]
    if len(decoys) <= 100:
        continue
    if len(decoys) > 1000:
        decoys = [decoys[i] for i in rng2.choice(len(decoys), 1000, replace=False)]
    sub = ab.filter(pl.col('pdb.id') == pid)
    pos = np.array(sub['pos.to'].to_list())
    tcr = np.array([AIDX.get(a, -1) for a in sub['residue.aa.from'].to_list()])
    keep = (pos >= 0) & (pos < len(cog)) & (tcr >= 0)
    pos, tcr = pos[keep], tcr[keep]
    if len(pos) == 0:
        continue
    allp = np.array([[AIDX[a] for a in cog]] + [[AIDX[a] for a in d] for d in decoys])  # (N+1, L)
    contact_aa = allp[:, pos]
    sc = M[tcr[None, :], contact_aa].sum(axis=1)
    rank_pct = float((sc[1:] < sc[0]).mean() * 100)
    rank_rows_iedb.append({'structure_id': pid, 'potential': POTENTIAL, 'cognate_rank_pct': rank_pct})
    score_rows_iedb.append({'structure_id': pid, 'is_cognate': 1, 'score': float(sc[0])})
    for s in sc[1:]:
        score_rows_iedb.append({'structure_id': pid, 'is_cognate': 0, 'score': float(s)})

pl.DataFrame(rank_rows_iedb).write_csv(f'{R}/raw_cognate_iedb.csv')
pl.DataFrame(score_rows_iedb).write_csv(f'{R}/raw_cognate_iedb_scores.csv')

ranks_i = np.array([r['cognate_rank_pct'] for r in rank_rows_iedb])
print('\n=== raw_cognate_iedb.csv (real IEDB decoys) ===')
print(f'  n structures         : {len(ranks_i)}          (target 204)')
print(f'  median rank_pct      : {np.median(ranks_i):.1f}%      (target 5%)')
print(f'  mean rank_pct        : {np.mean(ranks_i):.1f}%      (target 13%)')
print(f'  rank-AUC = 1-mean/100: {1 - np.mean(ranks_i) / 100:.3f}    (target 0.87)')
print(f'  raw_cognate_iedb_scores.csv rows: {len(score_rows_iedb)}')
