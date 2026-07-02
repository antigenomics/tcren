"""Extract per-(TCR, peptide) RAW yeast-display (Birnbaum) benchmark scores.

Reproduces 03_benchmark_yeast_display.ipynb cells 1-4 exactly:
  * leave-cluster-out (LCO) TCRen potential derived per the notebook (drop every non-redundant
    structure sharing a Damerau-Levenshtein cluster with any of the 3 Birnbaum test PDBs),
  * top-50 yeast-enriched binders per TCR as positives (y=1),
  * 1000 anchor-preserving random decoys per TCR as negatives (y=0),
  * score = -energy (higher = stronger binder), threaded onto the test structure's contacts.

NOTE on potential choice: the notebook's per-TCR AUC in benchmark_birnbaum.csv comes from the
'TCRen (LCO)' potential (derived on the fly via derive_tcren), NOT from TCRen_2026_LOO.csv.
This script reproduces that exact LCO potential so the AUCs match the committed CSV.

Output (in results_new/):
  raw_birnbaum.csv   tcr, peptide, score, y

Run from notebooks/natcompsci2022/.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, polars as pl
from tcren.potential import derive_tcren
from rapidfuzz.distance import DamerauLevenshtein
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import roc_auc_score

R = 'results_new'
AA = list('LFIMVWYCHAGPTSQNDERK'); AIDX = {a: i for i, a in enumerate(AA)}
TCR2PDB = {'2b4': '3qib', '226': '3qiu', '5cc7': '4p2r'}

contacts = pl.read_csv(f'{R}/contacts_2026.csv')
markup = pl.read_csv(f'{R}/markup_2026.csv')
loo = pl.read_csv(f'{R}/TCRen_2026_LOO.csv')
nonred = sorted(loo['pdb.id'].unique().to_list())
ab = contacts.filter(pl.col('chain.type.from').is_in(['TRA', 'TRB']))

birnbaum = pl.read_csv('data_legacy/Birnbaum.tsv.gz', separator='\t')

# --- Leave-cluster-out set (notebook cell 2) ------------------------------------------------
ab_ids = [r['pdb.id'] for r in contacts.group_by('pdb.id').agg(
    pl.col('chain.type.from').unique().alias('x')).iter_rows(named=True)
    if set(r['x']) <= {'TRA', 'TRB'}]
mk = markup.filter(pl.col('pdb.id').is_in(ab_ids)).with_columns(
    [pl.col(c).fill_null('') for c in ['cdr3a', 'cdr3b', 'peptide']])
ids = mk['pdb.id'].to_list(); seqs = {c: mk[c].to_list() for c in ['cdr3a', 'cdr3b', 'peptide']}
n = len(ids); D = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        d = sum(DamerauLevenshtein.distance(seqs[c][i], seqs[c][j]) for c in seqs)
        D[i, j] = D[j, i] = d
cl = fcluster(linkage(squareform(D), method='complete'), t=6, criterion='distance')
cluster_of = dict(zip(ids, cl))
excluded_clusters = {cluster_of[p] for p in TCR2PDB.values()}
include_lco = [p for p in nonred if cluster_of.get(p) not in excluded_clusters]


def mat_from(df, vcol='value'):
    m = np.zeros((20, 20))
    for r in df.iter_rows(named=True):
        i, j = AIDX.get(r['residue.aa.from']), AIDX.get(r['residue.aa.to'])
        if i is not None and j is not None:
            m[i, j] = r[vcol]
    return m


M_lco = mat_from(derive_tcren(contacts, include=include_lco).matrix)


# --- Scoring (notebook cell 3) --------------------------------------------------------------
def contacts_for(pid):
    sub = ab.filter(pl.col('pdb.id') == pid)
    pos = np.array(sub['pos.to'].to_list())
    tcr = np.array([AIDX.get(a, -1) for a in sub['residue.aa.from'].to_list()])
    keep = (pos >= 0) & (tcr >= 0)
    return pos[keep], tcr[keep]


def energies(pid, idx, M):
    pos, tcr = contacts_for(pid)
    keep = pos < idx.shape[1]; pos, tcr = pos[keep], tcr[keep]
    return M[tcr[None, :], idx[:, pos]].sum(axis=1)


def to_idx(peps):
    return np.array([[AIDX[a] for a in p] for p in peps])


N_DECOY = 1000
rng = np.random.default_rng(0)  # matches notebook cell 3/4 seed
real_peptides = {}
for tcr in TCR2PDB:
    peps = (birnbaum.filter(pl.col('TCR') == tcr).sort('round_5', descending=True)['peptide'].to_list())
    real_peptides[tcr] = [p for p in peps if all(a in AIDX for a in p)][:50]

# --- ROC / raw scores (notebook cell 4, TCRen (LCO) branch) --------------------------------
rows = []
per_tcr_auc = {}
for tcr, pid in TCR2PDB.items():
    real = to_idx(real_peptides[tcr]); L = real.shape[1]
    cog = real[0]
    dec = rng.integers(0, 20, size=(N_DECOY, L)); dec[:, 1] = cog[1]; dec[:, L - 1] = cog[L - 1]
    idx = np.vstack([real, dec]); y = np.r_[np.ones(len(real)), np.zeros(N_DECOY)]
    score = -energies(pid, idx, M_lco)            # higher = better binder
    per_tcr_auc[tcr] = roc_auc_score(y, score)
    # emit per-item rows: positives first (the top-50 real peptides), then decoys.
    for p, s in zip(real_peptides[tcr], score[:len(real)]):
        rows.append({'tcr': tcr, 'peptide': p, 'score': float(s), 'y': 1})
    for k, s in enumerate(score[len(real):]):
        decoy_seq = ''.join(AA[a] for a in dec[k])
        rows.append({'tcr': tcr, 'peptide': decoy_seq, 'score': float(s), 'y': 0})

pl.DataFrame(rows).write_csv(f'{R}/raw_birnbaum.csv')

# --- Verify against benchmark_birnbaum.csv (TCRen (LCO) rows) --------------------------------
target = {r['TCR']: r['AUC'] for r in
          pl.read_csv(f'{R}/benchmark_birnbaum.csv').filter(pl.col('potential') == 'TCRen (LCO)').iter_rows(named=True)}
print('=== raw_birnbaum.csv (potential: TCRen LCO) ===')
print(f'  rows: {len(rows)}  ({sum(r["y"] for r in rows)} positives / {sum(1 - r["y"] for r in rows)} decoys)')
print('  per-TCR AUC (reproduced-from-raw  vs  committed benchmark_birnbaum.csv):')
for tcr in TCR2PDB:
    print(f'    {tcr:5s}: {per_tcr_auc[tcr]:.3f}  vs  {target.get(tcr)}')
mean_repro = np.mean(list(per_tcr_auc.values()))
print(f'  mean per-TCR AUC: {mean_repro:.3f}  vs committed {np.mean(list(target.values())):.3f}  (task said ~0.89)')

# Sanity: recompute AUC straight from the saved raw CSV to confirm it is self-contained.
raw = pl.read_csv(f'{R}/raw_birnbaum.csv')
print('  self-check AUC recomputed from saved raw_birnbaum.csv:')
for tcr in TCR2PDB:
    sub = raw.filter(pl.col('tcr') == tcr)
    print(f'    {tcr:5s}: {roc_auc_score(sub["y"].to_list(), sub["score"].to_list()):.3f}')
