"""Extract per-pair RAW real-vs-shuffled-structure (Fig 3) benchmark scores.

Reproduces 05_benchmark_shuffle_structures.ipynb cells 1-3 exactly for the TCRen 2026 potential:
score every same-shape structure pair with structure-1's geometry+peptide and structure-2's CDR3
residues; real=1 when the CDR3 donor is the native structure, real=0 for a shuffled donor. eT is
the total TRA+TRB energy; eA/eB the per-chain energies.

Output (in results_new/):
  raw_shuffle.csv   pair_id, real, e_total, e_tra, e_trb

Run from notebooks/natcompsci2022/.
"""
import warnings; warnings.filterwarnings('ignore')
from collections import defaultdict
import numpy as np, polars as pl
from sklearn.metrics import roc_auc_score

R = 'results_new'
AA = list('LFIMVWYCHAGPTSQNDERK'); AIDX = {a: i for i, a in enumerate(AA)}

contacts = pl.read_csv(f'{R}/contacts_2026.csv')
markup = pl.read_csv(f'{R}/markup_2026.csv')
nonred = set(pl.read_csv(f'{R}/TCRen_2026_LOO.csv')['pdb.id'].unique().to_list())


def mat_from(df, vcol='value'):
    m = np.zeros((20, 20))
    for r in df.iter_rows(named=True):
        i, j = AIDX.get(r['residue.aa.from']), AIDX.get(r['residue.aa.to'])
        if i is not None and j is not None:
            m[i, j] = r[vcol]
    return m


M = mat_from(pl.read_csv(f'{R}/TCRen_2026.csv'))  # TCRen only (target = TCRen 0.73)


def _idx(seq):
    return np.array([AIDX[a] for a in seq]) if seq and all(a in AIDX for a in seq) else None


cdr3, lens = {}, {}
for r in markup.iter_rows(named=True):
    if r['pdb.id'] not in nonred:
        continue
    a, b, p = _idx(r['cdr3a']), _idx(r['cdr3b']), r['peptide']
    if a is None or b is None or not p:
        continue
    cdr3[r['pdb.id']] = {'TRA': a, 'TRB': b}
    lens[r['pdb.id']] = (len(a), len(b), len(p))

# CDR3 contact template per structure (notebook cell 2).
cdr3c = contacts.filter((pl.col('region.type.from') == 'CDR3') &
                        (pl.col('chain.type.from').is_in(['TRA', 'TRB'])) &
                        (pl.col('pdb.id').is_in(list(cdr3))))
templ = {pid: {'TRA': ([], []), 'TRB': ([], [])} for pid in cdr3}
for r in cdr3c.iter_rows(named=True):
    pid, ch = r['pdb.id'], r['chain.type.from']
    if r['pos.from'] is None or r['residue.aa.to'] not in AIDX:
        continue
    templ[pid][ch][0].append(r['pos.from']); templ[pid][ch][1].append(AIDX[r['residue.aa.to']])
templ = {pid: {ch: (np.array(v[0]), np.array(v[1])) for ch, v in d.items()} for pid, d in templ.items()}

groups = defaultdict(list)
for pid, shape in lens.items():
    groups[shape].append(pid)
groups = {k: v for k, v in groups.items() if len(v) >= 2}

# Score every same-shape pair (notebook cell 3). Row order preserved -> pair_id = row index.
rows = []
for members in groups.values():
    A = np.stack([cdr3[p]['TRA'] for p in members])   # (k, La) CDR3-alpha residues per donor
    B = np.stack([cdr3[p]['TRB'] for p in members])   # (k, Lb)
    for i, pid1 in enumerate(members):                # pid1 = geometry+peptide provider
        (posA, pepA), (posB, pepB) = templ[pid1]['TRA'], templ[pid1]['TRB']
        ea = M[A[:, posA], pepA[None, :]].sum(axis=1) if len(posA) else np.zeros(len(members))
        eb = M[B[:, posB], pepB[None, :]].sum(axis=1) if len(posB) else np.zeros(len(members))
        et = ea + eb
        for j, pid2 in enumerate(members):            # pid2 = CDR3 donor
            rows.append({
                'pair_id': len(rows),
                'template_id': pid1,
                'cdr3_donor_id': pid2,
                'real': int(j == i),
                'e_total': float(et[j]),
                'e_tra': float(ea[j]),
                'e_trb': float(eb[j]),
            })

pl.DataFrame(rows).write_csv(f'{R}/raw_shuffle.csv')

real = np.array([r['real'] for r in rows])
eT = np.array([r['e_total'] for r in rows])
eA = np.array([r['e_tra'] for r in rows])
eB = np.array([r['e_trb'] for r in rows])
target = {r['potential']: r['AUC'] for r in pl.read_csv(f'{R}/benchmark_shuffle_auc.csv').iter_rows(named=True)}
print('=== raw_shuffle.csv (potential: TCRen 2026) ===')
print(f'  pairs: {len(rows)}  ({int(real.sum())} real / {int((1 - real).sum())} shuffled)')
print(f'  AUC(-e_total) : {roc_auc_score(real, -eT):.3f}   vs committed benchmark_shuffle_auc.csv TCRen {target.get("TCRen")}')
print(f'  AUC(-e_tra)   : {roc_auc_score(real, -eA):.3f}   (Fig 3c TRA)')
print(f'  AUC(-e_trb)   : {roc_auc_score(real, -eB):.3f}   (Fig 3c TRB)')
