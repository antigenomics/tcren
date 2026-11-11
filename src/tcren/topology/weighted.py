"""Fixed, fit-free energy-weighted functionals of the recognition graph ladder.

Historical defaults use TCRen2 TCR/MHC to peptide and MJ TCR to MHC.
Explicit TP/TM matrices override those defaults; an explicit PM matrix uses
peptide to MHC. Goodness is max(-energy, 0). Unsupported edges remain present
and withhold weighted values.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.linalg import solve_triangular
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, maximum_flow, minimum_spanning_tree, shortest_path

from .graph import _contact_graph_counts, _pielou


def _cut(adjacency, receptor):
    """All-terminal cut; finite super-terminal capacity exceeds the whole graph."""
    n = len(adjacency)
    capacity = np.zeros((n+2, n+2), dtype=np.int64)
    capacity[:n, :n] = adjacency
    ceiling = int(adjacency.sum())+1
    if ceiling > np.iinfo(np.int32).max:
        raise ValueError("flow capacity exceeds SciPy int32 support")
    capacity[n, np.flatnonzero(receptor)] = ceiling
    capacity[np.flatnonzero(~receptor), n+1] = ceiling
    return float(maximum_flow(csr_matrix(capacity), n, n+1).flow_value)


def _betweenness(adjacency):
    """Normalized Brandes dependency reduction on SciPy reciprocal-weight distances."""
    n = len(adjacency)
    result = np.zeros(n)
    lengths = np.divide(1., adjacency, out=np.zeros_like(adjacency), where=adjacency > 0)
    distances = shortest_path(csr_matrix(lengths), directed=False)
    for source in range(n):
        order = np.argsort(distances[source], kind='stable')
        order = order[np.isfinite(distances[source, order])]
        distances_source = distances[source, order]
        selected = lengths[np.ix_(order, order)]
        # One native array comparison per source, not one per source/vertex pair.
        tight = ((selected > 0) & (distances_source[:, None] < distances_source[None, :]) &
                 np.isclose(distances_source[:, None]+selected, distances_source[None, :],
                            rtol=1e-12, atol=1e-12))
        count = len(order)
        seed = np.zeros(count)
        seed[0] = 1
        sigma = solve_triangular(np.eye(count)-tight.T, seed, lower=True, unit_diagonal=True)
        ratio = tight*sigma[:, None]/sigma[None, :]
        dependency = solve_triangular(np.eye(count)-ratio, ratio.sum(axis=1),
                                      lower=False, unit_diagonal=True)
        dependency[0] = 0
        result[order] += dependency
    return result/((n-1)*(n-2)) if n > 2 else result


def _weighted_features(weights, tcr, loops, targets, chains=None):
    """Reduce one canonical edge→goodness mapping; keys may be physical tuples."""
    edges = set(weights)
    vertices = sorted({v for edge in edges for v in edge})
    e, n, c, beta = _contact_graph_counts(edges)
    missing = sum(w is None or not np.isfinite(w) for w in weights.values())
    row = dict(edges=e, vertices=n, beta0=c, beta1=beta,
               available=bool(set(vertices) & tcr), missing_weight=missing,
               density=2*e/(n*(n-1)) if n > 1 else float('nan'),
               edges_per_vertex=e/n if n else float('nan'),
               edges_per_vertex_bounded=e/(e+n) if n else float('nan'))
    names = ('w_per_vertex w_per_vertex_bounded w_total w_mean w_max s_max s_mean s_even w_cyclo giant_frac_q25 '
             'giant_frac_q50 giant_frac_q75 graph_energy graph_energy_c graph_energy_bounded '
             'wmincut_tcr_pmhc wmincut_frac alg_conn_w alg_conn_w_bounded '
             'btw_tcr btw_pmhc partcoef_tcr_w partcoef_pmhc_w').split()
    row.update(dict.fromkeys(names, float('nan')))
    index = {v: i for i, v in enumerate(vertices)}
    a = np.zeros((n, n))
    for u, v in sorted(edges):
        a[index[u], index[v]] = a[index[v], index[u]] = 1
    degree = a.sum(axis=1)
    triangles_at = np.diag(a@a@a)/2
    paths_at = degree*(degree-1)/2
    local = np.divide(triangles_at, paths_at, out=np.zeros(n), where=paths_at > 0)
    row.update(clustering=float(local.mean()) if n else float('nan'),
               transitivity=float(triangles_at.sum()/paths_at.sum()) if paths_at.sum() else 0.)
    receptor = np.array([v in tcr for v in vertices], dtype=bool)
    supported = receptor.any() and (~receptor).any()
    row['mincut_tcr_pmhc'] = _cut(a.astype(np.int64), receptor) if supported else float('nan')
    denominator = min(degree[receptor].sum(), degree[~receptor].sum()) if supported else 0
    row['mincut_frac'] = row['mincut_tcr_pmhc']/denominator if denominator else float('nan')
    if missing or not n:
        return row
    for (u, v), w in sorted(weights.items()):
        if w < 0:
            raise ValueError('goodness weights must be nonnegative')
        a[index[u], index[v]] = a[index[v], index[u]] = w
    strength = a.sum(axis=1)
    total = float(strength.sum()/2)
    energy = float(abs(np.linalg.eigvalsh(a)).sum())
    row.update(w_total=total, w_per_vertex=total/n,
               w_per_vertex_bounded=total/(total+n), w_mean=total/e, w_max=float(a.max()),
               s_max=float(strength.max()), s_mean=float(strength.mean()),
               s_even=_pielou(strength, n), graph_energy=energy,
               graph_energy_bounded=energy/(2*total) if total else float('nan'))
    if chains is None:
        chains = {v: v for v in vertices}
    intra = 0.
    for chain in sorted(set(chains[v] for v in vertices)):
        mask = np.array([chains[v] == chain for v in vertices])
        intra += float(abs(np.linalg.eigvalsh(a[np.ix_(mask, mask)])).sum())
    inter = sum(chains[u] != chains[v] for u, v in edges)
    row['graph_energy_c'] = (energy-intra)/inter if inter else float('nan')
    # Negative costs make the minimum forest maximize retained positive strength.
    tree = minimum_spanning_tree(csr_matrix(-a))
    row['w_cyclo'] = (total+float(tree.sum()))/total if total else float('nan')
    values = np.array([weights[edge] for edge in sorted(edges)])
    for q in (25, 50, 75):
        threshold = np.percentile(values, q)
        keep = (a > 0) & (a >= threshold)
        # Zero-weight edges at a zero quantile remain in the unweighted graph.
        if threshold == 0:
            keep = np.zeros((n, n), dtype=bool)
            for u, v in edges:
                keep[index[u], index[v]] = keep[index[v], index[u]] = True
        _, labels = connected_components(csr_matrix(keep), directed=False)
        row[f'giant_frac_q{q}'] = float(np.bincount(labels).max()/n)
    if a.max() > (np.iinfo(np.int64).max//max(1, 4*e)-1)/100:
        raise ValueError('weighted flow capacities overflow int64')
    capacity = np.rint(100*a).astype(np.int64)
    rounded_strength = capacity.sum(axis=1)
    denominator = min(rounded_strength[receptor].sum(), rounded_strength[~receptor].sum()) if supported else 0
    cut = _cut(capacity, receptor) if supported else float('nan')
    row['wmincut_tcr_pmhc'] = cut/100
    row['wmincut_frac'] = cut/denominator if denominator else float('nan')
    _, labels = connected_components(csr_matrix(a), directed=False)
    keep = labels == np.bincount(labels).argmax()
    sub = a[np.ix_(keep, keep)]
    if len(sub) >= 3:
        d = sub.sum(axis=1)
        laplacian = np.eye(len(sub))-sub/np.sqrt(d[:, None]*d[None, :])
        value = float(np.linalg.eigvalsh(laplacian)[1])
        row.update(alg_conn_w=value, alg_conn_w_bounded=value/2)
    centrality = _betweenness(a)
    for name, mask in [('tcr', receptor), ('pmhc', ~receptor)]:
        if mask.any():
            row['btw_'+name] = float(centrality[mask].mean())
    # Participation modules match literature.py: target class / six CDR loops.
    for name, side, modules in [('tcr', tcr, targets), ('pmhc', set(vertices)-tcr, loops)]:
        participation = []
        for vertex in vertices:
            if vertex not in side:
                continue
            sums = {}
            for other in vertices:
                if other in modules:
                    module = modules[other]
                    sums[module] = sums.get(module, 0.)+a[index[vertex], index[other]]
            s = np.array(list(sums.values()))
            if s.sum() > 0:
                participation.append(float(1-np.square(s/s.sum()).sum()))
        if participation:
            row[f'partcoef_{name}_w'] = float(np.mean(participation))
    return row


def weighted_graph_features(structure, contact_map, *, cutoff=5.0,
                            peptide_potential=None, mhc_potential=None, pmhc_potential=None,
                            weighting="positive", beta=1.0):
    """Existing ladder graphs, explicit role-directed matrices and incident V.

    PM is peptide→MHC when supplied, otherwise the historical MHC→peptide
    lookup. Matrices set scalar weights on undirected physical edges.
    """
    from ..exports import _recognition_graphs
    from ..potential import mj, tcren2

    if weighting not in {"positive", "boltzmann"} or not np.isfinite(beta) or beta <= 0:
        raise ValueError("graph weighting requires positive|boltzmann and finite beta > 0")
    graphs, parts, _, _, _ = _recognition_graphs(structure, contact_map, cutoff=cutoff)
    residues = {(c.chain_id, r.seq_index): r for c in structure.chains for r in c.residues}
    types = {c.chain_id: c.chain_type for c in structure.chains}
    loops = {(c.chain_id, r.seq_index): c.chain_type+':'+reg.region_type
             for c in structure.chains if c.chain_type in ('TRA', 'TRB')
             for reg in c.regions if reg.region_type in ('CDR1', 'CDR2', 'CDR3') for r in reg.residues}
    targets = {v: 'peptide' if types[v[0]] == 'PEPTIDE' else 'mhc'
               for v in residues if types[v[0]] not in ('TRA', 'TRB')}
    # Build tiny 20x20 lookup maps once, avoiding a Polars filter per residue edge.
    potentials = {}
    for name, potential in [
        ('tp', peptide_potential if peptide_potential is not None else tcren2()),
        ('tm', mhc_potential if mhc_potential is not None else mj()),
        ('pm', pmhc_potential if pmhc_potential is not None else
         (peptide_potential if peptide_potential is not None else tcren2()))]:
        table = potential.matrix
        values = table['value'].to_numpy()
        finite = np.isfinite(values)
        goodness = np.full(values.shape, np.nan)
        if weighting == 'boltzmann':
            try:
                with np.errstate(over='raise', under='raise', invalid='ignore'):
                    goodness[finite] = np.exp(-beta*values[finite])
            except FloatingPointError as exc:
                raise ValueError('Boltzmann weights overflow/underflow for supplied beta and matrix') from exc
        else:
            goodness[finite] = np.maximum(-values[finite], 0.)
        potentials[name] = dict(zip(table.select('residue.aa.from', 'residue.aa.to').iter_rows(), goodness))
    weights = {}
    for interface, edges in parts.items():
        for u, v in sorted(edges):
            if interface == 'peptide_mhc':
                if pmhc_potential is not None:
                    source, target = (u, v) if types[u[0]] == 'PEPTIDE' else (v, u)
                else:
                    source, target = (v, u) if types[u[0]] == 'PEPTIDE' else (u, v)
            else:
                source, target = (u, v) if u in loops else (v, u)
            value = potentials[{'tcr_mhc':'tm', 'peptide_mhc':'pm'}.get(interface,'tp')].get(
                (residues[source].aa, residues[target].aa))
            weights[(u, v)] = float(value) if value is not None and np.isfinite(value) else None
    chains = {v: v[0] for v in residues}
    return pl.DataFrame([dict(graph=name, **_weighted_features(
        {edge: weights.get(edge) for edge in sorted(edges)}, set(loops), loops, targets, chains))
        for name, edges in graphs], infer_schema_length=None)
