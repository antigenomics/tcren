"""Finite H0/H1 Vietoris–Rips persistence on recognition-interface residues."""
from __future__ import annotations

import numpy as np

from ..contactmap import ContactMap
from .graph import _recognition_contact_edges


def _persistence_summary(points, max_edge):
    try:
        import gudhi
    except ImportError as error:
        raise ImportError('persistence requires tcren[topology] (gudhi)') from error
    points = np.asarray(points, dtype=float)
    if not np.isfinite(max_edge) or max_edge <= 0:
        raise ValueError('max_edge must be finite and positive')
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError('points must be finite N x 3 coordinates')
    tree = gudhi.RipsComplex(points=points, max_edge_length=max_edge).create_simplex_tree(max_dimension=2)
    tree.persistence(homology_coeff_field=2, persistence_dim_max=tree.dimension() < 2)
    row = {}
    for dimension in (0, 1):
        intervals = tree.persistence_intervals_in_dimension(dimension)
        finite = intervals[np.isfinite(intervals[:, 1])] if len(intervals) else np.empty((0, 2))
        lengths = finite[:, 1]-finite[:, 0]
        count = len(lengths)
        total = float(lengths.sum())
        probabilities = lengths/total if total else np.empty(0)
        entropy = float(-(probabilities*np.log(probabilities)).sum()/np.log(count)) if count > 1 else float('nan')
        row.update({f'pers{dimension}_bars': count,
                    f'pers{dimension}_infinite': len(intervals)-count,
                    f'pers{dimension}_total': total,
                    f'pers{dimension}_max': float(lengths.max()) if count else float('nan'),
                    f'pers{dimension}_mean': total/count if count else float('nan'),
                    f'pers{dimension}_entropy': entropy,
                    f'pers{dimension}_total_bounded': total/(len(points)*max_edge) if len(points) else float('nan'),
                    f'pers{dimension}_max_bounded': float(lengths.max()/max_edge) if count else float('nan')})
    return row


def persistence_features(structure, *, atom='CA', max_edge=12.0, contact_map=None):
    """CA and CB (Gly CA) finite bars through 12 Å, with missing atoms counted.

    Interface vertices are endpoints of CDR:peptide + peptide:MHC + CDR:MHC.
    Infinite bars remain censored, never assigned a fabricated death time.
    ``atom='CB'`` makes CB primary and supplies a CA companion suffix instead.
    """
    if atom not in ('CA', 'CB'):
        raise ValueError('atom must be CA or CB')
    cm = contact_map if contact_map is not None else ContactMap.from_structure(structure, cutoff=5.)
    parts = _recognition_contact_edges(cm)
    keys = sorted({v for edges in parts.values() for edge in edges for v in edge})
    residues = {(c.chain_id, r.seq_index): r for c in structure.chains for r in c.residues}
    row = dict(vertices=len(keys), max_edge_A=max_edge, available=bool(parts['tcr_peptide'] or parts['tcr_mhc']))
    for selected in (atom, 'CB' if atom == 'CA' else 'CA'):
        coordinates = []
        fallback = 0
        for key in keys:
            residue = residues[key]
            point = residue.ca if selected == 'CA' else residue.cb
            if selected == 'CB' and point is None and residue.aa == 'G':
                point = residue.ca
                fallback += int(point is not None)
            if point is not None:
                coordinates.append(point)
        suffix = '' if selected == atom else '_'+selected.lower()
        row.update({name+suffix: value for name, value in _persistence_summary(
            np.asarray(coordinates).reshape(-1, 3), max_edge).items()})
        row.update({'points'+suffix: len(coordinates), 'missing_atoms'+suffix: len(keys)-len(coordinates),
                    'gly_ca_fallback'+suffix: fallback})
    return row
