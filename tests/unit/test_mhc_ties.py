"""A cropped groove cannot resolve alleles differing outside its sequence."""
from types import SimpleNamespace

import polars as pl

from tcren.mhc.mapper import _COLUMNS, _best_hits, calls_from_hits


def test_mhc_ties_preserved_independent_of_search_order(tmp_path):
    rows = []
    for allele, bits in [('HLA-A*02:608N', 392), ('HLA-A*02:01', 392),
                         ('HLA-A*02:02', 380), ('HLA-A*02:01', 392)]:
        rows.append(['D', f'{allele}|A|MHCI|MHCa|human', 1, 180, 25, 204,
                     180, 365, 180, 0, 0, '180M', 'A'*180, 'A'*180, 1e-100,
                     bits, 100.0])
    path = tmp_path / 'hits.tsv'
    for values in (rows, rows[::-1]):
        pl.DataFrame(values, schema=_COLUMNS, orient='row').write_csv(
            path, separator='\t', include_header=False)
        call, = calls_from_hits([SimpleNamespace(chain_id='D')], _best_hits(path))
        assert call.allele == 'HLA-A*02:01'
        assert call.allele_candidates == ('HLA-A*02:01', 'HLA-A*02:608N')
        assert (call.qstart, call.qend, call.tstart, call.tend) == (1, 180, 25, 204)
