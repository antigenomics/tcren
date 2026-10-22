"""A sampled surface patch is buried once, regardless of competing atoms."""
import numpy as np
import pytest

from tcren import _geom


def test_surface_partition_conservation_and_sasa_oracle():
    xyz = np.array([[0., 0., 0.], [0., 3., 0.], [3., 0., 0.], [3., 3., 0.], [40., 0., 0.]])
    radii = np.array([1.7, 1.52, 1.55, 1.7, 1.7])
    roles = np.array([0, 0, 1, 1, 0], dtype=np.int32)
    d = _geom.surface_partition(xyz, radii, roles)
    assert np.allclose(d['bound_A2'], _geom.shrake_rupley(xyz, radii), rtol=0, atol=1e-12)
    for side in (0, 1):
        mask = roles == side
        assert np.allclose(d['isolated_A2'][mask], _geom.shrake_rupley(xyz[mask], radii[mask]), rtol=0, atol=1e-12)
    allocated = np.bincount(d['from'], weights=d['area_A2'], minlength=len(xyz))
    assert np.allclose(allocated, d['isolated_A2']-d['bound_A2'], rtol=0, atol=1e-12)
    assert all(roles[i] != roles[j] for i, j in zip(d['from'], d['to']))
    assert allocated[-1] == 0
    reverse = _geom.surface_partition(xyz, radii, 1-roles)
    assert np.array_equal(reverse['area_A2'], d['area_A2'])
    assert reverse['from'] == d['from'] and reverse['to'] == d['to']


def test_surface_partition_exact_ties_and_duplicate_partner():
    # One sphere point at (3.1,0,0): target surfaces at equal signed gap.
    xyz = np.array([[0., 0., 0.], [4., 1., 0.], [4., -1., 0.]])
    rad = np.full(3, 1.7); roles = np.array([0, 1, 1], dtype=np.int32)
    d = _geom.surface_partition(xyz, rad, roles, n_points=1)
    amounts = [a for i, a in zip(d['from'], d['area_A2']) if i == 0]
    assert amounts == pytest.approx([4*np.pi*3.1**2/2]*2)
    one = _geom.surface_partition(xyz[:2], rad[:2], roles[:2], n_points=1)
    assert sum(amounts) == pytest.approx(sum(a for i,a in zip(one['from'],one['area_A2']) if i==0))
    duplicate = _geom.surface_partition(np.vstack([xyz[:2],xyz[1]]),rad,roles,n_points=1)
    assert sum(a for i,a in zip(duplicate['from'],duplicate['area_A2']) if i==0) == pytest.approx(sum(amounts))
    with pytest.raises(TypeError):
        _geom.surface_partition(xyz, rad, np.array([0.5, 1., 1.]))
    for args in ((np.zeros((2,2)),rad,roles), (xyz,rad,roles[:2]),
                 (xyz,rad,np.array([0,1,2],np.int32)), (xyz,rad*np.nan,roles)):
        with pytest.raises(ValueError):
            _geom.surface_partition(*args)


def test_surface_export_keeps_noncontact_denominator(tmp_path, monkeypatch):
    import json
    import polars as pl
    import tcren.annotation
    import tcren.annotation.batch
    import tcren.mhc
    import tcren.structure
    from tcren.exports import export_batch
    from test_exports import typed_structure, PDB

    monkeypatch.setattr(tcren.structure, 'parse_structure', lambda *a,**k: typed_structure())
    monkeypatch.setattr(tcren.annotation.batch, 'annotate_batch', lambda ss,**k: [{} for _ in ss])
    monkeypatch.setattr(tcren.annotation, 'classify_chains', lambda *a,**k: None)
    monkeypatch.setattr(tcren.mhc, 'annotate_mhc_batch', lambda *a,**k: None)
    monkeypatch.setattr('tcren.exports.annotation_lineage', lambda *a: {})
    output = tmp_path/'surface.parquet'
    ledger = export_batch([str(PDB)], output, quantity='surface-contacts',
                          peptide_potential=None, mhc_potential=None)
    assert ledger['status'].to_list() == ['ok'], ledger
    edges = pl.read_parquet(output)
    atoms = pl.read_parquet(str(output)+'.atoms.parquet')
    assert atoms.filter(pl.col('partner_burial_A2') == 0).height > 0
    assert edges['role.from'].unique().to_list() == ['binder']
    assert edges['role.to'].unique().to_list() == ['target']
    for role, column in [('binder','binder_area_A2'),('target','target_area_A2')]:
        assert edges[column].sum() == pytest.approx(atoms.filter(pl.col('role')==role)['partner_burial_A2'].sum())
    meta = json.loads((tmp_path/'surface.parquet.provenance.json').read_text())
    assert meta['surface_allocation']['sphere'] == 'native Fibonacci-y'
    assert meta['potentials'] == {}


def test_accessible_boundary_excludes_enclosed_atom_and_preserves_defaults():
    from tcren.contacts.exposure import _boundary_allocation

    xyz = np.array([[-2., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
    radii = np.array([1., 2., .5])
    roles = np.array([0, 1, 1], dtype=np.int32)
    old = _geom.surface_partition(xyz, radii, roles, probe=0.)
    sampled = _geom.surface_partition(xyz, radii, roles, probe=0., samples=True)
    for key in old:
        assert np.array_equal(old[key], sampled[key])
    source, target, area = _boundary_allocation(sampled, roles)
    assert sampled['isolated_A2'][2] == 0 and 2 not in target
    assert np.all(roles[source] != roles[target])
    assert np.allclose(np.bincount(source, weights=area, minlength=3),
        sampled['isolated_A2']-sampled['bound_A2'], rtol=0, atol=1e-12)
    assert np.all(sampled['isolated_A2'][target] > 0)


def test_accessible_boundary_exact_ties_and_no_contact():
    from tcren.contacts.exposure import _boundary_allocation

    xyz = np.array([[0., 0., 0.], [1.5, 0., 0.], [1.5, 0., 0.]])
    roles = np.array([0, 1, 1], dtype=np.int32)
    result = _geom.surface_partition(xyz, np.ones(3), roles, probe=0., n_points=1, samples=True)
    source, target, area = _boundary_allocation(result, roles)
    assert source.tolist() == [0, 0] and target.tolist() == [1, 2]
    assert area == pytest.approx([2*np.pi, 2*np.pi])
    xyz[1:, 0] = 20.
    result = _geom.surface_partition(xyz, np.ones(3), roles, probe=0., n_points=1, samples=True)
    assert all(len(x)==0 for x in _boundary_allocation(result, roles))
