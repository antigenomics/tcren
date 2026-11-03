"""Both orientation routes pass one CPU budget to their batched searches."""
import importlib
import threading
from types import SimpleNamespace

import pytest

from tcren.docking import pipeline


@pytest.mark.parametrize('mode', ['orient', 'superimpose'])
def test_orientation_forwards_threads_to_both_searches(tmp_path, monkeypatch, mode):
    monkeypatch.setattr('tcren.runtime.available_cpus', lambda: 2)
    structs = [SimpleNamespace(chains=[], complex_species='Human') for _ in range(2)]
    files = [tmp_path / '1ao7.pdb', tmp_path / '1bd2.pdb']
    seen = []

    def tcr_batch(items, arda, threads=1):
        seen.append(('tcr', threads, len(items)))
        return [{} for _ in items]

    def mhc_batch(items, threads=1):
        seen.append(('mhc', threads, len(items)))

    monkeypatch.setattr('tcren.annotation.batch._batch_annotate', tcr_batch)
    monkeypatch.setattr('tcren.annotation.arda_adapter._import_arda', lambda: object())
    monkeypatch.setattr('tcren.annotation.classify_chains', lambda *a, **k: None)
    monkeypatch.setattr('tcren.mhc.annotate_mhc_batch', mhc_batch)
    monkeypatch.setattr('tcren.structure.io.structure_paths', lambda path: files)
    monkeypatch.setattr('tcren.structure.io.import_structure', lambda *a, **k: structs[0])
    monkeypatch.setattr('tcren.paths.canonical_dir', lambda: tmp_path)
    monkeypatch.setattr(pipeline, '_write_metadata', lambda *a: None)
    monkeypatch.setattr(pipeline, '_finish_orient', lambda s, pid, *a: {'pdb.id': pid, 'status': 'ok'})
    monkeypatch.setattr(pipeline, '_resolve_inputs', lambda *a: list(zip(['1ao7', '1bd2'], structs)))
    monkeypatch.setattr(pipeline, 'check_oriented_complex', lambda s: (True, None))
    monkeypatch.setattr('tcren.structure.io.write_structure', lambda *a: None)
    superimpose = importlib.import_module('tcren.docking.superimpose')
    parent = threading.get_ident()
    ref_calls = []

    def reference(*args, **kw):
        assert threading.get_ident() == parent, 'reference annotation started in a geometry worker'
        ref_calls.append(kw.get('threads'))
        return []

    monkeypatch.setattr('tcren.docking.align._reference_structure', reference)
    monkeypatch.setattr(superimpose, '_canonical_references', reference)
    monkeypatch.setattr(superimpose, 'superimpose', lambda s, **kw: (
        s, SimpleNamespace(reference_id='reference', rmsd=0.0, n_anchor_atoms=1)))
    function = pipeline.run_folder if mode == 'orient' else pipeline.run_superimpose
    assert function(tmp_path, tmp_path / 'output', threads=2).height == 2
    assert seen == [('tcr', 2, 2), ('mhc', 2, 2)]
    assert ref_calls == ([None] if mode == 'orient' else [2])


def test_canonical_reference_cache_separates_organisms(tmp_path, monkeypatch):
    module = importlib.import_module('tcren.docking.superimpose')
    monkeypatch.setattr(module, '_DB_CACHE', {})
    seen = []

    def load(path, cls, species, organism, threads):
        seen.append((organism, threads))
        return [object()]

    monkeypatch.setattr(module, '_load_canonical_references', load)
    human = module._canonical_references(tmp_path, 'MHCI', 'Human', 'human', threads=2)
    assert module._canonical_references(tmp_path, 'MHCI', 'human', 'human') is human
    mouse = module._canonical_references(tmp_path, 'MHCI', 'Human', 'mouse', threads=2)
    assert mouse is not human
    assert seen == [('human', 2), ('mouse', 2)]
