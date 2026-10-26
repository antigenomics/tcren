"""Declared alleles are context, not arbitrary search representatives."""
from types import SimpleNamespace

import pytest

from tcren.mhc.mapper import MhcCall
from tcren.mhc.metadata import apply_metadata, set_supplied_allele
from tcren.metadata import read_metadata


def call(role='MHCa', cls='MHCI', allele='HLA-A*02:608N'):
    return MhcCall('D', role, cls, allele, 'A', 'human', 100, 20,
                   1, 3, 1, 3, '3M', (allele,))


def test_preserve_declared_resolution_with_exact_cropped_support():
    result = call()
    set_supplied_allele(result, 'A*02', 'ABC', [
        ('HLA-A*02:01|A|MHCI|MHCa|human', 'XABCY')])
    assert result.allele == 'A*02'
    assert result.representative_allele == 'HLA-A*02:608N'
    assert result.allele_source == 'supplied'
    assert result.supplied_support == 'exact_observed_sequence'


def test_mismatch_and_role_conflict_are_exposed():
    result = call()
    with pytest.warns(UserWarning, match='observed_sequence_not_exact'):
        set_supplied_allele(result, 'HLA-A*02:01', 'XYZ', [
            ('HLA-A*02:01|A|MHCI|MHCa|human', 'ABC')])
    with pytest.raises(ValueError, match='conflicts'):
        set_supplied_allele(result, 'DRB1*01:01', 'ABC', [])


def test_class_ii_alias_metadata_and_conflict(tmp_path, monkeypatch):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>HLA-DRA*01:01|DRA|MHCII|MHCa|human\nABC\n'
                     '>HLA-DRB1*01:01|DRB1|MHCII|MHCb|human\nABC\n')
    monkeypatch.setattr('tcren.mhc.metadata.reference.reference_fasta', lambda: fasta)
    a, b = call('MHCa', 'MHCII', 'HLA-DRA*01:01'), call('MHCb', 'MHCII', 'HLA-DRB1*01:01')
    b.chain_id = 'E'
    chains = {key: SimpleNamespace(sequence=lambda: 'ABC', allele_info=None) for key in ('D', 'E')}
    structure = SimpleNamespace(pdb_id='001', mhc_calls=[a, b], chain=chains.__getitem__)
    path = tmp_path / 'metadata.tsv'
    path.write_text('complex.id\tmhc\n001\tDRA*01:01/DRB1*01:01\n')
    assert read_metadata(tmp_path)['id'].to_list() == ['001']
    apply_metadata([structure], tmp_path)
    assert a.allele == 'DRA*01:01'
    assert b.allele == 'DRB1*01:01'
    path.write_text('id\tmhc.b\n001\tDRA*01:01\n')
    with pytest.raises(ValueError, match='conflicts'):
        apply_metadata([structure], tmp_path)


def test_reannotation_does_not_turn_representative_into_declaration():
    from tcren.mhc.mapper import apply_mhc_calls
    old = call()
    chain = SimpleNamespace(chain_id='D', allele_info=old.allele)
    structure = SimpleNamespace(chains=[chain], mhc_calls=[old])
    new = call(allele='HLA-A*02:01')
    apply_mhc_calls(structure, [new])
    assert new.allele_source == 'sequence_representative'
    assert chain.allele_info == 'HLA-A*02:01'


def test_metadata_rejects_conflicting_key_aliases(tmp_path):
    (tmp_path / 'metadata.tsv').write_text('id\tcomplex.id\n001\t002\n')
    with pytest.raises(ValueError, match='conflicting id'):
        read_metadata(tmp_path)


def test_invalid_allele_is_not_authoritative():
    with pytest.raises(ValueError, match='malformed'):
        set_supplied_allele(call(), 'unknown annotation', 'ABC', [])


def test_reference_headers_parsed_once_per_batch(tmp_path, monkeypatch):
    from tcren.mhc import reference
    fasta=tmp_path/'ref.fa'
    fasta.write_text('>HLA-A*02:01|A|MHCI|MHCa|human\nABC\n')
    monkeypatch.setattr(reference,'reference_fasta',lambda:fasta)
    original=reference.parse_header
    parsed=[]
    monkeypatch.setattr(reference,'parse_header',lambda header:(parsed.append(header),original(header))[1])
    (tmp_path/'metadata.tsv').write_text('id\tmhc\n' + ''.join(f'{i}\tA*02\n' for i in range(10)))
    structures=[]
    for i in range(10):
        chain=SimpleNamespace(sequence=lambda:'ABC',allele_info=None)
        structures.append(SimpleNamespace(pdb_id=str(i),mhc_calls=[call()],chain=lambda _,c=chain:c))
    apply_metadata(structures,tmp_path)
    assert len(parsed)==1
