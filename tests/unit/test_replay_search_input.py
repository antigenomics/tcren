import json
import pytest
from tcren.validation.prediction_inputs import prepare_af3_search_input


def test_partial_search_reuses_only_identical_declared_chain(tmp_path):
    sequences={'A':'ACDEFG','B':'HIKLMN','C':'PQRSTV','D':'WYADFH'}
    source=tmp_path/'source.json'
    source.write_text(json.dumps(dict(dialect='alphafold3',version=3,name='fixture',modelSeeds=[42],sequences=[
        dict(protein=dict(id=k,sequence=v,unpairedMsa='>query\n'+v+'\n',pairedMsa='>query\n'+v+'\n',templates=[]))
        for k,v in sequences.items()])))
    expected=dict(sequences,A='ACDEFA',C='PQRST')
    output=tmp_path/'next.json'
    report=prepare_af3_search_input(source,output,expected,name='next',seeds=[42])
    entities={x['protein']['id']:x['protein'] for x in json.loads(output.read_text())['sequences']}
    assert set(entities['A'])=={'id','sequence'}
    assert set(entities['B'])=={'id','sequence'}
    assert entities['C']['unpairedMsa']==entities['C']['pairedMsa']==''
    assert entities['D']['unpairedMsa']=='>query\nWYADFH\n'
    assert report['decisions'][-1]['decision']=='reused_exact_sequence_prepared_chain'
    expected['D']='WYADFA'
    other=tmp_path/'other.json'
    prepare_af3_search_input(source,other,expected,name='other',seeds=[42])
    assert set(json.loads(other.read_text())['sequences'][-1]['protein'])=={'id','sequence'}
    with pytest.raises(FileExistsError):
        prepare_af3_search_input(source,output,expected,name='next',seeds=[42])
