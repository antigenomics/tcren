# Separate table scoring from coordinate/contact reconstruction on the existing fixture.
import cProfile, pstats, timeit
from tcren.structure.io import import_structure
from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.contactmap import ContactMap
from tcren.pipeline import run, _resolve_potentials
from tcren.energetics.scoring import _interface_energy
s=import_structure('tests/assets/pdb/1ao7.pdb')
classify_chains(s)
s.mhc_calls=annotate_mhc(s)
cm=ContactMap.from_structure(s)
pots=_resolve_potentials(None)
tables={i:cm.interface(i) for i in ('tcr_peptide','tcr_mhc','peptide_mhc')}
def table_scores():
    return {i:_interface_energy(t,pots[i]) for i,t in tables.items()}
def pipeline():
    return run(s,typed=True,superimpose=False,reference_aa='A',reference_mode='structural')
print('contacts',cm.contacts.height, 'table_scores_ms', min(timeit.repeat(table_scores,number=1000,repeat=3))/1000*1000,flush=True)
print('pipeline_ms',min(timeit.repeat(pipeline,number=30,repeat=3))/30*1000,flush=True)
p=cProfile.Profile();p.runcall(lambda:[pipeline() for _ in range(30)])
pstats.Stats(p).strip_dirs().sort_stats('cumulative').print_stats(24)
