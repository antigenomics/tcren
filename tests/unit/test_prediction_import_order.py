"""A real file-based call must work after an immutable caller pool initializes."""
import os
from pathlib import Path
import subprocess
import sys


def test_file_inspection_preserves_preinitialized_polars_pool(tmp_path):
    fixture=Path(__file__).parents[1]/'assets/pdb/1ao7.pdb'
    manifest=tmp_path/'manifest.tsv'
    manifest.write_text('complex.id\tstructure\tgenerator\npublic_fixture\t'+str(fixture)+'\tnative\n')
    code='''
import polars as pl
assert pl.thread_pool_size()==2
from tcren.prediction import inspect_models_command
from pathlib import Path
import sys
inspect_models_command(Path(sys.argv[1]),Path(sys.argv[2]),threads=1)
assert pl.thread_pool_size()==2
rows=pl.read_csv(Path(sys.argv[2])/'models.tsv',separator='\\t')
assert rows['complex.id'].to_list()==['public_fixture']
assert rows['descriptor_status'].to_list()==['available']
assert (Path(sys.argv[2])/'contacts.tsv').stat().st_size>100
'''
    subprocess.run([sys.executable,'-c',code,str(manifest),str(tmp_path/'inspection')],
                   env=dict(os.environ,POLARS_MAX_THREADS='2'),check=True,timeout=180)
