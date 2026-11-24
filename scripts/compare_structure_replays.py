"""Sequence-verified replay manifest self-check; scientific inputs stay local."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'POLARS_MAX_THREADS'):
    os.environ[name] = '1'

import argparse
import hashlib
import json
from pathlib import Path

from tcren.validation.replay import inspect_replay


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a fresh self-check report')
    manifest = json.loads(args.manifest.read_text())
    ids = [row['id'] for row in manifest['comparisons']]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError('Require a nonempty manifest with unique comparison IDs')
    records = []
    for row in manifest['comparisons']:
        result = inspect_replay(row['reference'], row['model'], row['sequences'],
                                alerts=manifest['alerts'], **row.get('mapping', {}))
        records.append(dict(id=row['id'], result=result))
    packet = dict(manifest_sha256=hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
                  comparisons=records, producer='tcren.validation.replay',
                  scope='Diagnostic checks, not scientific acceptance or rank selection')
    with args.output.open('x') as output:
        json.dump(packet, output, indent=2, allow_nan=False)
        output.write('\n')


if __name__ == '__main__':
    main()
