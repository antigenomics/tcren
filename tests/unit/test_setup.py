"""Bootstrap must install test requirements before invoking the test runner."""
import os
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.mark.parametrize('tests', [False, True])
def test_setup_selects_declared_test_extra(tmp_path, tests):
    root = tmp_path / 'checkout'
    root.mkdir()
    shutil.copyfile(Path(__file__).resolve().parents[2] / 'setup.sh', root / 'setup.sh')
    tools = tmp_path / 'tools'
    tools.mkdir()
    log = tmp_path / 'commands.txt'
    for name in ['uv', 'cc', 'python']:
        script = tools / name
        script.write_text('''#!/bin/sh
printf '%s\\n' "$0 $*" >> "$BOOTSTRAP_LOG"
if [ "$1" = venv ]; then
    mkdir -p "$2/bin"
    touch "$2/bin/activate"
fi
''')
        script.chmod(0o755)
    env = dict(os.environ, PATH=str(tools) + ':' + os.environ['PATH'],
               BOOTSTRAP_LOG=str(log), TCREN_NO_FETCH='1')
    subprocess.run(['bash', str(root / 'setup.sh'), *(['--tests'] if tests else [])],
                   env=env, check=True, capture_output=True, text=True, timeout=10)
    commands = log.read_text().splitlines()
    assert any('pip install -e ' + str(root) + ('[test]' if tests else '') == row.split(' ', 1)[1]
               for row in commands)
    assert any('import tcren' in row and str(root / 'tests') in row for row in commands) == tests
