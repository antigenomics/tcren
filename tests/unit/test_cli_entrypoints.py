"""The two entry points must expose the same commands.

`tcren` (the console script) imports the module and then calls `main()`, so every `@app.command`
in the file is registered before the app runs. `python -m tcren.cli` executes the module top to
bottom and hits the `__main__` guard where it sits, so any command defined BELOW the guard is
invisible to it.

That has now happened twice, with `potts` and then with `explain`, and neither was caught because
every test and every benchmark producer used one entry point or the other, never both. This test
compares them, which is the only check that fails on a third recurrence rather than shipping it.
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest


def _commands(help_text: str) -> set[str]:
    """The command names out of a Typer/Rich help screen.

    Rich draws the command panels as box-drawing tables, so a command line looks like
    `│ score   Score candidate epitopes ... │`. Take the first token of each such row and drop the
    rows that are continuations of a previous description (they are indented past the name column).
    """
    names: set[str] = set()
    for line in help_text.splitlines():
        m = re.match(r"^[│|]\s(\S+)\s\s+\S", line)
        if m and not m.group(1).startswith("-"):
            names.add(m.group(1))
    return names


def _run(argv: list[str]) -> str:
    env_help = subprocess.run(argv + ["--help"], capture_output=True, text=True,
                              env={"COLUMNS": "200", "TERM": "dumb", "PATH": "/usr/bin:/bin"})
    return env_help.stdout + env_help.stderr


def _console_script() -> Path | None:
    """The `tcren` that belongs to THIS interpreter, not the first one on PATH.

    `shutil.which` is wrong here and finding that out is half the value of this test: on the
    development machine it resolves to a stale `uv tool` install several versions behind the
    checkout, so the comparison reports a difference between two DIFFERENT builds rather than
    between two entry points of one build. The release checklist warns about the same stale
    install for the same reason.
    """
    candidate = Path(sys.executable).parent / "tcren"
    return candidate if candidate.exists() else None


def test_module_and_console_script_expose_the_same_commands():
    console = _console_script()
    if console is None:
        pytest.skip("no console script beside this interpreter; install the package to run this")
    console = str(console)
    module_cmds = _commands(_run([sys.executable, "-m", "tcren.cli"]))
    script_cmds = _commands(_run([console]))
    assert module_cmds, "parsed no commands from the module form; the help layout changed"
    assert module_cmds == script_cmds, (
        "the two entry points disagree. Commands missing from `python -m tcren.cli`: "
        f"{sorted(script_cmds - module_cmds)}. This means the __main__ guard in cli.py sits above "
        "a later @app.command; move it back to the very end of the module."
    )


def test_guard_is_the_last_statement_in_cli():
    """A direct check on the cause, so the failure names the fix even without a console script."""
    from pathlib import Path

    import tcren.cli

    src = Path(tcren.cli.__file__).read_text()
    guard = src.index('if __name__ == "__main__":')
    assert "@app.command" not in src[guard:], (
        "a command is defined below the __main__ guard in cli.py, so `python -m tcren.cli` cannot "
        "see it; move the guard to the very end of the module"
    )
