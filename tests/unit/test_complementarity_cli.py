import pytest
import typer
from typer.testing import CliRunner

from tcren.complementarity.command import complementarity_command


@pytest.mark.parametrize("arguments", [
    ["--task", "unknown"], ["--origin", "unknown"], ["--cutoff", "nan"],
    ["--cutoff", "0"], ["--cutoff", "11"], ["--threads", "-1"],
    ["--legacy", "--cutoff", "2"], ["--scratch", "/nonexistent/task-scratch"],
])
def test_invalid_options_fail_before_export(tmp_path, arguments):
    app = typer.Typer()
    app.command()(complementarity_command)
    output = tmp_path / "new"
    result = CliRunner().invoke(app, ["-s", "unused.pdb", "-o", str(output), *arguments])
    assert result.exit_code != 0
    assert not output.exists()


def test_existing_output_is_never_overwritten(tmp_path):
    app = typer.Typer()
    app.command()(complementarity_command)
    marker = tmp_path / "preserved.txt"
    marker.write_text("existing user output")
    result = CliRunner().invoke(app, ["-s", "unused.pdb", "-o", str(tmp_path)])
    assert result.exit_code != 0
    assert marker.read_text() == "existing user output"
