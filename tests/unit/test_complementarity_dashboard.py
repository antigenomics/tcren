import gzip
import json

import pytest

from tcren.complementarity.dashboard import write_dashboard


def test_streamed_scene_identity_is_never_a_path(tmp_path):
    scene = dict(id="../unsafe</script>", frame=dict(reference="ref"), warnings=[])
    out = tmp_path / "browser"
    path = write_dashboard(iter([scene]), out)
    assert path.is_file()
    assert json.loads(gzip.decompress((out / "scene-0.json.gz").read_bytes())) == scene
    assert json.loads((out / "index.json").read_text())[0]["id"] == scene["id"]
    assert "../unsafe" not in path.read_text()
    with pytest.raises(FileExistsError):
        write_dashboard([], out)


def test_nan_export_fails_explicitly(tmp_path):
    scene = dict(id="bad", frame=dict(reference="ref"), warnings=[], cutoff=float("nan"))
    with pytest.raises(ValueError):
        write_dashboard([scene], tmp_path / "browser")
