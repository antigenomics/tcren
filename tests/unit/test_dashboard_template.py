from io import BytesIO

import pytest

from tcren.complementarity import template


def test_bootstrap_reads_the_pinned_revision(tmp_path, monkeypatch):
    lock = template._lock()
    paths = []

    def fetch(url, *, timeout):
        name = url.rsplit("/", 1)[-1]
        assert f"/{lock['revision']}/assets/" in url
        paths.append(name)
        return BytesIO((template.ASSETS / name).read_bytes())

    monkeypatch.setattr(template, "urlopen", fetch)
    output = tmp_path / "template"
    assert template.bootstrap_template(output) == lock["revision"]
    assert set(paths) == set(lock["sha256"])
    template.verify_template(output)
    with pytest.raises(ValueError, match="already exists"):
        template.bootstrap_template(output)


def test_bootstrap_rejects_mutated_assets_atomically(tmp_path, monkeypatch):
    monkeypatch.setattr(template, "urlopen", lambda *args, **kwargs: BytesIO(b"invalid"))
    output = tmp_path / "template"
    with pytest.raises(ValueError, match="checksum mismatch"):
        template.bootstrap_template(output)
    assert not output.exists()


def test_bootstrap_from_checkout(tmp_path):
    output = tmp_path / "template"
    template.bootstrap_template(output, template.ASSETS)
    template.verify_template(output)


def test_shared_style_is_installed_for_a_prepared_viewer(tmp_path):
    from tcren.complementarity.dashboard import install_viewer_assets

    original = tmp_path / "prepared.html"
    original.write_text("<style>body{margin:8px}</style><h1>Prepared models</h1>")
    source_style = tmp_path / "prepared.css"
    source_style.write_text("/* prepared */")
    output = tmp_path / "output"
    output.mkdir()
    (output / "dashboard.css").symlink_to(source_style)
    install_viewer_assets(output, template=original)
    assert source_style.read_text() == "/* prepared */"
    assert (output / "dashboard.css").read_bytes() == (template.ASSETS / "dashboard.css").read_bytes()
    html = (output / "index.html").read_text()
    assert html.count('href="dashboard.css"') == 1
    assert "<style>body{margin:8px}</style>" in html
    assert 'id="view-controls"' in html


def test_reinstalling_prepared_shared_viewer_does_not_duplicate_chrome(tmp_path):
    from tcren.complementarity.dashboard import install_viewer_assets

    install_viewer_assets(tmp_path)
    index = tmp_path / "index.html"
    original = index.read_bytes()
    install_viewer_assets(tmp_path, template=index)
    assert index.read_bytes() == original
