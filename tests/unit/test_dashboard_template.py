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
