"""Structural scans validate allocation budgets before packing or native work."""
import pytest

from tcren.energetics import background


@pytest.mark.parametrize('threads', [-1])
def test_scan_rejects_invalid_budget_before_packing(monkeypatch, threads):
    monkeypatch.setattr('tcren.runtime.available_cpus', lambda **kwargs: 2)
    monkeypatch.setattr(background, '_pack', lambda *a: pytest.fail('packed before validation'))
    with pytest.raises(ValueError):
        background.scan_structure(None, n_threads=threads)
