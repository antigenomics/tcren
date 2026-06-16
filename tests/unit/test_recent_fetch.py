"""RCSB fetch into pdb_recent: download + 5-chain validation (network-gated)."""

from __future__ import annotations

import pytest

requests = pytest.importorskip("requests")
pytest.importorskip("arda")  # fetch_ids validates via batched annotation (needs arda/mmseqs)

from tcren.recent import _download_cif_gz, fetch_ids


def test_download_and_validate_complete_complex(tmp_path):
    # 1ao7 is a complete MHCI TCR-pMHC complex; must download as .cif.gz and validate.
    p = _download_cif_gz("1ao7", tmp_path)
    if p is None:
        pytest.skip("RCSB unreachable")
    assert p.exists() and p.name == "1ao7.cif.gz"

    summary = fetch_ids(["1ao7"], dest=tmp_path)
    assert summary["downloaded"] == 1
    assert summary["complete"] == 1 and summary["kept"] == 1  # has all 5 required chains
    assert (tmp_path / "1ao7.cif.gz").exists()


def test_extended_id_uses_cif(tmp_path):
    # mmCIF endpoint must be used (PDB deprecates split .pdb); a bogus id yields no file.
    assert _download_cif_gz("zzzz9", tmp_path) is None
