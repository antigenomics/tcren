"""Unit tests for MHC reference parsing (no network)."""

from __future__ import annotations

from tcren.mhc import imgt
from tcren.mhc.imgt import _two_field


def test_two_field_collapse():
    assert _two_field("A*02:01:01:01") == "A*02:01"
    assert _two_field("DRB1*15:01:01") == "DRB1*15:01"
    assert _two_field("B*07") == "B*07"


def test_parse_human_classifies_loci(tmp_path):
    fasta = tmp_path / "hla.fasta"
    fasta.write_text(
        ">HLA:HLA00001 A*02:01:01:01 365 bp\nMAVMAPRTLLLLL\n"
        ">HLA:HLA09999 DRB1*15:01:01 266 bp\nGDTRPRFLWQ\n"
        ">HLA:HLA08888 DQA1*05:01 255 bp\nEDIVADHVAS\n"
        ">HLA:HLA07777 MICA*001 100 bp\nXXXXX\n"  # non-classical, dropped
    )
    alleles = imgt.parse_human(fasta)
    by_allele = {a.allele: a for a in alleles}
    assert "HLA-A*02:01" in by_allele
    assert by_allele["HLA-A*02:01"].mhc_class == "MHCI"
    assert by_allele["HLA-A*02:01"].chain_role == "MHCa"
    assert by_allele["HLA-DRB1*15:01"].chain_role == "MHCb"
    assert by_allele["HLA-DQA1*05:01"].chain_role == "MHCa"
    assert not any(a.locus == "MICA" for a in alleles)  # MICA dropped


def test_parse_mouse_roles(tmp_path):
    mouse = tmp_path / "mouse.fasta"
    mouse.write_text(
        ">sp|P01901|HA1B_MOUSE H-2 class I histocompatibility antigen, K-B alpha chain "
        "OS=Mus musculus OX=10090 GN=H2-K1 PE=1 SV=1\nMVPCTLLLLLAA\n"
        ">sp|P14483|HB2A_MOUSE H-2 class II histocompatibility antigen, A beta chain "
        "OS=Mus musculus OX=10090 GN=H2-Ab1 PE=1 SV=1\nARDSPEDFV\n"
        ">sp|P14434|HA2B_MOUSE H-2 class II histocompatibility antigen, A-B alpha chain "
        "OS=Mus musculus OX=10090 GN=H2-Aa PE=1 SV=2\nEDDIEADHV\n"
        ">sp|P01887|B2MG_MOUSE Beta-2-microglobulin OS=Mus musculus OX=10090 GN=B2m\nIQKTPQIQV\n"
    )
    human_b2m = tmp_path / "hb2m.fasta"
    human_b2m.write_text(">sp|P61769|B2MG_HUMAN Beta-2-microglobulin GN=B2M\nIQRTPKIQV\n")

    alleles = imgt.parse_mouse(mouse, human_b2m)
    roles = {a.locus: (a.mhc_class, a.chain_role) for a in alleles}
    assert roles["H2-K1"] == ("MHCI", "MHCa")
    assert roles["H2-Ab1"] == ("MHCII", "MHCb")
    assert roles["H2-Aa"] == ("MHCII", "MHCa")
    assert roles["B2m"] == ("MHCI", "B2M")
    assert any(a.species == "human" and a.chain_role == "B2M" for a in alleles)


def test_download_is_atomic_and_concurrent_callers_share_publication(tmp_path, monkeypatch):
    import threading
    import time
    from types import SimpleNamespace

    dest = tmp_path / 'input.fasta'
    calls, errors = [], []

    def fetch(*args, **kwargs):
        calls.append(1)
        assert not dest.exists()
        time.sleep(0.05)
        return SimpleNamespace(text='>synthetic\nACDE\n', raise_for_status=lambda: None)

    monkeypatch.setattr(imgt.requests, 'get', fetch)

    def worker():
        try:
            assert imgt._download('synthetic', dest, False).read_text() == '>synthetic\nACDE\n'
        except Exception as exc:
            errors.append(exc)

    workers = [threading.Thread(target=worker) for _ in range(4)]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join(timeout=5)
    assert not any(worker.is_alive() for worker in workers)
    assert not errors
    assert len(calls) == 1


def test_bad_download_preserves_previous_reference(tmp_path, monkeypatch):
    from types import SimpleNamespace
    import pytest

    dest = tmp_path / 'input.fasta'
    dest.write_text('>previous\nACDE\n')
    monkeypatch.setattr(imgt.requests, 'get', lambda *a, **k: SimpleNamespace(
        text='<html>service failure</html>', raise_for_status=lambda: None))
    with pytest.raises(ValueError, match='not FASTA'):
        imgt._download('synthetic', dest, True)
    assert dest.read_text() == '>previous\nACDE\n'


def test_reference_build_stages_on_destination_filesystem(tmp_path, monkeypatch):
    from tcren.mhc import reference
    import os

    allele = imgt.MhcAllele('HLA-A*01:01', 'A', 'MHCI', 'MHCa', 'human', 'ACDE')
    monkeypatch.setattr(imgt, 'download_human', lambda *a, **k: tmp_path / 'raw.fasta')
    monkeypatch.setattr(imgt, 'parse_human', lambda p: [allele])
    replaced = []
    real_replace = os.replace

    def replace(source, destination):
        assert source.parent.parent == tmp_path
        replaced.append(destination.name)
        real_replace(source, destination)

    monkeypatch.setattr(reference.os, 'replace', replace)
    fasta = reference.build(species=('human',), out_dir=tmp_path)
    assert fasta.read_text() == '>HLA-A*01:01|A|MHCI|MHCa|human\nACDE\n'
    assert replaced == ['metadata.tsv', 'alleles.aa.fasta']
    assert not list(tmp_path.glob('.reference-*'))
