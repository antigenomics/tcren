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
