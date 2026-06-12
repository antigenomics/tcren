"""Unit tests for the SVG complementarity-map renderer (fast — synthetic tables)."""

from __future__ import annotations

import xml.etree.ElementTree as ET

import polars as pl
import pytest

from tcren.viz.svg2d import render_complementarity_map

_NS = "{http://www.w3.org/2000/svg}"


def _markup():
    return pl.DataFrame(
        {
            "structure_id": ["x"] * 4,
            "structure_chain": ["P", "P", "D", "M"],
            "complex_chain": ["peptide", "peptide", "tra", "mhca"],
            "complex_region": ["peptide", "peptide", "cdr3", "mhc_helix_a1"],
            "residue_index": [1, 2, 95, 59],
            "aa_index": [0, 1, 94, 58],
            "aa_len": [9, 9, 110, 180],
            "aa": ["G", "I", "Y", "E"],
            "u": [0.0, 5.0, 2.0, 8.0],
            "v": [0.0, 0.0, 6.0, 3.0],
        }
    )


def _contacts():
    return pl.DataFrame(
        {
            "structure_chain_1": ["D"],
            "structure_chain_2": ["P"],
            "aa_index_1": [94],
            "aa_index_2": [1],
            "min_dist": [3.1],
            "contact_type": ["hydrogen_bond"],
        }
    )


def _ca_contacts():
    return pl.DataFrame(
        {
            "structure_chain_1": ["D"],
            "structure_chain_2": ["M"],
            "aa_index_1": [94],
            "aa_index_2": [58],
            "ca_dist": [7.2],
        }
    )


def test_svg_is_well_formed_with_one_rect_per_residue():
    svg = render_complementarity_map(_markup(), contacts=_contacts(), ca_contacts=_ca_contacts())
    root = ET.fromstring(svg)  # raises if malformed
    rects = root.findall(f".//{_NS}g[@class='residues']/{_NS}g/{_NS}rect")
    assert len(rects) == 4


def test_dashed_contacts_carry_metadata():
    svg = render_complementarity_map(_markup(), contacts=_contacts())
    root = ET.fromstring(svg)
    dashed = root.findall(f".//{_NS}g[@class='contacts-dashed']/{_NS}line")
    assert len(dashed) == 1
    assert dashed[0].get("data-contact-type") == "hydrogen_bond"
    assert dashed[0].get("data-min-dist") == "3.10"


def test_bold_ca_contacts_present():
    svg = render_complementarity_map(_markup(), ca_contacts=_ca_contacts())
    root = ET.fromstring(svg)
    bold = root.findall(f".//{_NS}g[@class='contacts-ca']/{_NS}line")
    assert len(bold) == 1
    assert bold[0].get("data-ca-dist") == "7.20"


def test_residue_metadata_and_title():
    svg = render_complementarity_map(_markup())
    root = ET.fromstring(svg)
    groups = root.findall(f".//{_NS}g[@class='residues']/{_NS}g")
    g0 = groups[0]
    assert g0.get("data-complex-chain") == "peptide"
    assert g0.get("data-aa") == "G"
    assert g0.find(f"{_NS}title") is not None


def test_show_chains_filters_residues():
    # Hide the MHC chain → only peptide + TCR residues remain.
    svg = render_complementarity_map(_markup(), show_chains=["peptide", "tra"])
    root = ET.fromstring(svg)
    rects = root.findall(f".//{_NS}g[@class='residues']/{_NS}g/{_NS}rect")
    assert len(rects) == 3  # 2 peptide + 1 tra (mhca dropped)
    chains = {g.get("data-complex-chain")
              for g in root.findall(f".//{_NS}g[@class='residues']/{_NS}g")}
    assert "mhca" not in chains


def test_backbone_connects_sequence_adjacent_residues():
    # Peptide residues aa_index 0 and 1 are adjacent → one backbone line; the others are not.
    svg = render_complementarity_map(_markup(), draw_backbone=True)
    root = ET.fromstring(svg)
    bb = root.findall(f".//{_NS}g[@class='backbone']/{_NS}line")
    assert len(bb) == 1
    assert bb[0].get("data-chain") == "P"

    assert root.findall(f".//{_NS}g[@class='backbone']") and not ET.fromstring(
        render_complementarity_map(_markup(), draw_backbone=False)
    ).findall(f".//{_NS}g[@class='backbone']")


def test_empty_projection_raises():
    empty = _markup().with_columns(pl.lit(None, dtype=pl.Float64).alias("u"))
    with pytest.raises(ValueError):
        render_complementarity_map(empty)
