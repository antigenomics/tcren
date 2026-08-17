"""Unit tests for pMHC surface topology (:mod:`tcren.surface`)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.surface import (
    DEFAULT_EXTENT,
    KYTE_DOOLITTLE,
    SOURCE_CODES,
    SurfaceMap,
    _groove_frame,
    surface_distance,
    surface_map,
    surface_stats,
    surface_table,
)

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"


@pytest.fixture(scope="module")
def annotated():
    """Chain-typed, MHC-annotated fixtures: one class I, one class II."""
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.structure import parse_structure

    out = {}
    for pid in ("1ao7", "4ozg"):
        s = parse_structure(PDB_DIR / f"{pid}.pdb")
        classify_chains(s, organism="human", autodetect_species=True)
        annotate_mhc(s)
        s.pdb_id = pid
        out[pid] = s
    return out


# --- the ray-cast height kernel -------------------------------------------------------------------
def test_height_candidates_recover_a_single_sphere_exactly():
    """One atom at the origin: the height above the centre is its radius, and the profile is the
    sphere's, not an approximation of it."""
    from tcren.surface import _height_candidates

    grid, extent = (40, 40), (-10.0, 10.0, -10.0, 10.0)
    radius = 4.0
    cell, _owner, height = _height_candidates(np.zeros((1, 3)), np.array([radius]), grid, extent)

    n_x = grid[1]
    cx = _centres_ref(extent[0], extent[1], n_x)[cell % n_x]
    cy = _centres_ref(extent[2], extent[3], grid[0])[cell // n_x]
    expected = np.sqrt(radius ** 2 - (cx ** 2 + cy ** 2))
    assert np.allclose(height, expected, atol=1e-12)
    assert height.max() == pytest.approx(radius, abs=0.1)     # nearest cell centre to (0, 0)


def test_height_candidates_ignore_columns_that_miss_the_sphere():
    from tcren.surface import _height_candidates

    cell, _owner, _height = _height_candidates(np.zeros((1, 3)), np.array([1.0]),
                                               (40, 40), (-10.0, 10.0, -10.0, 10.0))
    assert 0 < len(cell) < 40 * 40                            # a small footprint, not the whole grid


def test_height_candidates_handle_an_empty_structure():
    from tcren.surface import _height_candidates

    cell, owner, height = _height_candidates(np.zeros((0, 3)), np.zeros(0),
                                             (8, 8), (-4.0, 4.0, -4.0, 4.0))
    assert len(cell) == len(owner) == len(height) == 0


def _centres_ref(lo, hi, n):
    edges = np.linspace(lo, hi, n + 1)
    return 0.5 * (edges[:-1] + edges[1:])


# --- the frame -----------------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("pdb_id", ["1ao7", "4ozg"])
def test_groove_frame_is_orthonormal_and_right_handed(annotated, pdb_id):
    _origin, basis = _groove_frame(annotated[pdb_id])
    assert np.allclose(basis @ basis.T, np.eye(3), atol=1e-9)
    assert np.linalg.det(basis) == pytest.approx(1.0, abs=1e-9)


@pytest.mark.slow
@pytest.mark.parametrize("pdb_id", ["1ao7", "4ozg"])
def test_groove_frame_puts_the_helices_on_opposite_sides(annotated, pdb_id):
    """x must be the groove *width*: the two helices straddle it, the peptide runs down the middle.

    This is the regression that matters. Taking x from the floor's principal axis instead put the
    helices diagonally across the map and the peptide off-centre.
    """
    s = annotated[pdb_id]
    origin, basis = _groove_frame(s)

    def local(pts):
        return (np.asarray(pts, float) - origin) @ basis.T

    def helix(region_type):
        return [r.ca for c in s.chains if c.chain_type in ("MHCa", "MHCb", "MHC")
                for reg in c.regions if reg.region_type == region_type
                for r in reg.residues if r.ca is not None]

    helices = [local(pts) for h in ("HELIX_A1", "HELIX_A2", "HELIX_B1") if (pts := helix(h))]
    x_means = sorted(float(h[:, 0].mean()) for h in helices)
    assert len(x_means) == 2, "expected exactly two groove helices"
    assert x_means[0] < -3.0 < 3.0 < x_means[1], f"helices not straddling x=0: {x_means}"

    pep = local([r.ca for c in s.chains if c.chain_type == "PEPTIDE"
                 for r in c.residues if r.ca is not None])
    assert abs(pep[:, 0].mean()) < 1.0                    # centred across the groove
    assert np.ptp(pep[:, 1]) > 3 * np.ptp(pep[:, 0])      # and elongated along it
    assert pep[:, 2].mean() > 0                           # above the floor


@pytest.mark.slow
def test_groove_frame_is_independent_of_input_orientation(annotated):
    """A rigid rotation of the input must not move the map — this is what removes SURFMAP's
    prealignment requirement."""
    import copy

    s = annotated["1ao7"]
    base = surface_map(s)

    theta = 0.9
    rot = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                    [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]])
    shifted = copy.deepcopy(s)
    for chain in shifted.chains:
        for res in chain.residues:
            for atom in res.atoms:
                atom.coord[:] = rot @ atom.coord + np.array([12.0, -5.0, 3.0])
    moved = surface_map(shifted)

    both = np.isfinite(base.channels["h"]) & np.isfinite(moved.channels["h"])
    assert both.sum() > 100
    assert np.abs(base.channels["h"][both] - moved.channels["h"][both]).max() < 0.5


# --- the map -------------------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("pdb_id", ["1ao7", "4ozg"])
def test_surface_map_shape_and_occupancy(annotated, pdb_id):
    m = surface_map(annotated[pdb_id])
    assert isinstance(m, SurfaceMap)
    assert m.extent == DEFAULT_EXTENT
    for name, arr in m.channels.items():
        assert arr.shape == m.grid, name
    assert m.source.shape == m.grid
    assert 0.5 < m.occupancy() <= 1.0
    assert m.peptide and len(m.peptide) >= 8


@pytest.mark.slow
def test_surface_map_cells_carry_their_residue_chemistry(annotated):
    """A cell's hydropathy must be its owning residue's, not an average of the neighbourhood."""
    m = surface_map(annotated["1ao7"], smooth=False)
    phobic = m.channels["phobic"]
    finite = phobic[np.isfinite(phobic)]
    assert finite.size > 100
    assert set(np.round(finite, 3)) <= {round(v, 3) for v in KYTE_DOOLITTLE.values()}


@pytest.mark.slow
def test_smoothing_never_invents_an_unoccupied_cell(annotated):
    rough = surface_map(annotated["1ao7"], smooth=False)
    smooth = surface_map(annotated["1ao7"], smooth=True)
    assert np.array_equal(np.isfinite(rough.channels["h"]), np.isfinite(smooth.channels["h"]))


@pytest.mark.slow
def test_both_helices_and_the_peptide_are_visible(annotated):
    m = surface_map(annotated["1ao7"])
    present = set(np.unique(m.source[np.isfinite(m.channels["h"])]))
    assert SOURCE_CODES["peptide"] in present
    assert {SOURCE_CODES["mhc_helix_a1"], SOURCE_CODES["mhc_helix_a2"]} <= present


@pytest.mark.slow
def test_to_frame_matches_the_occupied_cells(annotated):
    m = surface_map(annotated["4ozg"])
    df = m.to_frame()
    assert df.height == int(np.isfinite(m.channels["h"]).sum())
    assert set(("structure.id", "ix", "iy", "x", "y", "source", "h", "phobic", "charge")) <= set(df.columns)
    assert df["x"].min() >= m.extent[0] and df["x"].max() <= m.extent[1]


# --- scalars and comparison ----------------------------------------------------------------------
@pytest.mark.slow
def test_surface_stats_are_finite_and_bounded(annotated):
    st = surface_stats(surface_map(annotated["1ao7"]))
    assert st["relief"] > 0
    assert st["peak_to_valley"] >= st["relief"]
    assert 0.0 <= st["frac_above_ridge"] <= 1.0
    assert 0.0 < st["area_frac_peptide"] < 1.0
    assert min(KYTE_DOOLITTLE.values()) <= st["phobic_mean"] <= max(KYTE_DOOLITTLE.values())


@pytest.mark.slow
def test_surface_distance_is_a_metric_shape(annotated):
    maps = [surface_map(annotated[p]) for p in ("1ao7", "4ozg")]
    ids, d = surface_distance(maps)
    assert ids == ["1ao7", "4ozg"]
    assert d.shape == (2, 2)
    assert np.allclose(np.diag(d), 0.0)
    assert d[0, 1] == pytest.approx(d[1, 0])
    assert d[0, 1] > 0


@pytest.mark.slow
def test_surface_distance_to_self_is_zero(annotated):
    m = surface_map(annotated["1ao7"])
    _ids, d = surface_distance([m, m])
    assert d[0, 1] == pytest.approx(0.0, abs=1e-12)


def test_surface_distance_rejects_mismatched_grids():
    a = SurfaceMap("a", (4, 4), DEFAULT_EXTENT, {"h": np.zeros((4, 4))}, np.zeros((4, 4), int))
    b = SurfaceMap("b", (8, 8), DEFAULT_EXTENT, {"h": np.zeros((8, 8))}, np.zeros((8, 8), int))
    with pytest.raises(ValueError, match="same grid"):
        surface_distance([a, b])


def test_surface_distance_uses_only_cells_both_maps_reached():
    """A half-empty map must not read as "close to everything" through its missing cells."""
    grid, ext = (2, 2), DEFAULT_EXTENT
    full = SurfaceMap("full", grid, ext, {"h": np.array([[1.0, 1.0], [1.0, 1.0]])}, np.ones(grid, int))
    half = SurfaceMap("half", grid, ext,
                      {"h": np.array([[3.0, np.nan], [np.nan, np.nan]])}, np.ones(grid, int))
    _ids, d = surface_distance([full, half])
    assert d[0, 1] == pytest.approx(2.0)          # |1-3| over the single shared cell, not /4


@pytest.mark.slow
def test_surface_table_one_row_per_map(annotated):
    maps = [surface_map(annotated[p]) for p in ("1ao7", "4ozg")]
    tab = surface_table(maps)
    assert tab.height == 2
    assert {"structure.id", "peptide", "relief", "frac_above_ridge"} <= set(tab.columns)


# --- rendering -----------------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("channel", ["h", "phobic", "charge"])
def test_render_surface_map_is_well_formed_svg(annotated, channel):
    import xml.etree.ElementTree as ET

    from tcren.viz.surface2d import render_surface_map

    svg = render_surface_map(surface_map(annotated["1ao7"], grid=(16, 8)), channel)
    ET.fromstring(svg)                                     # raises if malformed
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert f'data-channel="{channel}"' in svg
    assert "<title>" in svg                                # metadata-bearing, per viz.svg2d


def test_signed_channels_get_a_zero_centred_ramp():
    """An all-hydrophobic surface must not paint its least-greasy cell as hydrophilic."""
    from tcren.viz.surface2d import _normaliser

    values = np.array([2.0, 3.0, 4.0, 4.5])
    to_t, _ramp, lo, hi = _normaliser(values, "phobic")
    assert lo < 0 < hi and lo == pytest.approx(-hi)
    assert to_t(0.0) == pytest.approx(0.5)                 # zero sits at the ramp's midpoint
    assert to_t(2.0) > 0.5                                 # the least greasy cell is still greasy
