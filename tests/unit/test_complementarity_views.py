import numpy as np
import pytest
import gzip
import hashlib
import json
import polars as pl
from polars.testing import assert_frame_equal

from tcren.complementarity.frame import CanonicalFrame
from tcren.complementarity.scene import build_scene
from tcren.complementarity.views import LOOPS, heatmap_table, topography_tables, validate_grid
from tcren.structure.model import Atom, Chain, Residue, Structure


def example(task="tcr", *, missing_ca=False):
    chains = []
    for name, role, xyz in [("A", "TRA", [0., 0., 5.]), ("P", "PEPTIDE", [3., 4., 5.])]:
        atom_name = "N" if missing_ca and name == "A" else "CA"
        atom = Atom(atom_name, "C", np.array(xyz))
        chains.append(Chain(name, [Residue(0, 10, "", "A", "ALA", (atom,))], role))
    scene = build_scene(Structure("native", chains), CanonicalFrame(np.eye(3), np.zeros(3), "reference", "apex"))
    scene["task"] = task
    scene["residues"][0]["region"] = "CDR3"
    return scene


def test_topography_uses_scene_frame_and_existing_kd_values():
    scene = example()
    residues, cells = topography_tables(scene, grid=(1, 1), extent=(-.5, .5, -.5, .5))
    assert residues["ca_z"].to_list() == [5.]
    assert residues["face_apex_z"].to_list() == [5.]
    assert cells["surface_z"][0] == pytest.approx(5.-1.7-1.4)
    assert cells["owner_atom"][0] == 0
    assert cells["hydropathy"][0] == 1.8
    assert cells["reference"][0] == "reference"
    _, peptide = topography_tables(example("epitope"), grid=(1, 1), extent=(2.5, 3.5, 3.5, 4.5))
    assert peptide["surface_z"][0] == pytest.approx(5.+1.7+1.4)
    assert peptide["owner_atom"][0] == 1


def test_empty_cells_and_unknown_hydropathy_stay_null():
    scene = example()
    scene["residues"][0]["aa"] = "X"
    _, occupied = topography_tables(scene, grid=(1, 1), extent=(-.5, .5, -.5, .5))
    assert occupied["hydropathy"][0] is None
    _, empty = topography_tables(scene, grid=(2, 3), extent=(100., 110., 100., 110.))
    assert empty.height == 6
    assert empty["occupied"].sum() == 0
    assert empty["surface_z"].null_count() == 6
    assert empty["owner_residue"].null_count() == 6


def test_heatmap_known_distance_multiplicity_and_full_facet_order():
    scene = example()
    table, facets = heatmap_table(scene)
    assert table.height == 1
    assert table["ca_distance"][0] == 5.
    assert table["atom_contact_count"][0] == 1
    assert table["contacting"][0]
    assert table["loop_facet"][0] == "alpha:CDR3"
    assert len(facets) == 18
    assert [f["loop_facet"] for f in facets[::3]] == list(LOOPS)
    assert facets[-1]["partner_residues"] == 0


def test_contacting_pair_without_ca_has_null_distance_not_zero():
    table, _ = heatmap_table(example(missing_ca=True))
    assert table["ca_distance"][0] is None
    assert not table["has_both_ca"][0]
    assert table["atom_contact_count"][0] == 1
    assert table["contacting"][0]


@pytest.mark.parametrize("grid,extent,probe", [
    ((0, 1), (-1, 1, -1, 1), 1.4), ((True, 1), (-1, 1, -1, 1), 1.4),
    ((257, 1), (-1, 1, -1, 1), 1.4), ((1, 1), (1, -1, -1, 1), 1.4),
    ((1, 1), (-1, 1, -1, float("nan")), 1.4), ((1, 1), (-1, 1, -1, 1), -1),
    ((1, 1), (-1e308, 1e308, -1, 1), 1.4),
])
def test_invalid_grid_bounds_and_probe(grid, extent, probe):
    with pytest.raises(ValueError):
        validate_grid(grid, extent, probe)


def test_heatmap_rejects_nonfinite_resolved_ca():
    scene = example()
    scene["atoms"][1]["xyz"][0] = float("inf")
    with pytest.raises(ValueError, match="finite"):
        heatmap_table(scene)


def test_class_ii_second_helix_and_explicit_missing_facets():
    scene = example()
    scene["residues"][1].update(role="mhc", region="HELIX_B1")
    table, facets = heatmap_table(scene)
    assert table["partner_facet"].to_list() == ["helix_2"]
    assert table["partner_region"].to_list() == ["HELIX_B1"]
    assert next(f for f in facets if f["loop_facet"] == "alpha:CDR3" and f["partner_facet"] == "peptide")["partner_residues"] == 0
    scene["residues"].append({**scene["residues"][1], "region": "HELIX_A2"})
    with pytest.raises(ValueError, match="ambiguous"):
        heatmap_table(scene)


def test_export_owns_formats_cache_summary_and_frames(tmp_path):
    from tcren.complementarity.dashboard import write_dashboard
    from tcren.complementarity.view_export import export_views
    from tcren.complementarity.view_plot import read_structure_tables

    source, output = tmp_path / "source", tmp_path / "views"
    scene = example()
    write_dashboard([scene], source)
    original = (source / "scene-0.json.gz").read_bytes()
    records = export_views(source, output, grid=(1, 1), extent=(-.5, .5, -.5, .5))
    assert original == (source / "scene-0.json.gz").read_bytes()
    cache = json.loads(gzip.decompress((output / records[0]["cache"]).read_bytes()))
    tables, _, manifest = read_structure_tables(output, "native")
    for file in manifest["files"]:
        name = file["table"]
        table = tables[name]
        restored = pl.read_csv(output / file["tsv_file"], separator="\t", schema=table.schema)
        assert_frame_equal(table, restored)
        assert table.to_dicts() == cache["tables"][name]
        assert cache["source_sha256"][name] == hashlib.sha256((output / file["file"]).read_bytes()).hexdigest()
    summary = pl.read_parquet(output / "summary.parquet").row(0, named=True)
    assert summary["occupied_cells"] == 1
    assert summary["occupied_boundary_cells"] == 1
    assert summary["surface_may_be_clipped"]
    assert summary["atom_contact_witnesses"] == 1
    assert summary["surface_z_min"] == pytest.approx(1.9)
    with pytest.raises(FileExistsError):
        export_views(source, output)
    with pytest.raises(ValueError, match="unknown structure"):
        read_structure_tables(output, "missing")


def test_unknown_chemistry_summary_remains_typed(tmp_path):
    from tcren.complementarity.dashboard import write_dashboard
    from tcren.complementarity.view_export import export_views

    scene = example()
    scene["residues"][0]["aa"] = "X"
    source, output = tmp_path / "source", tmp_path / "views"
    write_dashboard([scene], source)
    export_views(source, output, grid=(3, 3), extent=(-10., 10., -10., 10.))
    summary = pl.read_parquet(output / "summary.parquet")
    assert summary.schema["occupied_cell_mean_hydropathy"] == pl.Float64
    assert summary["occupied_cell_mean_hydropathy"][0] is None
    assert not summary["surface_may_be_clipped"][0]


@pytest.mark.parametrize("identity", ["", "  ", None, 3])
def test_export_refuses_invalid_source_id_before_writing(tmp_path, identity):
    from tcren.complementarity.view_export import export_views

    source, output = tmp_path / "source", tmp_path / "views"
    source.mkdir()
    (source / "index.json").write_text(json.dumps([{"id": identity, "file": "absent.json.gz"}]))
    with pytest.raises(ValueError, match="unique nonempty"):
        export_views(source, output)
    assert not output.exists()


def test_views_cli_exports_and_rejects_invalid_grid_before_writing(tmp_path):
    import typer
    from typer.testing import CliRunner
    from tcren.complementarity.dashboard import write_dashboard
    from tcren.complementarity.view_export import views_command

    source = tmp_path / "source"
    write_dashboard([example()], source)
    app = typer.Typer()
    app.command()(views_command)
    runner = CliRunner()
    for arguments in [["--grid", "0"], ["--grid", "257"], ["--probe", "nan"],
                      ["--extent", "1", "-1", "-1", "1"]]:
        output = tmp_path / "invalid"
        result = runner.invoke(app, ["--exports", str(source), "-o", str(output), *arguments])
        assert result.exit_code != 0
        assert not output.exists()
    output = tmp_path / "valid"
    result = runner.invoke(app, ["--exports", str(source), "-o", str(output),
                                 "--grid", "3", "--probe", "0", "--extent", "-10", "10", "-10", "10"])
    assert result.exit_code == 0, result.output
    assert "1 structures" in result.output
    assert pl.read_parquet(output / "summary.parquet")["grid_cells"][0] == 9
    again = runner.invoke(app, ["--exports", str(source), "-o", str(output)])
    assert again.exit_code != 0


def test_static_plot_explains_zero_contact_facets_and_uses_navigation_label(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    from matplotlib.figure import Figure
    from tcren.complementarity.dashboard import write_dashboard
    from tcren.complementarity.view_export import export_views
    from tcren.complementarity.view_plot import plot_views

    source, output = tmp_path / "source", tmp_path / "views"
    scene = example()
    scene["residue_contacts"] = []
    write_dashboard([scene], source)
    index = json.loads((source / "index.json").read_text())
    index[0]["display_label"] = "Known peptide · native"
    (source / "index.json").write_text(json.dumps(index))
    export_views(source, output, grid=(3, 3))
    captions = []
    original = Figure.savefig

    def capture(figure, *args, **kwargs):
        captions.append((figure._suptitle.get_text(),
                         [text.get_text() for axis in figure.axes for text in axis.texts]))
        return original(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", capture)
    paths = plot_views(output, "native", tmp_path / "plots", contact_mode="pairs")
    assert all(path.is_file() for path in paths)
    assert all("Known peptide · native" in title for title, _ in captions)
    assert "No contacting pairs" in captions[1][1]


@pytest.mark.parametrize("extent,probe", [
    ((-1., 1., -1., 1.), 1e100),
    ((0., 1e-300, 0., 1e-300), 1.4),
    ((0., 1e-320, 0., 1e-320), 1.4),
])
def test_topography_refuses_unbounded_offset_work_before_sampler(monkeypatch, extent, probe):
    def forbidden(*args, **kwargs):
        pytest.fail("unbounded request reached ray caster")

    monkeypatch.setattr("tcren.complementarity.views._height_candidates", forbidden)
    with pytest.raises(ValueError, match="atom/offset evaluations"):
        topography_tables(example(), grid=(256, 256), extent=extent, probe=probe)
