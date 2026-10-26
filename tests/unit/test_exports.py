"""Exact exports preserve existing quantities, physical identity and failed inputs."""
import json
import tarfile
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal
from typer.testing import CliRunner

from tcren import exports
from tcren.cli import app
from tcren.contactmap import ContactMap
from tcren.potential import mj, tcren as classic
from tcren.structure import parse_structure
from tcren.structure.model import RegionMarkup

PDB = Path(__file__).resolve().parents[1] / "assets/pdb/1ao7.pdb"


def typed_structure():
    structure = parse_structure(PDB)
    for chain in structure.chains:
        chain.chain_type = {"A": "MHCa", "B": "B2M", "C": "PEPTIDE",
                            "D": "TRA", "E": "TRB"}[chain.chain_id]
        chain.chain_supertype = "MHCI" if chain.chain_type in {"MHCa", "B2M"} else chain.chain_type
        if chain.chain_type == "PEPTIDE":
            chain.regions = [RegionMarkup("PEPTIDE", 0, len(chain.residues)-1, chain.sequence(), chain.residues)]
        if chain.chain_type in {"TRA", "TRB"}:
            residues = chain.residues
            chain.regions = [RegionMarkup("CDR3", 0, len(residues)-1, chain.sequence(), residues)]
    return structure


def test_exact_export_quantities_and_author_identity():
    from tcren.energetics.scoring import position_profile
    from tcren.pipeline import _interface_energy
    from tcren.viz.pymol import residue_importance

    structure = typed_structure()
    pot, other = classic(), mj()
    cm = ContactMap.from_structure(structure)
    before = cm.contacts.clone()
    row = exports.receptor_reference(cm, pot, other)
    for interface, suffix, potential in (("tcr_peptide", "tcr_pep", pot), ("tcr_mhc", "tcr_mhc", other)):
        contacts = cm.interface(interface)
        expected = _interface_energy(contacts, potential) - _interface_energy(
            contacts.with_columns(pl.lit("A").alias("residue.aa.from")), potential)
        assert row[f"dFT_{suffix}_all"] == expected
        assert row[f"dFT_{suffix}_cdr3"] == expected
        assert row[f"dFT_{suffix}_cdr12"] == 0
    for name in exports.RECEPTOR_REGIONS:
        assert row[f"dFT_total_{name}"] == row[f"dFT_tcr_pep_{name}"] + row[f"dFT_tcr_mhc_{name}"]
    assert_frame_equal(cm.contacts, before)
    for quantity, expected in (("position", position_profile(cm, pot, structure=structure)),
                               ("residue", residue_importance(structure, potential=pot))):
        actual = exports.export_structure(structure, quantity, pot, other)
        assert_frame_equal(actual.select(expected.columns).sort(expected.columns),
                           expected.sort(expected.columns))
        assert actual["residue.index.pdb"].null_count() == 0
        assert "insertion_code" in actual.columns
    identity = exports.residue_identity(structure)
    assert identity.height == sum(len(c.residues) for c in structure.chains)
    placement = exports.export_structure(structure, "placement", pot, other)
    assert placement.height == 1
    assert placement["tcr_low"][0] <= placement["tcr_high"][0]


@pytest.mark.parametrize("archived", [False, True])
@pytest.mark.parametrize("quantity", ["contacts", "contact-atoms", "atom-inventory"])
def test_batch_accounts_for_parse_failure_and_zero_contact_rows(tmp_path, monkeypatch, archived, quantity):
    import tcren.annotation
    import tcren.annotation.batch
    import tcren.mhc
    import tcren.structure

    good = tmp_path / "good_full_stem.pdb"
    good.write_bytes(PDB.read_bytes())
    bad = tmp_path / "missing.pdb"
    calls = []
    original_parse = tcren.structure.parse_structure

    def parse(path, pdb_id, **kwargs):
        if pdb_id == "missing":
            raise ValueError("invalid fixture")
        original_parse(path, pdb_id=pdb_id, **kwargs)
        structure = typed_structure()
        if quantity == "contact-atoms":
            for chain in structure.by_type("MHCa"):
                chain.regions = [RegionMarkup("HELIX_A1", 0, len(chain.residues)-1,
                                             chain.sequence(), chain.residues)]
        structure.pdb_id = pdb_id
        return structure

    monkeypatch.setattr(tcren.structure, "parse_structure", parse)
    monkeypatch.setattr(tcren.annotation.batch, "annotate_batch",
                        lambda structures, **kw: calls.append(len(structures)) or [{} for s in structures])
    monkeypatch.setattr(tcren.annotation, "classify_chains", lambda *a, **kw: None)
    monkeypatch.setattr(tcren.mhc, "annotate_mhc_batch", lambda *a, **kw: None)
    monkeypatch.setattr(exports, "annotation_lineage", lambda organisms: {"organisms": list(organisms)})
    out = tmp_path / ("reference.parquet" if quantity in {"contact-atoms", "atom-inventory"} else "reference.tsv")
    sources = [good, bad]
    if archived:
        bad.write_text("broken structure")
        archive = tmp_path / "inputs.tar.gz"
        with tarfile.open(archive, "w:gz") as stream:
            for path in sources:
                stream.add(path, arcname=path.name)
        sources = [archive]
    ledger = exports.export_batch(sources, out, quantity=quantity, peptide_potential=classic(),
                                   mhc_potential=mj(), cutoff=0.00001)
    assert ledger["status"].to_list() == ["ok", "error"]
    if quantity == "atom-inventory":
        assert ledger["n_rows"][0] > 0 and ledger["n_rows"][1] == 0
    else:
        assert ledger["n_rows"].to_list() == [0, 0]
    assert calls == [1]
    assert ledger["pdb.id"].to_list() == ["good_full_stem", "missing"]
    if archived:
        assert ledger["archive_sha256"].null_count() == 0
        assert ledger["source"].str.contains("inputs.tar.gz#").all()
    meta = json.loads(Path(str(out) + ".provenance.json").read_text())
    assert (meta["n_input"], meta["n_ok"], meta["n_error"]) == (2, 1, 1)
    if quantity == "contact-atoms":
        assert pl.read_parquet(out).is_empty()
        assert meta["potentials"] == {}
        assert meta["contact_observations"]["conformer_policy"]
        assert meta["contact_summary_sha256"]
        assert pl.read_csv(str(out) + ".summary.tsv", separator="\t").is_empty()
    elif quantity == "atom-inventory":
        inventory = pl.read_parquet(out)
        assert inventory.height == ledger["n_rows"][0]
        assert meta["potentials"] == {}
        assert meta["atom_inventory"]["conformer_policy"]
        assert meta["atom_inventory"]["noncontact_atoms_retained"]
        assert meta["interface"] == "all"
        assert meta["residue_phi_accounting"].startswith("not applicable")
    else:
        assert meta["potentials"]["peptide"]["sha256"]
    assert Path(str(out) + ".residues.tsv").exists()


@pytest.mark.parametrize("quantity", ["contacts", "contact-atoms", "atom-inventory"])
def test_all_failed_cli_retains_accounting(tmp_path, quantity):
    out = tmp_path / "failed.tsv"
    annotation_path = Path(str(out) + ".annotation.json")
    annotation_path.write_text('{"old": "stale"}')
    result = CliRunner().invoke(app, ["export", "-s", str(tmp_path / "absent.pdb"),
                                    "-o", str(out), "--quantity", quantity,
                                    "--potential", "karnaukhov2022"])
    assert result.exit_code == 1, result.output
    ledger = pl.read_csv(str(out) + ".accounting.tsv", separator="\t")
    assert ledger["status"].to_list() == ["error"]
    assert "FileNotFoundError" in ledger["error"][0]
    assert json.loads(annotation_path.read_text()) == {}


def test_export_rejects_invalid_cutoff():
    with pytest.raises(ValueError, match="cutoff"):
        exports.export_structure(typed_structure(), "contacts", classic(), mj(), cutoff=float("nan"))


def test_contact_atoms_rejects_modified_peptide_before_shortening(tmp_path, monkeypatch):
    import tcren.annotation
    import tcren.annotation.batch
    import tcren.mhc
    lines = PDB.read_text().splitlines(keepends=True)
    peptide_residue = next(line[22:27] for line in lines if line.startswith("ATOM  ") and line[21] == "C")
    modified = tmp_path / "modified.pdb"
    modified.write_text("".join("HETATM" + line[6:17] + "MSE" + line[20:]
        if line.startswith("ATOM  ") and line[21] == "C" and line[22:27] == peptide_residue
        else line for line in lines))
    assert "X" in parse_structure(modified, keep_modified=True).chain("C").sequence()
    assert "X" not in parse_structure(modified).chain("C").sequence()
    def classify(s, **kwargs):
        for chain in s.chains:
            chain.chain_type = {"A": "MHCa", "B": "B2M", "C": "PEPTIDE", "D": "TRA", "E": "TRB"}[chain.chain_id]
            chain.chain_supertype = "MHCI" if chain.chain_type in {"MHCa", "B2M"} else chain.chain_type
            region = "HELIX_A1" if chain.chain_type == "MHCa" else "CDR3"
            chain.regions = [RegionMarkup(region, 0, len(chain.residues)-1, chain.sequence(), chain.residues)]
    monkeypatch.setattr(tcren.annotation.batch, "annotate_batch", lambda structures, **kw: [{} for _ in structures])
    monkeypatch.setattr(tcren.annotation, "classify_chains", classify)
    monkeypatch.setattr(tcren.mhc, "annotate_mhc_batch", lambda *a, **kw: None)
    monkeypatch.setattr(exports, "annotation_lineage", lambda *a: {})
    ledger = exports.export_batch([modified], tmp_path / "atoms.parquet", quantity="contact-atoms",
                                 peptide_potential=None, mhc_potential=None, interface="all")
    assert ledger["status"].to_list() == ["error"]
    assert "standard-amino-acid peptide" in ledger["error"][0]


def test_structural_reference_rebuilds_contacts_and_preserves_input():
    """Directly truncated coordinates are the oracle; the virtual map must differ."""
    from tcren.pipeline import _interface_energy, run, score_row
    from tcren.refine.substitute import substitute_peptide, substitute_residues

    s = typed_structure()
    s.mhc_calls = []
    cm = ContactMap.from_structure(s)
    before = cm.contacts.clone()
    pot, other = classic(), mj()
    peptide = next(c.sequence() for c in s.chains if c.chain_type == 'PEPTIDE')
    ref = substitute_peptide(s, 'A' * len(peptide))
    rc = ContactMap.from_structure(ref)
    assert rc.interface('tcr_peptide').height < cm.interface('tcr_peptide').height
    structural = {}
    for weight in ('residue', 'atomic'):
        actual = run(s, typed=True, superimpose=False, reference_aa='A',
                     reference_mode='structural', contact_weight=weight,
                     potentials={'tcr_peptide': 'karnaukhov2022'})
        real_map = ContactMap.from_structure(s, count_atoms=weight == 'atomic')
        ref_map = ContactMap.from_structure(ref, count_atoms=weight == 'atomic')
        expected = (_interface_energy(real_map.interface('tcr_peptide'), pot, contact_weight=weight)
                    - _interface_energy(ref_map.interface('tcr_peptide'), pot, contact_weight=weight))
        structural[weight] = expected
        assert actual.scores['delta_tcr_peptide'] == pytest.approx(expected)
        assert score_row(actual)['reference_mode'] == 'structural'
    virtual = run(s, typed=True, superimpose=False, reference_aa='A',
                  potentials={'tcr_peptide': 'karnaukhov2022'})
    assert virtual.scores['delta_tcr_peptide'] != pytest.approx(structural['residue'])
    targets = {(c.chain_id, r.seq_index): 'A' for c in s.chains
               if c.chain_type in ('TRA', 'TRB') for r in c.residues}
    receptor_map = ContactMap.from_structure(substitute_residues(s, targets))
    actual = exports.receptor_reference(cm, pot, other, structure=s)
    for interface, suffix, potential in [('tcr_peptide', 'tcr_pep', pot), ('tcr_mhc', 'tcr_mhc', other)]:
        expected = (_interface_energy(cm.interface(interface), potential)
                    - _interface_energy(receptor_map.interface(interface), potential))
        assert actual[f'dFT_{suffix}_all'] == pytest.approx(expected)
    assert actual['dFT_tcr_pep_cdr12'] == pytest.approx(0)
    assert actual['dFT_tcr_pep_cdr3'] == pytest.approx(actual['dFT_tcr_pep_all'])
    assert_frame_equal(ContactMap.from_structure(s).contacts, before)


def test_structural_reference_cli_requires_alanine_delta(tmp_path):
    result = CliRunner().invoke(app, ['scoring', '-s', str(PDB), '-o',
                                     str(tmp_path / 'scores.csv'),
                                     '--reference-mode', 'structural'])
    assert result.exit_code == 2
    assert 'requires' in result.output
    assert not (tmp_path / 'scores.csv').exists()


def test_separate_contact_topology():
    # K2,2 has one component and one graph cycle; disjoint edges have two components.
    structure = typed_structure()
    cm = ContactMap.from_structure(structure)
    from types import SimpleNamespace
    schema = {"chain.id.from": pl.String, "residue.index.from": pl.Int64,
              "chain.id.to": pl.String, "residue.index.to": pl.Int64,
              "region.type.from": pl.String}
    cycle = pl.DataFrame([("D", i, "C", j, "CDR3") for i in range(2) for j in range(2)],
                         schema=schema, orient="row")
    separate = pl.DataFrame([("C", i, "A", i, "PEPTIDE") for i in range(2)],
                            schema=schema, orient="row")
    frames = {"tcr_peptide": cycle, "tcr_mhc": cycle.clear(), "peptide_mhc": separate}
    rows = exports.contact_topology(SimpleNamespace(interface=frames.__getitem__)).to_dicts()
    assert [(r["vertices"], r["edges"], r["beta0"], r["beta1"], r["available"])
            for r in rows] == [(4,4,1,1,True), (0,0,0,0,False), (4,2,2,0,True)]
    actual = exports.export_structure(structure, "contact-topology", None, None)
    assert actual["interface"].to_list() == list(frames)
    assert_frame_equal(exports.contact_topology(cm), actual)


def test_backbone_export_reuses_geometry_and_is_rigid_invariant():
    import numpy as np
    from tcren.geometry import cdr3_internal_coords
    s = typed_structure()
    row = exports.export_structure(s, "backbone", None, None)
    coords = cdr3_internal_coords(s, 'cdr3a')
    bonds = row.filter((pl.col('chain.type') == 'TRA') & (pl.col('metric') == 'ca_bond'))
    np.testing.assert_allclose(bonds['value'], coords.bonds)
    for c in s.chains:
        for r in c.residues:
            for a in r.atoms:
                a.coord[:] = a.coord @ np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]) + np.array([4., -3., 2.])
    moved = exports.backbone_table(s)
    np.testing.assert_allclose(row['value'], moved['value'], atol=1e-10)
    assert 'residue.index.pdb.1' in row.columns


def test_contact_context_graph_oracle_and_exact_energy():
    import numpy as np
    from tcren.pipeline import _interface_energy
    structure = typed_structure()
    potential = classic()
    table = exports.export_structure(structure, 'contact-context', potential, mj())
    keys = [(r['chain.id.from'], r['residue.index.from'], r['chain.id.to'], r['residue.index.to'])
            for r in table.to_dicts()]
    edges = set(keys)
    for key, count in zip(keys, table['four_cycles']):
        a, i, b, j = key
        expected = sum(x != (a,i) and y != (b,j) and (*x,*y) in edges
                       for x in {(k[0],k[1]) for k in keys if k[2:]==(b,j)}
                       for y in {(k[2],k[3]) for k in keys if k[:2]==(a,i)})
        assert count == expected
    assert table['pair_energy'].sum() == pytest.approx(
        _interface_energy(ContactMap.from_structure(structure).interface('tcr_peptide'), potential))
    before = table.select('surface_h_A','surface_hydropathy','ridge_height_A')
    for c in structure.by_type('TRA','TRB'):
        for r in c.residues:
            for atom in r.atoms:
                atom.coord[:] += np.array([0.,0.,.01])
    after = exports.export_structure(structure, 'contact-context', potential, mj())
    np.testing.assert_allclose(before.row(0), after.select(before.columns).row(0), equal_nan=True)


def test_contact_summary_keeps_zero_models_and_composition_ties():
    from tcren.topology.contact_context import summarize_contacts
    table = exports.export_structure(typed_structure(), 'contact-context', classic(), mj())
    # Same composition in opposite contact order; retain a separate zero-contact model.
    table = pl.concat([table.with_columns(pl.lit('a').alias('pdb.id')),
                       table.reverse().with_columns(pl.lit('b').alias('pdb.id'))])
    n = table.height//2
    ledger = pl.DataFrame({'pdb.id':['a','b','zero'],'status':['ok']*3,'n_rows':[n,n,0]})
    result = summarize_contacts(table, ledger).sort('pdb.id')
    assert result['context_energy'][0] == result['context_energy'][1]
    assert result['cycle_energy'][0] == result['cycle_energy'][1]
    assert result['n_all'][2] == result['context_energy'][2] == 0
    assert result['n_missing_topography'][0] > 0  # Fixture has no MHC helix markup.
    assert result['bulge_energy'].is_nan()[0]
    with pytest.raises(ValueError, match='accounting'):
        summarize_contacts(table, ledger.with_columns(pl.lit(1).alias('n_rows')))



def test_contact_summary_accepts_header_only_csv():
    from tcren.topology.contact_context import summarize_contacts
    table = exports.export_structure(typed_structure(), 'contact-context', classic(), mj())
    empty = pl.read_csv(table.head(0).write_csv(separator="\t").encode(), separator="\t")
    ledger = pl.DataFrame({'pdb.id': ['a', 'b'], 'status': ['ok', 'ok'], 'n_rows': [0, 0]})
    result = summarize_contacts(empty, ledger)
    assert result['pdb.id'].to_list() == ['a', 'b']
    assert all(value == 0 for row in result.drop('pdb.id').rows() for value in row)
    with pytest.raises(ValueError, match='accounting'):
        summarize_contacts(empty, ledger.with_columns(pl.lit(1).alias('n_rows')))


def test_recognition_graph_ladder():
    from types import SimpleNamespace
    import numpy as np
    schema = {"chain.id.from": pl.String, "residue.index.from": pl.Int64,
              "chain.id.to": pl.String, "residue.index.to": pl.Int64,
              "region.type.from": pl.String}
    def frame(rows):
        return pl.DataFrame(rows, schema=schema, orient="row")
    frames = {"tcr_peptide": frame([("T",0,"P",0,"CDR3"), ("T",0,"P",1,"CDR3")]),
              "peptide_mhc": frame([("P",0,"M",0,"PEPTIDE"), ("P",1,"M",0,"PEPTIDE")]),
              "tcr_mhc": frame([("T",0,"M",0,"CDR3"), ("F",0,"M",0,"FR1")])}
    residues = [SimpleNamespace(seq_index=i, aa="A", atoms=[SimpleNamespace(name=n, element=n[0], coord=np.array(x))])
                for i,n,x in [(0,"C",[0.,0,0]), (1,"N",[1.33,0,0])]]
    structure = SimpleNamespace(by_type=lambda _: [SimpleNamespace(chain_id="P",residues=residues)])
    cm = SimpleNamespace(interface=frames.__getitem__)
    result = exports.recognition_topology(structure, cm)
    assert result["beta1"].to_list() == [0,1,2,3,1,2,2]
    assert result["edges"].to_list() == [3,4,5,6,4,5,5]
    assert result["sequence_breaks"].to_list() == [0]*7
    residues[1].atoms[0].coord[0] = 10
    assert exports.recognition_topology(structure, cm)["beta1"].to_list() == [0,1,1,2,1,2,2]
    # Spatial side-chain edges survive a broken backbone; duplicate atom pairs
    # count once. Hydrogen and backbone proximity cannot create these edges.
    for residue, x in zip(residues, [0., 5.]):
        residue.atoms.extend([SimpleNamespace(name="CB", element="C", coord=np.array([x, 3., 0.])),
                              SimpleNamespace(name="CG", element="C", coord=np.array([x, 3., 0.])),
                              SimpleNamespace(name="HB", element="H", coord=np.array([0., 3., 0.]))])
    spatial = exports.recognition_topology(structure, cm)
    assert spatial["beta1"].to_list() == [0,1,1,2,2,3,2]
    assert spatial["peptide_sidechain_edges"].to_list() == [1]*7
    assert exports.recognition_topology(structure, cm, cutoff=4.9)["peptide_sidechain_edges"].to_list() == [0]*7
    from tcren.topology.graph import recognition_graph_features
    production = recognition_graph_features(cm)
    assert production["g_recognition_cyclo_frac"] == 2/5
    assert production["g_recognition_vertices"] == 4
    tripartite = spatial.filter(pl.col("graph") == "tp_pm_tm")
    assert tripartite.select("edges", "simplicial_beta1", "simplicial_beta2", "triangles", "tetrahedra").row(0) == (5, 0, 0, 2, 0)
    frames["tcr_peptide"] = frames["tcr_peptide"].clear()
    frames["tcr_mhc"] = frames["tcr_mhc"].clear()
    assert not exports.recognition_topology(structure, cm)["available"].any()
    frames["peptide_mhc"] = frame([("P",i,"M",j,"PEPTIDE") for i in [0,1] for j in [0,1]])
    raw = exports.recognition_betti(cm)
    assert raw["complex_cyclo_frac"].item() == 1/4
    assert not raw["receptor_contacts_available"].item()
    assert raw["g_recognition_cyclo_frac"].is_nan().item()
    real = typed_structure()
    actual = exports.export_structure(real, "recognition-topology", None, None)
    compact = exports.export_structure(real, "recognition-betti", None, None)
    expected = actual.filter(pl.col("graph") == "tp_pm_tm")
    assert compact["g_recognition_edges"].item() == expected["edges"].item()
    assert compact["g_recognition_beta1"].item() == expected["beta1"].item()
    assert compact["g_recognition_cyclo_frac"].item() == expected["cyclo_frac"].item()
    assert compact["complex_cyclo_frac"].item() == expected["cyclo_frac"].item()
    from tcren.topology.graph import graph_features
    assert actual.filter(pl.col("graph") == "released")["cyclo_frac"].item() == graph_features(real)["g_cyclo_frac"]


def test_simplicial_homology_fills_triangles_and_tetrahedra():
    from itertools import combinations
    from tcren.topology.graph import _clique_homology

    def betti(edges, cycle_rank):
        result = _clique_homology(set(edges), cycle_rank)
        return result['simplicial_beta1'], result['simplicial_beta2']

    assert betti([], 0) == (0, 0)
    assert betti([(0, 1)], 0) == (0, 0)
    assert betti([(0, 1), (1, 2), (2, 3), (0, 3)], 1) == (1, 0)
    assert betti(combinations(range(3), 2), 1) == (0, 0)
    # The octahedral surface has one H2 generator despite having no tetrahedra.
    shell = {(a, b) for a, b in combinations(range(6), 2) if a//2 != b//2}
    result = _clique_homology(shell, 7)
    assert (result['triangles'], result['tetrahedra']) == (8, 0)
    assert betti(shell, 7) == (0, 1)
    # A filled tetrahedron has four triangles, but no H2 hole. K5 also has none:
    # the 3-skeleton suffices even though its clique complex has a 4-simplex.
    assert betti(combinations(range(4), 2), 3) == (0, 0)
    assert betti(combinations(range(5), 2), 6) == (0, 0)
    assert betti(shell | {(0, 1)}, 8) == (0, 0)
    with pytest.raises(ValueError, match='sorted endpoints'):
        betti([(1, 0)], 0)


@pytest.mark.parametrize('quantity', ['weighted-recognition-topology', 'persistence'])
def test_weighted_and_persistence_exports(quantity):
    if quantity == 'persistence':
        pytest.importorskip('gudhi')
    import numpy as np
    from tcren.topology.graph import _recognition_contact_edges
    structure = typed_structure()
    assert quantity in exports.QUANTITIES
    result = exports.export_structure(structure, quantity, classic(), mj())
    if quantity == 'weighted-recognition-topology':
        reference = exports.recognition_topology(structure, ContactMap.from_structure(structure))
        assert_frame_equal(result.select('graph', 'vertices', 'edges', 'beta0', 'beta1'),
                           reference.select('graph', 'vertices', 'edges', 'beta0', 'beta1'))
        assert result.filter(pl.col('graph').is_in(['released', 'tp_pm', 'tp_pm_tm']))['missing_weight'].sum() == 0
        assert np.isfinite(result.filter(pl.col('graph') == 'tp_pm_tm')['w_total'][0])
    else:
        parts = _recognition_contact_edges(ContactMap.from_structure(structure))
        vertices = {v for edges in parts.values() for edge in edges for v in edge}
        assert result['vertices'][0] == len(vertices)
        assert result['points'][0]+result['missing_atoms'][0] == len(vertices)
        assert 'pers1_bars_cb' in result.columns


@pytest.mark.parametrize('quantity', ['weighted-recognition-topology', 'persistence', 'cdr3-contacts'])
def test_new_topology_cli_requires_no_user_potential(monkeypatch, tmp_path, quantity):
    captured = {}
    def batch(*args, **kwargs):
        captured.update(kwargs)
        return pl.DataFrame({'status': ['ok']})
    monkeypatch.setattr(exports, 'export_batch', batch)
    result = CliRunner().invoke(app, ['export', '-s', str(PDB), '-o', str(tmp_path/'out.tsv'),
                                    '--quantity', quantity])
    assert result.exit_code == 0, result.output
    assert captured['quantity'] == quantity
    assert captured['peptide_potential'] is None


def test_cdr3_neighborhood_export_uses_complete_interfaces():
    structure = typed_structure()
    result = exports.export_structure(structure, 'cdr3-neighborhood', None, None)
    cm = ContactMap.from_structure(structure)
    assert result.height == 1
    assert result['peptide_edges'][0] == cm.interface('tcr_peptide').height
    assert result['mhc_edges'][0] == cm.interface('tcr_mhc').height


def test_cdr3_contacts_export_retains_all_atom_chemical_witnesses():
    structure = typed_structure()
    atoms = exports.export_structure(structure, 'cdr3-contacts', None, None)
    from tcren.topology.neighborhood import neighborhood_counts, neighborhood_features
    from tcren.contact_types import stacked_pairs
    assert atoms['region.type.from'].unique().to_list() == ['CDR3']
    assert set(atoms['interface']) == {'peptide', 'mhc'}
    assert 'is_hydrophobic' in atoms.columns and 'sidechain.from' in atoms.columns
    expected = neighborhood_features(structure)
    expected.pop('pdb.id')
    assert neighborhood_counts(atoms, stacking=stacked_pairs(structure)) == expected


@pytest.mark.parametrize("interface", ["tcr_mhc", "peptide_mhc"])
def test_mhc_derivation_contact_table_preserves_physical_edges(interface):
    from tcren.paper.helpers import contact_table
    from copy import deepcopy
    structure = typed_structure()
    # Coincident class-II alpha/beta sites must count as separate physical edges.
    duplicate = deepcopy(structure.by_type('MHCa')[0])
    duplicate.chain_id = 'Z'
    duplicate.chain_type = 'MHCb'
    structure.chains.append(duplicate)
    result = contact_table(structure, interface=interface)
    expected = ContactMap.from_structure(structure).interface(interface)
    keys = ['chain.id.from','residue.index.from','chain.id.to','residue.index.to']
    assert result.height == expected.height
    assert_frame_equal(result.select(keys).sort(keys),expected.select(keys).sort(keys))
    assert set(result['chain.type.to']) >= {'MHCa','MHCb'}


def test_weighted_cli_forwards_all_three_explicit_matrices(monkeypatch, tmp_path):
    captured = {}
    def batch(*args, **kwargs):
        captured.update(kwargs)
        return pl.DataFrame({'status':['ok']})
    monkeypatch.setattr(exports, 'export_batch', batch)
    paths=[]
    for label,value in [('tp',-1.),('tm',-2.),('pm',-3.)]:
        path=tmp_path/(label+'.csv')
        pl.DataFrame({'residue.aa.from':['A'],'residue.aa.to':['D'],'value':[value]}).write_csv(path)
        paths.append(str(path))
    result=CliRunner().invoke(app,['export','-s',str(PDB),'-o',str(tmp_path/'out.tsv'),
        '--quantity','weighted-recognition-topology','--potential',paths[0],
        '--mhc-potential',paths[1],'--pmhc-potential',paths[2],
        '--graph-weighting','boltzmann','--graph-beta','1.0'])
    assert result.exit_code==0,result.output
    assert captured['graph_weighting']=='boltzmann' and captured['graph_beta']==1.0
    for field,value in [('peptide_potential',-1.),('mhc_potential',-2.),('pmhc_potential',-3.)]:
        assert captured[field].matrix['value'].to_list()==[value]
