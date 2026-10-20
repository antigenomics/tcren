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
