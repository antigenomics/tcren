"""Exact coordinate correspondence must not guess historical residue numbering."""
import numpy as np
import polars as pl


def test_exact_mapping_rejects_ambiguous_atoms_and_recovers_frames():
    from tcren.contacts.prior_geometry import endpoint_frames
    raw = pl.DataFrame({
        "resname": ["ALA"] * 4, "atom": ["N", "CA", "CB", "CB"],
        "x": [0., 1., 2., 2.], "y": [0.] * 4, "z": [0.] * 4,
        "chain": ["A", "A", "A", "B"], "resid": ["8"] * 4,
        "icode": [""] * 4, "altloc": [""] * 4,
    }).with_columns(pl.col("x", "y", "z").cast(pl.Float32))
    endpoints = pl.DataFrame({"resname": ["ALA", "ALA", "ALA"],
        "atom": ["CA", "CB", "CA"], "x": [1., 2., 1.001],
        "y": [0.] * 3, "z": [0.] * 3,
    }).with_columns(pl.col("x", "y", "z").cast(pl.Float32))
    result = endpoint_frames(endpoints, raw)
    assert result["matched"].to_list() == [True, False, False]
    assert result["n.x"][0] == 0 and result["cb.x"][0] == 2
    assert np.isnan(result["cb.x"].to_numpy()[1:]).all()


def test_external_writer_produces_readable_gzip(tmp_path, monkeypatch):
    import gzip
    from tcren.contacts import prior_geometry
    table = {"PDB": ["test"], "Contact distance": [3.]}
    for side, chain, aa, x in (("Ab", "A", "ALA", 0.), ("Ag", "B", "GLY", 3.)):
        table.update({side+" chain ID": [chain], side+" resi ID": ["1"], side+" atom": ["CA"],
                      side+" resi aa": [aa], side+" atom coord": [[x, 0., 0.]]})
    pl.DataFrame(table).write_parquet(tmp_path/"Output_contact.parquet")
    raw = pl.DataFrame(dict(resname=["ALA", "GLY"], atom=["CA", "CA"], x=[0., 3.],
        y=[0., 0.], z=[0., 0.], chain=["A", "B"], resid=["1", "1"], icode=["", ""], altloc=["", ""]))
    raw = raw.with_columns(pl.col("x", "y", "z").cast(pl.Float32))
    monkeypatch.setattr(prior_geometry, "raw_atoms", lambda _: raw)
    prior_geometry.external_geometry(tmp_path, tmp_path/"out", set())
    with gzip.open(tmp_path/"out/external_geometry.tsv.gz", "rt") as handle:
        text = handle.read()
    assert text.startswith("pdb.id\t") and "test" in text
