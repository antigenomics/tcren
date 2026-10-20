"""Coarse workers preserve every ordered score, surface row and input failure."""
from dataclasses import replace

import polars as pl
from polars.testing import assert_frame_equal
import pytest

from tcren.contacts.observations import contact_observations
from tcren.contacts.penalty_export import _score_chunk
from tcren.contacts.plausibility import fit_resource
from tcren.pipeline import _BUNDLED_POTENTIALS
from tcren.runtime import available_cpus, native_threads, parallel_map
from test_contact_observations import complex_fixture


@pytest.mark.parametrize("mode", ["scores", "geometry", "deepest_sphere", "accessible_boundary", "atom_count_and_chemical_core", "atom_pattern"])
def test_contact_penalty_coarse_workers_equal_serial(mode):
    if available_cpus() < 2:
        pytest.skip("requires two allocated CPUs")
    structures = [replace(complex_fixture(), pdb_id=f"synthetic-{i}") for i in range(3)]
    observations = contact_observations(structures[0], interface="tcr_peptide")
    reference = pl.concat([observations, observations.with_columns((pl.col("dist")*.95).alias("dist"))])
    resource = fit_resource(reference, observations)
    geometry = mode == "geometry"
    tables = {} if geometry else {name: factory() for name, factory in _BUNDLED_POTENTIALS.items()
                                 if name in ("tcren2", "karnaukhov2022")}
    atom_count = chemical_core = None
    if mode in ("atom_count_and_chemical_core", "atom_pattern"):
        from tcren.contacts.atom_counts import TYPES, derive_atom_counts
        from tcren.contacts.chemical_core import derive_chemical_core
        swap = pl.col("role.from") != "tcr"
        counts = observations.select(*[
            pl.when(swap).then(pl.col("chem."+other)).otherwise(pl.col("chem."+side))
            .str.split_exact(":", 1).struct.field(field).alias(dest+"_"+side)
            for side, other in [("from", "to"), ("to", "from")]
            for field, dest in [("field_0", "aa"), ("field_1", "atom")]
        ]).group_by(TYPES).len(name="count")
        if mode == 'atom_pattern':
            from tcren.contacts.atom_patterns import derive_atom_patterns
            data = observations.select(
                pl.struct('chain.id.from','residue.index.from','chain.id.to','residue.index.to').alias('edge'), *[
                    pl.when(swap).then(pl.col('chem.'+other)).otherwise(pl.col('chem.'+side))
                    .str.split_exact(':',1).struct.field(field).alias(dest+'_'+side)
                    for side,other in [('from','to'),('to','from')]
                    for field,dest in [('field_0','aa'),('field_1','atom')]])
            atom_count = derive_atom_patterns(data,data.head(0))
        else:
            chemical_core = derive_chemical_core(counts, counts.head(0))
            atom_count = derive_atom_counts(counts, counts.head(0))
    options = (tables, resource, geometry, None, None,
               mode in ("deepest_sphere", "accessible_boundary"),
               mode if mode in ("deepest_sphere", "accessible_boundary") else "deepest_sphere", atom_count, chemical_core)
    # Unequal contiguous chunks include an ordinary per-structure failure.
    work = [(structures[:2], *options), ([None, structures[2]], *options)]
    with native_threads():
        serial = parallel_map(_score_chunk, work, 1)
        parallel = parallel_map(_score_chunk, work, 2)
    for expected_chunk, actual_chunk in zip(serial, parallel, strict=True):
        for expected, actual in zip(expected_chunk, actual_chunk, strict=True):
            for left, right in zip(expected[:1], actual[:1], strict=True):
                if left is None:
                    assert right is None
                else:
                    assert_frame_equal(left, right, check_exact=True)
            for index in (1, 3, 4):
                assert len(expected[index]) == len(actual[index])
                for left, right in zip(expected[index], actual[index], strict=True):
                    assert_frame_equal(left, right, check_exact=True)
            assert expected[2] == actual[2]
            assert expected[5] == actual[5]
    assert parallel[1][0][-1].startswith("score: AttributeError:")
    if chemical_core:
        assert serial[0][0][-1] == ""
        summaries = serial[0][0][2]
        assert len(summaries) == 2
        assert summaries[0]["chemical_core_increment"] == summaries[1]["chemical_core_increment"]
        assert all(row["chemical_core_pairs"] == observations.height for row in summaries)
        for table, summary in zip(serial[0][0][1], summaries, strict=True):
            assert {"atom_count_increment", "atom_count_contribution", "chemical_core_increment",
                    "chemical_core_contribution"} <= set(table.columns)
            assert summary["atom_count_pairs"] == observations.height
            assert summary["atom_count"] == summary["original"]+summary["atom_count_increment"]
            assert summary["chemical_core"] == summary["original"]+summary["chemical_core_increment"]


def test_parallel_export_keeps_artifacts_and_propagates_pool_errors(tmp_path, monkeypatch):
    import gzip
    import json
    import tcren.annotation
    import tcren.annotation.batch
    import tcren.contacts.penalty_export as export
    import tcren.mhc
    import tcren.structure

    if available_cpus() < 2:
        pytest.skip("requires two allocated CPUs")
    structures = complex_fixture()
    observations = contact_observations(structures, interface="tcr_peptide")
    reference = pl.concat([observations, observations.with_columns((pl.col("dist")*.95).alias("dist"))])
    resource = tmp_path/"resource.json"
    resource.write_text(json.dumps(fit_resource(reference, observations)))
    sources = [tmp_path/f"synthetic-{i}.pdb" for i in range(4)]
    for path in sources:
        path.write_text("synthetic parser input")

    def parse(path, *, pdb_id, **kwargs):
        if pdb_id == "synthetic-1":
            raise ValueError("deliberate parser failure")
        return replace(structures, pdb_id=pdb_id)

    monkeypatch.setattr(tcren.structure, "parse_structure", parse)
    calls = []
    def annotate(structures, **kwargs):
        calls.append(len(structures))
        return [{} for _ in structures]
    monkeypatch.setattr(tcren.annotation.batch, "annotate_batch", annotate)
    monkeypatch.setattr(tcren.annotation, "classify_chains", lambda *a, **k: None)
    monkeypatch.setattr(tcren.mhc, "annotate_mhc_batch", lambda *a, **k: None)
    for threads in (1, 2):
        export.penalty_export([str(path) for path in sources], resource,
                              tmp_path/str(threads), threads=threads, chunk=3,
                              surface_context=True, surface_ownership="accessible_boundary")
    assert calls == [2, 1, 2, 1]
    for name in ("contacts.tsv.gz", "residues.tsv.gz"):
        with gzip.open(tmp_path/"1"/name, "rb") as left, gzip.open(tmp_path/"2"/name, "rb") as right:
            assert left.read() == right.read()
    for name in ("scores.tsv", "accounting.tsv", "annotation.jsonl", "tcren2.potential.tsv", "karnaukhov2022.potential.tsv"):
        assert (tmp_path/"1"/name).read_bytes() == (tmp_path/"2"/name).read_bytes()
    for path in (tmp_path/"1"/"surface").glob("*.parquet"):
        assert_frame_equal(pl.read_parquet(path), pl.read_parquet(tmp_path/"2"/"surface"/path.name), check_exact=True)
    def fail(*args):
        raise RuntimeError("deliberate pool failure")
    monkeypatch.setattr(export, "parallel_map", fail)
    with pytest.raises(RuntimeError, match="deliberate pool failure"):
        export.penalty_export([str(path) for path in sources], resource,
                              tmp_path/"failed", threads=2)
