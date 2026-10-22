"""Bounded driver accounting, fresh raw hashing and spawn/serial consistency."""
import importlib.util
import json
from hashlib import sha256

import polars as pl
import pytest
from polars.testing import assert_frame_equal


def test_audit_module_exists():
    assert importlib.util.find_spec("tcren.contacts.class_audit") is not None


def setup_pipeline(monkeypatch, tmp_path, fail=False):
    from test_contact_class_features import fixture
    from tcren import annotation, mhc, provenance, structure
    from tcren.annotation import batch
    calls = []
    def parse(path, *, pdb_id, keep_modified):
        assert keep_modified is True
        if fail and pdb_id == "m0":
            raise ValueError("synthetic parse failure")
        return fixture(pdb_id)
    def annotate(structures, *, threads):
        calls.append(len(structures))
        return [[] for _ in structures]
    monkeypatch.setattr(structure, "parse_structure", parse)
    monkeypatch.setattr(batch, "annotate_batch", annotate)
    monkeypatch.setattr(annotation, "classify_chains", lambda *args, **kwargs: None)
    monkeypatch.setattr(mhc, "annotate_mhc_batch", lambda *args, **kwargs: None)
    monkeypatch.setattr(provenance, "annotation_indexes", lambda structures, previous=None: {"synthetic": "fingerprint"})
    rows = []
    for i in range(2):
        raw = tmp_path / f"raw{i}.pdb"
        raw.write_text("synthetic raw input\n")
        rows.append(dict(id=f"m{i}", label=i, receptor_cluster=f"r{i}", dataset="synthetic",
                         cohort="all", epitope="AAA", mhc="MHC", mhc_class="I",
                         template_coverage="covered", source=str(raw), sha256=sha256(raw.read_bytes()).hexdigest()))
    manifest = tmp_path / "manifest.tsv"
    pl.DataFrame(rows).write_csv(manifest, separator="\t")
    return manifest, calls


def test_fresh_bounded_driver_and_complete_evidence(monkeypatch, tmp_path):
    from tcren.contacts.class_audit import class_audit
    manifest, calls = setup_pipeline(monkeypatch, tmp_path)
    output = tmp_path / "audit"
    class_audit(manifest, output, chunk=1, threads=1)
    assert calls == [1, 1]
    assert set(pl.read_csv(output / "accounting.tsv", separator="\t")["status"]) == {"ok"}
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["complete"] and provenance["source_count"] == 2
    assert provenance["searched_indexes"] == {"synthetic": "fingerprint"}
    assert pl.read_parquet(output / "case001" / "model_features.parquet")["id"].n_unique() == 2
    assert (output / "case001" / "inventory" / "batch000.parquet").is_file()
    assert (output / "case001" / "pairs" / "batch001.parquet").is_file()
    assert set(pl.read_parquet(output / "overall.parquet")["view"]) == {"chemistry"}
    assert set(pl.read_parquet(output / "case001" / "statistics.parquet")["view"]) == {"site", "chemistry"}
    with pytest.raises(FileExistsError):
        class_audit(manifest, output, chunk=1, threads=1)


def test_parse_failure_retains_ledger_and_aborts_all_statistics(monkeypatch, tmp_path):
    from tcren.contacts.class_audit import class_audit
    manifest, calls = setup_pipeline(monkeypatch, tmp_path, fail=True)
    output = tmp_path / "failed"
    with pytest.raises(ValueError, match="failed.*accounting"):
        class_audit(manifest, output, chunk=2, threads=1)
    accounting = pl.read_csv(output / "accounting.tsv", separator="\t")
    assert accounting.filter(pl.col("id") == "m0")["error"][0].startswith("parse:")
    assert not json.loads((output / "provenance.json").read_text())["complete"]
    assert not list(output.glob("**/statistics.parquet"))
    assert not (output / "overall.parquet").exists()


def test_raw_hash_mismatch_is_not_parsed(monkeypatch, tmp_path):
    from tcren.contacts.class_audit import class_audit
    manifest, calls = setup_pipeline(monkeypatch, tmp_path)
    frame = pl.read_csv(manifest, separator="\t").with_columns(pl.lit("0" * 64).alias("sha256"))
    frame.write_csv(manifest, separator="\t")
    output = tmp_path / "failed"
    with pytest.raises(ValueError, match="failed.*accounting"):
        class_audit(manifest, output, threads=1)
    assert calls == []
    assert set(pl.read_csv(output / "accounting.tsv", separator="\t")["stage"]) == {"source"}


def test_worker_real_spawn_matches_serial():
    from test_contact_class_features import fixture
    from tcren.contacts.class_audit import _audit_chunk
    from tcren.runtime import available_cpus, parallel_map
    if available_cpus() < 2:
        pytest.skip("two allocated CPUs required")
    work = [[fixture("one")], [fixture("two")]]
    serial = parallel_map(_audit_chunk, work, 1)
    spawned = parallel_map(_audit_chunk, work, 2)
    for left, right in zip(serial, spawned, strict=True):
        for (lid, ltables, lerror), (rid, rtables, rerror) in zip(left, right, strict=True):
            assert (lid, lerror) == (rid, rerror)
            for a, b in zip(ltables, rtables, strict=True):
                assert_frame_equal(a, b)


def test_global_chemistry_union_retains_other_case_eligible_zeros(monkeypatch, tmp_path):
    from dataclasses import replace
    import numpy as np
    from test_contact_class_features import fixture
    from tcren import structure
    from tcren.contacts.class_audit import class_audit
    manifest, _ = setup_pipeline(monkeypatch, tmp_path)
    frame = pl.read_csv(manifest, separator="\t")
    extra = frame.with_columns((pl.col("id") + "extra").alias("id"), pl.lit("BBB").alias("epitope"))
    pl.concat([frame, extra]).write_csv(manifest, separator="\t")
    def parse(path, *, pdb_id, keep_modified):
        result = fixture(pdb_id)
        if pdb_id.endswith("extra"):
            chain = result.chain("T")
            r = chain.residues[1]
            chain.residues[1] = replace(r, atoms=(replace(r.atoms[0], coord=np.array([2., 0., 0.])),))
        return result
    monkeypatch.setattr(structure, "parse_structure", parse)
    output = tmp_path / "audit"
    class_audit(manifest, output, chunk=2, threads=1)
    first = pl.read_parquet(output / "case001" / "statistics.parquet")
    framework = first.filter((pl.col("view") == "chemistry") & (pl.col("binder.region.type") == "UNASSIGNED"))
    assert framework.height > 0
    assert framework["sum_C_binder"].sum() == framework["sum_C_nonbinder"].sum() == 0
    assert framework["n_eligible_binder"].min() == framework["n_eligible_nonbinder"].min() == 1
    assert set(framework["set_membership"]) == {"neither"}


def test_reduction_feature_blocks_bound_calls_and_preserve_exact_statistics(monkeypatch, tmp_path):
    from importlib import import_module
    from test_contact_class_features import fixture
    from tcren.contacts.class_features import features_for_structure, opportunity_rows
    from tcren.contacts.class_statistics import STRATUM, class_statistics
    audit = import_module("tcren.contacts.class_audit")
    manifest, _ = setup_pipeline(monkeypatch, tmp_path)
    monkeypatch.setattr(audit, "_FEATURE_BLOCK_SIZE", 2)
    seen = []
    def checked_statistics(rows, meta):
        seen.append(rows["feature"].n_unique())
        assert seen[-1] <= 2
        assert meta.height == 2
        return class_statistics(rows, meta)
    monkeypatch.setattr(audit, "class_statistics", checked_statistics)
    output = tmp_path / "blocked"
    audit.class_audit(manifest, output, threads=1)
    assert len(seen) > 2
    outputs = [features_for_structure(fixture(f"m{i}")) for i in range(2)]
    c, b, t, dictionary = [pl.concat([result[i] for result in outputs]).unique() for i in (2, 3, 4, 5)]
    chemistry = dictionary.filter(pl.col("view") == "chemistry")
    counts = c.join(chemistry.select("feature"), on="feature", how="inner")
    rows = opportunity_rows(counts, b, t, chemistry, ["m0", "m1"]).rename({"pdb.id": "id"})
    metadata = pl.read_csv(manifest, separator="\t")
    expected_per, expected_overall = class_statistics(rows, metadata)
    expected_overall = expected_overall.join(chemistry, on="feature", how="left")
    actual_overall = pl.read_parquet(output / "overall.parquet")
    assert_frame_equal(actual_overall.sort("feature"), expected_overall.sort("feature"), check_exact=False)
    actual_rows = pl.read_parquet(output / "chemistry_model_features.parquet")
    assert_frame_equal(actual_rows.sort("id", "feature"), rows.sort("id", "feature"))
    actual_per = pl.read_parquet(output / "case001" / "statistics.parquet").filter(pl.col("view") == "chemistry")
    expected_per = expected_per.join(chemistry, on="feature", how="left")
    assert_frame_equal(actual_per.sort(STRATUM + ["feature"]), expected_per.sort(STRATUM + ["feature"]), check_exact=False)


def test_feature_reduction_window_real_spawn_matches_serial(monkeypatch, tmp_path):
    from importlib import import_module
    from test_contact_class_features import fixture
    from tcren.contacts.class_features import features_for_structure
    from tcren.runtime import available_cpus
    if available_cpus() < 2:
        pytest.skip("two allocated CPUs required")
    audit = import_module("tcren.contacts.class_audit")
    manifest, _ = setup_pipeline(monkeypatch, tmp_path)
    meta = pl.read_csv(manifest, separator="\t")
    outputs = [features_for_structure(fixture(f"m{i}")) for i in range(2)]
    c, b, t, dictionary = [pl.concat([result[i] for result in outputs]).unique() for i in (2, 3, 4, 5)]
    dictionary = dictionary.filter(pl.col("view") == "chemistry").sort("feature")
    cases = [(tmp_path, meta, c, b, t, dictionary)]
    monkeypatch.setattr(audit, "_FEATURE_BLOCK_SIZE", 2)
    serial = list(audit._reduction_windows(cases, dictionary, meta, 1))
    spawned = list(audit._reduction_windows(cases, dictionary, meta, 2))
    assert len(serial) > 2
    for (i, d1, tables1), (j, d2, tables2) in zip(serial, spawned, strict=True):
        assert i == j
        assert_frame_equal(d1, d2)
        for a, b in zip(tables1, tables2, strict=True):
            assert_frame_equal(a, b, check_exact=False)
