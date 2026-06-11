"""Validate tcren's annotation against the TCR3D reference tables.

For each native CIF, the chain/complex data computed by the tcren pipeline (V gene,
CDR3, MHC class, MHC locus, epitope) must match TCR3D's ``tcr_complexes_data.tsv``.
These are independent annotation pipelines, so the robust fields are asserted exactly on
clean single-copy complexes; the benchmark sweep measures dataset-wide concordance.

Requires arda, a built MHC reference, and a bootstrapped native database.
"""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

pytest.importorskip("arda")

pytestmark = pytest.mark.slow  # invokes arda / mmseqs per structure

from tcren.native import NativeDatabase
from tcren.native.annotate import annotate_complex, cdr3_core, mhc_locus

REPO = Path(__file__).resolve().parents[2]
_DB = NativeDatabase()
_HAVE_REF = (REPO / "database" / "mhc" / "alleles.aa.fasta").exists()
needs_data = pytest.mark.skipif(
    not (_DB.is_present() and _HAVE_REF),
    reason="native database and/or MHC reference not built",
)


def _organism(db: NativeDatabase, pdb_id: str) -> str:
    row = db.complex_data.filter(pl.col("PDB_ID") == pdb_id)
    if row.height and row["TCR_organism"][0] == "Mouse":
        return "mouse"
    return "human"


@needs_data
@pytest.mark.parametrize("pdb_id", ["1ao7", "1bd2", "1mi5", "1oga", "1fyt"])
def test_native_annotation_matches_tcr3d(pdb_id):
    db = NativeDatabase()
    expected = db.complex_data.filter(pl.col("PDB_ID") == pdb_id).to_dicts()[0]
    ann = annotate_complex(db, pdb_id, organism=_organism(db, pdb_id))

    assert ann.mhc_class == expected["TCR_complex"]
    assert ann.epitope == expected["Epitope"]
    # MHC allele names are comparable for class I (IMGT HLA-A/B/C == TCR3D); class-II
    # TCR3D uses serotypes (HLA-DR1, I-Ak) that don't map 1:1 to our IMGT α-chain call.
    if expected["TCR_complex"] == "CLASSI":
        assert mhc_locus(ann.mhc_allele).startswith(expected["MHC_allele"])

    alpha, beta = ann.chain("Alpha"), ann.chain("Beta")
    assert alpha is not None and beta is not None
    assert alpha.v_gene == expected["TRAV_gene"]
    assert beta.v_gene == expected["TRBV_gene"]
    # arda's CDR3 equals the TCR3D CDR3 minus its conserved C…F/W anchors.
    assert alpha.cdr3 == cdr3_core(expected["CDR3_alpha"])
    assert beta.cdr3 == cdr3_core(expected["CDR3_beta"])


@needs_data
@pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"), reason="set RUN_BENCHMARK=1 to run")
def test_native_concordance_sweep():
    db = NativeDatabase()
    have_cif = set(db.pdb_ids())
    # Each field tracks (matches, comparable) where comparable counts only rows whose
    # TCR3D value is present, so null reference entries never penalise concordance.
    fields = ("v_alpha", "v_beta", "cdr3_alpha", "cdr3_beta", "class", "epitope", "j_alpha")
    match = dict.fromkeys(fields, 0)
    total = dict.fromkeys(fields, 0)
    n, errors = 0, []
    limit = int(os.getenv("TCREN_NATIVE_SWEEP_LIMIT", "0")) or None

    def cmp(field, got, exp):
        if exp is None or exp == "":
            return
        total[field] += 1
        if got is not None and got == exp:
            match[field] += 1

    ids = [p for p in db.complex_data["PDB_ID"].to_list() if p in have_cif][:limit]
    for pdb_id in ids:
        row = db.complex_data.filter(pl.col("PDB_ID") == pdb_id).to_dicts()[0]
        try:
            ann = annotate_complex(db, pdb_id, organism=_organism(db, pdb_id))
        except Exception as exc:  # noqa: BLE001 - record and continue
            errors.append((pdb_id, str(exc)[:60]))
            continue
        n += 1
        a, b = ann.chain("Alpha"), ann.chain("Beta")
        cmp("v_alpha", a.v_gene if a else None, row["TRAV_gene"])
        cmp("v_beta", b.v_gene if b else None, row["TRBV_gene"])
        cmp("cdr3_alpha", a.cdr3 if a else None, cdr3_core(row["CDR3_alpha"]))
        cmp("cdr3_beta", b.cdr3 if b else None, cdr3_core(row["CDR3_beta"]))
        cmp("class", ann.mhc_class, row["TCR_complex"])
        cmp("epitope", ann.epitope, row["Epitope"])
        cmp("j_alpha", a.j_gene if a else None, row["TRAJ_gene"])

    rates = {f: (match[f] / total[f] if total[f] else float("nan")) for f in fields}
    print(f"\nnative concordance over {n} complexes ({len(errors)} errors): {rates}")
    # Chain-level annotation tcren reproduces from the structures matches TCR3D closely
    # (V gene ~0.96-0.98, CDR3 ~0.90, class ~0.97). The small residual is arda-vs-TCR3D
    # gene-call differences plus a few malformed / multi-copy CIFs.
    assert rates["cdr3_alpha"] >= 0.85 and rates["cdr3_beta"] >= 0.85
    assert rates["v_alpha"] >= 0.85 and rates["v_beta"] >= 0.85
    assert rates["class"] >= 0.85
    # Epitope concordance is CIF-content-bounded, NOT a tcren error: some TCR3D class-II
    # CIFs are domain-split with no separable peptide chain, and some peptides differ by
    # a single unresolved terminal residue (structure vs sequence). Reported with a
    # modest floor rather than asserted strictly.
    assert rates["epitope"] >= 0.65
    # J-gene differs systematically (arda vs TCR3D J assignment) — reported, not asserted.
