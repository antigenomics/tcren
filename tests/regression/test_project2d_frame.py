"""Projection-frame tests (slow — needs arda, MHC reference, Native2026 references)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("arda")

pytestmark = pytest.mark.slow  # invokes arda / mmseqs per structure

from scipy.spatial.distance import pdist

from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.paths import native_dir
from tcren.project2d import project_structure
from tcren.structure import parse_structure

REPO = Path(__file__).resolve().parents[2]
_HAVE_REF = (REPO / "database" / "mhc" / "alleles.aa.fasta").exists()
_HAVE_1AO7 = native_dir().is_dir() and any(native_dir().glob("1ao7.*"))
needs_data = pytest.mark.skipif(
    not (_HAVE_1AO7 and _HAVE_REF), reason="Native2026 references / MHC reference not built"
)


def _annotated_1ao7():
    path = next(native_dir().glob("1ao7.*"))
    s = parse_structure(path, pdb_id="1ao7")
    classify_chains(s, organism="human")
    annotate_mhc(s)
    return s


@needs_data
def test_native_and_pca_frames_agree_up_to_rigid():
    s = _annotated_1ao7()
    native = project_structure(s)
    pca = project_structure(s, force_pca=True)
    assert native.frame == "native" and pca.frame == "pca"
    assert native.keys == pca.keys and len(native.keys) > 50
    # Both are rigid transforms of the same Cα set → identical pairwise distances.
    d_native = pdist(native.coords3d)
    d_pca = pdist(pca.coords3d)
    assert np.allclose(d_native, d_pca, atol=1e-6)


@needs_data
def test_peptide_between_helices_and_below_cdrs():
    s = _annotated_1ao7()
    proj = project_structure(s)
    height = dict(zip(proj.keys, proj.height))
    pep = [height[(c.chain_id, r.seq_index)] for c in s.chains if c.chain_type == "PEPTIDE"
           for r in c.residues if (c.chain_id, r.seq_index) in height]
    cdr = []
    for c in s.chains:
        if c.chain_type in ("TRA", "TRB"):
            for reg in c.regions:
                if reg.region_type.startswith("CDR"):
                    cdr += [height[(c.chain_id, r.seq_index)] for r in reg.residues
                            if (c.chain_id, r.seq_index) in height]
    # In the canonical frame the TCR sits above the peptide (TCR on top, pMHC below).
    assert np.mean(cdr) > np.mean(pep)
