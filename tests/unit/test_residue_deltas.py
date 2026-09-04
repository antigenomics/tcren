"""Leave-one-residue-out attribution: the masking contract and the zero-contact invariant.

These check the CONSTRUCTION, never binder discrimination. The one that catches a real break is
`test_noncontacting_residue_moves_nothing`: a residue that touches no interface must leave every
descriptor where it was, so its delta is zero. If masking ever leaks -- a region markup that still
lists the dropped residue, a cached contact map, a mutated original -- that invariant is the first
thing to fail.
"""
from pathlib import Path

import numpy as np
import pytest

from tcren.annotation import classify_chains
from tcren.score.explain import _drop_residue, _interface_residues, residue_deltas
from tcren.structure.io import import_structure

PDB = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


@pytest.fixture(scope="module")
def s1ao7():
    pytest.importorskip("arda")            # classify_chains needs the arda backend
    s = import_structure(str(PDB))
    classify_chains(s, organism="human")
    return s


def test_drop_removes_one_residue_and_its_markup(s1ao7):
    c = next(c for c in s1ao7.chains if c.chain_type == "PEPTIDE")
    victim = c.residues[2]
    out = _drop_residue(s1ao7, c.chain_id, victim.seq_index)
    kept = out.chain(c.chain_id)
    assert len(kept.residues) == len(c.residues) - 1
    assert all(r.seq_index != victim.seq_index for r in kept.residues)
    assert all(r.seq_index != victim.seq_index for g in kept.regions for r in g.residues)
    # surviving indices keep their original numbering, gaps and all
    assert [r.seq_index for r in kept.residues] == \
           [r.seq_index for r in c.residues if r.seq_index != victim.seq_index]


def test_drop_does_not_mutate_the_original(s1ao7):
    c = next(c for c in s1ao7.chains if c.chain_type == "PEPTIDE")
    before = len(c.residues)
    _drop_residue(s1ao7, c.chain_id, c.residues[0].seq_index)
    assert len(s1ao7.chain(c.chain_id).residues) == before
    assert len(s1ao7.chains) == len(_drop_residue(s1ao7, c.chain_id, 0).chains)


def test_other_chains_are_shared_not_copied(s1ao7):
    """Only the edited chain is rebuilt; the rest are the same objects, which is what keeps the
    per-residue pass cheap enough to run over a whole interface."""
    c = next(c for c in s1ao7.chains if c.chain_type == "PEPTIDE")
    out = _drop_residue(s1ao7, c.chain_id, c.residues[0].seq_index)
    for a, b in zip(s1ao7.chains, out.chains):
        assert (a is b) == (a.chain_id != c.chain_id)


def test_noncontacting_residue_moves_nothing(s1ao7):
    """A residue on no interface cannot change a descriptor, so its delta is exactly zero.

    Run on the peptide score, which is a sum over TCR:peptide contacts: a receptor residue that
    contacts nothing contributes no term to it.
    """
    from tcren.contactmap import ContactMap
    from tcren.score import peptide_score
    from tcren.descriptors.table import _featurise_families
    import polars as pl

    touching = {t for t in _interface_residues(ContactMap.from_structure(s1ao7), None,
                                               ("tcr_peptide", "tcr_mhc", "peptide_mhc"))}
    tra = next(c for c in s1ao7.chains if c.chain_type == "TRA")
    idle = next(r.seq_index for r in tra.residues if (tra.chain_id, r.seq_index) not in touching)

    fam = ["placement", "interface", "topology", "energetics", "potts", "kinetics"]
    rows = [_featurise_families("full", s1ao7, "human", fam, (7.0, 8.0)),
            _featurise_families("cut", _drop_residue(s1ao7, tra.chain_id, idle), "human", fam, (7.0, 8.0))]
    v = np.asarray(peptide_score(pl.DataFrame(rows, infer_schema_length=None)), float)
    assert v[0] == pytest.approx(v[1], abs=1e-9), \
        f"masking non-contacting {tra.chain_id}{idle} moved the peptide score by {v[0] - v[1]:.3e}"


@pytest.mark.slow
def test_deltas_table_shape_and_identity(s1ao7):
    """delta is exactly full - without, one row per targeted residue, baseline excluded."""
    d = residue_deltas(s1ao7, score="peptide", chain_types=("PEPTIDE",))
    pep = next(c for c in s1ao7.chains if c.chain_type == "PEPTIDE")
    assert d.height == len({(pep.chain_id, r.seq_index) for r in pep.residues
                            if (pep.chain_id, r.seq_index)} & set(
        _interface_residues(__import__("tcren.contactmap", fromlist=["ContactMap"])
                            .ContactMap.from_structure(s1ao7), ("PEPTIDE",),
                            ("tcr_peptide", "tcr_mhc", "peptide_mhc"))))
    assert d["score.full"].n_unique() == 1
    got = np.asarray(d["delta"], float)
    want = float(d["score.full"][0]) - np.asarray(d["score.without"], float)
    assert np.allclose(got, want, atol=0, rtol=0)
