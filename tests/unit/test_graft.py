"""Unit tests for TCR grafting (:func:`tcren.docking.substitute_tcr`).

Fast synthetic tests build two toy complexes where the donor is a rigid-body rotation of the host,
so a correct MHC- or TCR-anchored superposition must map the donor's TCR back exactly onto the host
TCR frame (RMSD ≈ 0). This exercises the superposition, transform, chain assembly, and id-relabelling
without needing arda. A slow, arda-gated self-graft checks the real annotation path end to end.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.docking import substitute_tcr
from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure

# A fixed rotation (90° about z) + translation; the donor is the host under this rigid move.
_R0 = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
_T0 = np.array([10.0, 20.0, 30.0])


def _res(i: int, aa: str, xyz) -> Residue:
    return Residue(seq_index=i, pdb_index=i + 1, insertion_code="", aa=aa, resname="ALA",
                   atoms=(Atom("CA", "C", np.asarray(xyz, float)),))


def _chain(cid: str, ctype: str, seq: str, coords: np.ndarray, groove: bool = False) -> Chain:
    residues = [_res(i, seq[i], coords[i]) for i in range(len(seq))]
    ch = Chain(chain_id=cid, residues=residues, chain_type=ctype)
    if groove:  # MHC-mode superposition reads Cα from groove regions
        ch.regions = [RegionMarkup("GROOVE_FLOOR", 0, len(residues) - 1, seq, residues)]
    return ch


def _complex(pdb_id: str, transform=None) -> Structure:
    """A toy TCR:pMHC: peptide + MHCa(groove) + TRA + TRB, with well-spread (non-degenerate) Cα.

    ``transform`` (rot, tran) rigidly moves every atom — used to make the donor a moved host copy.
    """
    rng = np.random.default_rng(0)
    spec = [("C", "PEPTIDE", "SLYNTVATL", False), ("A", "MHCa", "GSHSMRYFYT", True),
            ("D", "TRA", "CAASFGDNSK", False), ("E", "TRB", "CASSPGQGAY", False)]
    chains = []
    for cid, ctype, seq, groove in spec:
        coords = rng.standard_normal((len(seq), 3)) * 6.0 + rng.standard_normal(3) * 20.0
        if transform is not None:
            rot, tran = transform
            coords = coords @ rot + tran
        chains.append(_chain(cid, ctype, seq, coords, groove=groove))
    return Structure(pdb_id=pdb_id, chains=chains, complex_species="Human")


def _ca(structure: Structure, ctype: str) -> np.ndarray:
    ch = next(c for c in structure.chains if c.chain_type == ctype)
    return np.array([r.ca for r in ch.residues])


@pytest.mark.parametrize("by", ["mhc", "tcr"])
def test_graft_recovers_host_frame_for_rigid_donor(by):
    # host and a rigidly-moved donor built from the SAME RNG stream -> identical local geometry.
    host = _complex("host")
    donor = _complex("donor", transform=(_R0, _T0))

    chimera = substitute_tcr(host, donor, by=by)

    # Result keeps host peptide + MHC and adds the donor TCR (2 chains) -> 4 chains, unique ids.
    types = sorted(c.chain_type for c in chimera.chains)
    assert types == ["MHCa", "PEPTIDE", "TRA", "TRB"]
    ids = [c.chain_id for c in chimera.chains]
    assert len(ids) == len(set(ids))

    # Both anchors recover the rigid transform, so the grafted TCR lands on the host TCR frame.
    for ctype in ("TRA", "TRB"):
        assert np.allclose(_ca(chimera, ctype), _ca(host, ctype), atol=1e-6), ctype
    # The pMHC side is the untouched host.
    for ctype in ("PEPTIDE", "MHCa"):
        assert np.allclose(_ca(chimera, ctype), _ca(host, ctype), atol=1e-9), ctype


def test_grafted_tcr_ids_are_relabelled_on_collision():
    host = _complex("host")
    donor = _complex("donor", transform=(_R0, _T0))
    # Force a collision: rename the donor TRA to the host peptide's id "C".
    next(c for c in donor.chains if c.chain_type == "TRA").chain_id = "C"

    chimera = substitute_tcr(host, donor, by="tcr")
    ids = [c.chain_id for c in chimera.chains]
    assert len(ids) == len(set(ids))                      # no duplicates
    assert {c.chain_id for c in chimera.chains if c.chain_type == "PEPTIDE"} == {"C"}  # host kept "C"


def test_invalid_by_raises():
    host = _complex("host")
    with pytest.raises(ValueError, match="by must be"):
        substitute_tcr(host, host, by="sideways")


def test_missing_chains_raise():
    host = _complex("host")
    no_pep = Structure("nopep", [c for c in host.chains if c.chain_type != "PEPTIDE"])
    with pytest.raises(ValueError, match="no peptide"):
        substitute_tcr(no_pep, host, by="tcr")
    no_tcr = Structure("notcr", [c for c in host.chains if c.chain_type not in ("TRA", "TRB")])
    with pytest.raises(ValueError, match="no TCR"):
        substitute_tcr(host, no_tcr, by="tcr")


# --- slow: real annotation path (self-graft on a crystal complex) -----------------------------------

_ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


@pytest.mark.slow
@pytest.mark.parametrize("by", ["mhc", "tcr"])
def test_self_graft_reproduces_crystal(by):
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.structure import parse_structure

    s = parse_structure(_ASSET, pdb_id="1ao7")
    classify_chains(s, organism="human")
    annotate_mhc(s)

    chimera = substitute_tcr(s, s, by=by)  # self-graft = identity superposition
    # 1ao7 = peptide + MHCa + B2M + TRA + TRB (5 chains); the TCR should land back on itself.
    assert len(chimera.chains) == len(s.chains)
    for ctype in ("TRA", "TRB"):
        host = next(c for c in s.chains if c.chain_type == ctype)
        graft = next(c for c in chimera.chains if c.chain_type == ctype)
        a = np.array([r.ca for r in host.residues if r.ca is not None])
        b = np.array([r.ca for r in graft.residues if r.ca is not None])
        assert np.sqrt(((a - b) ** 2).sum(1).mean()) < 0.5, ctype
