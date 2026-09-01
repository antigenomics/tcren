"""Peptide coverage, normalized so class I and class II are on one scale.

A raw contacted-residue count and a fixed N/central/C band split are both length-confounded: three
contacted positions mean different things on a class I 8-mer and a class II 12-mer, and a third of
one is not a third of the other. Every column here divides by the peptide's own length, and the
anchors are found from the coordinates rather than assumed at fixed positions -- there is no
position index anywhere in the definition.

The check is that this recovers the canonical registers: **class I buries its two termini** and
**class II buries a gapped core**, both read off contacts alone.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.contactmap import ContactMap
from tcren.footprint import PEPTIDE_COVERAGE_FEATURES, footprint_features

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
#: class I at three lengths, and the class II gliadin complex.
CLASS_I = ("1fo0", "1ao7", "5jhd")
CLASS_II = ("4ozg",)


def _load(pdb: str):
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.structure import parse_structure

    s = parse_structure(PDB_DIR / f"{pdb}.pdb")
    classify_chains(s, organism="human", autodetect_species=True)
    annotate_mhc(s)
    s.pdb_id = pdb
    return s


def _accessibility(structure) -> np.ndarray:
    """Per position, the share of its contacts that face the receptor rather than the groove."""
    pep = next(c for c in structure.chains if c.chain_type == "PEPTIDE")
    length = len(pep.residues)
    cm = ContactMap.from_structure(structure, cutoff=5.0)

    def per_position(frame, column):
        out = np.zeros(length)
        if frame.height:
            tally = frame.drop_nulls(column).group_by(column).len()
            for pos, n in zip(tally[column].to_list(), tally["len"].to_list()):
                if pos is not None and 0 <= pos < length:
                    out[int(pos)] = float(n)
        return out

    n_tcr = per_position(cm.interface("tcr_peptide"), "pos.to")
    n_mhc = per_position(cm.interface("peptide_mhc"), "pos.from")
    both = n_tcr + n_mhc
    return np.divide(n_tcr, both, out=np.zeros(length), where=both > 0)


@pytest.fixture(scope="module")
def structures():
    return {p: _load(p) for p in CLASS_I + CLASS_II}


@pytest.mark.parametrize("pdb", CLASS_I + CLASS_II)
def test_every_coverage_column_is_finite_and_bounded(structures, pdb):
    """Lengths 8, 9, 10 and 12 all land on the same [0, 1] scale."""
    row = footprint_features(structures[pdb])
    for col in PEPTIDE_COVERAGE_FEATURES:
        assert np.isfinite(row[col]), f"{pdb}: {col} is not finite"
        assert 0.0 <= row[col] <= 1.0, f"{pdb}: {col} = {row[col]} is outside [0, 1]"


@pytest.mark.parametrize("pdb", CLASS_I)
def test_class_i_buries_both_termini(structures, pdb):
    """P1 and P-omega sit in the A and F pockets, so their contacts go to the groove, not the TCR.

    Found from the contact map: no position index enters the definition, and the canonical
    class I register comes out anyway.
    """
    a = _accessibility(structures[pdb])
    assert a[0] < a.mean(), f"{pdb}: first position is not groove-facing"
    assert a[-1] < a.mean(), f"{pdb}: last position is not groove-facing"


@pytest.mark.parametrize("pdb", CLASS_I)
def test_class_i_presents_its_middle_to_the_receptor(structures, pdb):
    """The solvent-exposed bulge between the anchors is what the TCR reads."""
    a = _accessibility(structures[pdb])
    peak = int(np.argmax(a))
    assert 0 < peak < len(a) - 1, f"{pdb}: accessibility peaks at a terminus ({peak})"


@pytest.mark.parametrize("pdb", CLASS_II)
def test_class_ii_buries_a_gapped_core_not_a_terminal_block(structures, pdb):
    """Class II holds its peptide by periodic pockets, so buried positions sit *between* free ones.

    The open groove lets the peptide run through and out at both ends, anchoring it at pockets
    spaced along the register rather than at the two termini. The signature is interleaving: at
    least one groove-held position with a receptor-facing position on each side.
    """
    a = _accessibility(structures[pdb])
    free = a > a.mean()
    interleaved = [
        i for i in range(1, len(a) - 1)
        if not free[i] and free[:i].any() and free[i + 1:].any()
    ]
    assert interleaved, f"{pdb}: no groove-held position lies between two receptor-facing ones"


@pytest.mark.parametrize("pdb", CLASS_II)
def test_class_ii_leaves_the_receptor_less_peptide_than_class_i(structures, pdb):
    """Its groove is open at both ends and grips more of the chain, so pep_free_frac runs lower."""
    class_ii = footprint_features(structures[pdb])["pep_free_frac"]
    class_i = [footprint_features(structures[p])["pep_free_frac"] for p in CLASS_I]
    assert class_ii < max(class_i)


def test_coverage_does_not_depend_on_a_fixed_band_split(structures):
    """Peptides of three different lengths give a centre near the middle, not a length artefact."""
    for pdb in CLASS_I:
        centre = footprint_features(structures[pdb])["pep_cov_centre"]
        assert 0.25 < centre < 0.85, f"{pdb}: centre {centre} is not on the peptide's middle"
