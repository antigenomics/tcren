"""The per-position recognition matrix — the CPL/motif-matrix generalisation of score_peptides.

The load-bearing correctness property: the matrix decomposes the full interface score exactly, so
summing each position's native-amino-acid energy must equal :func:`score_peptides` on the native
sequence. Everything else (shapes, sides) follows from that.
"""
from pathlib import Path

import numpy as np
import pytest

from tcren import recognition_matrix, score_peptides
from tcren.annotation import classify_chains
from tcren.contactmap import ContactMap
from tcren.potential import tcren as tcren_pot
from tcren.structure.io import import_structure

PDB = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def _cm():
    pytest.importorskip("arda")  # classify_chains needs the arda backend; CI installs tcren without it
    s = import_structure(str(PDB))
    classify_chains(s, organism="human")
    native = "".join(r.aa for c in s.chains if c.chain_type == "PEPTIDE" for r in c.residues)
    return ContactMap.from_structure(s), native


def test_peptide_side_matrix_decomposes_the_full_score():
    cm, native = _cm()
    pot = tcren_pot()
    rm = recognition_matrix(cm, pot, side="to")                       # peptide side = CPL-matrix analog
    aa = list(rm.aa)
    got = float(np.nansum([rm.energy[i, aa.index(k[3])]              # matrix marginals == full score
                           for i, k in enumerate(rm.positions) if k[3] in aa]))
    want = float(score_peptides(cm, [native], pot)["score"][0])
    assert got == pytest.approx(want, abs=1e-9)                       # exact decomposition (float order)


def test_default_side_is_tcr_for_tcr_peptide():
    cm, _ = _cm()
    rm = recognition_matrix(cm, tcren_pot())                          # default scans the TCR side
    assert rm.side == "from"
    assert rm.energy.shape[1] == 20 and len(rm.positions) == rm.energy.shape[0]
    # entries are NaN only for amino acids the potential leaves undefined (Cys pairs); all others finite
    cys = list(rm.aa).index("C")
    assert np.isfinite(np.delete(rm.energy, cys, axis=1)).all()


def test_cdr_filter_restricts_positions():
    cm, _ = _cm()
    full = recognition_matrix(cm, tcren_pot(), side="from", tcr_regions="all")
    cdr = recognition_matrix(cm, tcren_pot(), side="from", tcr_regions="cdr")
    assert len(cdr.positions) <= len(full.positions)
    assert all(k[1] in {"CDR1", "CDR2", "CDR3"} for k in cdr.positions)
