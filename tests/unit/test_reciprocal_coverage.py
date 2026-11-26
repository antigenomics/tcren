import numpy as np
import pytest

from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure


def test_reciprocal_coverage_counts_distinct_endpoints_and_preserves_empty_states():
    from tcren.contacts.reciprocal_coverage import reciprocal_coverage

    def residue(i, *atoms):
        return Residue(i, i, "", "A", "ALA", tuple(
            Atom(name, "C", np.array([x, 0., 0.])) for name, x in atoms
        ))

    tcr = residue(0, ("CA", 0), ("N", -10))
    peptide = [
        residue(0, ("CA", 4), ("CG", 5), ("C", 1)),
        residue(1, ("CA", 12)),
        residue(2, ("CA", 20), ("CG", 4.8)),
    ]
    mhc = residue(0, ("CA", 40))
    structure = Structure("pose", [
        Chain("Z", [tcr], "TRA", "TRAB",
              regions=[RegionMarkup("CDR3", 0, 0, "A", [tcr])]),
        Chain("B", [], "TRB", "TRAB"),
        Chain("A", peptide, "PEPTIDE", "PEPTIDE"),
        Chain("M", [mhc], "MHCa", "MHCI",
              regions=[RegionMarkup("HELIX_A1", 0, 0, "A", [mhc])]),
    ])
    counts = reciprocal_coverage(structure)["tcr_peptide"]
    assert (counts.from_available, counts.to_available) == (2, 4)
    assert (counts.from_engaged, counts.to_engaged) == (1, 2)
    assert (counts.contact_pairs, counts.in_shell_pairs, counts.clashing_pairs) == (4, 3, 1)
    assert counts.coverage == pytest.approx(0.5)

    assert counts.information_bits is None
    with_information = reciprocal_coverage(structure, partner_information=True)["tcr_peptide"]
    assert with_information.from_entropy_bits == with_information.to_entropy_bits == 0
    assert with_information.information_bits == pytest.approx(np.log2(4 / 3) / 2)

    structure.chain("A").residues = [residue(0, ("CA", 12))]
    empty = reciprocal_coverage(structure, partner_information=True)["tcr_peptide"]
    assert empty.information_bits is None
    assert reciprocal_coverage(structure)["tcr_peptide"].coverage == 0
    structure.chain("A").residues = [residue(0, ("CA", 13))]
    assert reciprocal_coverage(structure)["tcr_peptide"].coverage is None
