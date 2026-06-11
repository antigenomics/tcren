"""tcren — structure-based prediction of TCR recognition of epitopes.

A Python re-implementation of the TCRen method (Karnaukhov et al. 2024). Phase A
exposes the statistical-potential core; structure parsing, annotation, contacts and
scoring are added in later milestones.
"""

from . import potential
from .contactmap import ContactMap
from .contacts import all_atom_contacts, ca_distance_matrix
from .potential import Potential, derive_tcren, derive_tcren_loo
from .scoring import score_peptides, score_structures
from .structure import Structure, import_structure, parse_structure

__version__ = "0.1.0"

__all__ = [
    "potential",
    "Potential",
    "derive_tcren",
    "derive_tcren_loo",
    "parse_structure",
    "import_structure",
    "Structure",
    "all_atom_contacts",
    "ca_distance_matrix",
    "ContactMap",
    "score_peptides",
    "score_structures",
    "__version__",
]
