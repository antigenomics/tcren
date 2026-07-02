"""tcren — structure-based prediction of TCR recognition of epitopes.

A Python re-implementation of the TCRen method (Karnaukhov et al. 2024), extended to
the full TCR-pMHC picture: structure parsing, TCR/MHC annotation, canonical orientation,
contacts, and a configurable per-interface statistical potential, plus percentile rank,
fast ΔΔG, and a one-call oracle facade composing these for the paper notebooks.
"""

from . import potential
from .binder import BINDER_MODEL, binder_score
from .clashes import ClashReport, has_clashes, interface_clashes
from .contactmap import ContactMap
from .contacts import all_atom_contacts, ca_distance_matrix
from .ddg import alanine_scan, ddg, neoantigen_ddg
from .mechanics import coupling_residues, interface_springs, rupture, stiffness_tensor
from .oracle import summarize_structure
from .pipeline import PipelineResult
from .pipeline import run as run_pipeline
from .potential import Potential, derive_tcren, derive_tcren_loo
from .refine import check_register, fix_register, refine_peptide, substitute_peptide
from .scoring import score_peptides, score_structures
from .scoring_rank import background_peptides, percentile_rank
from .structure import Structure, import_structure, parse_structure

__version__ = "2.1.2"

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
    "percentile_rank",
    "background_peptides",
    "ddg",
    "alanine_scan",
    "neoantigen_ddg",
    "binder_score",
    "BINDER_MODEL",
    "interface_springs",
    "stiffness_tensor",
    "rupture",
    "coupling_residues",
    "interface_clashes",
    "has_clashes",
    "ClashReport",
    "summarize_structure",
    "run_pipeline",
    "PipelineResult",
    "substitute_peptide",
    "refine_peptide",
    "check_register",
    "fix_register",
    "__version__",
]
