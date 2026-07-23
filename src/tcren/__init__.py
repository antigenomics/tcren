"""tcren — structure-based prediction of TCR recognition of epitopes.

A Python re-implementation of the TCRen method (Karnaukhov et al. 2024), extended to
the full TCR-pMHC picture: structure parsing, TCR/MHC annotation, canonical orientation,
contacts, and a configurable per-interface statistical potential, plus percentile rank,
fast ΔΔG, and a one-call oracle facade composing these for the paper notebooks.
"""

from . import potential
from .binder import BINDER_MODEL, binder_score, is_real_interface
from .paper.helpers import annotate_batch
from .refine.anchors import native_peptide
from .cohort import phi_bind, q_score, strain_z, zscore
from .clashes import ClashReport, has_clashes, interface_clashes
from . import geometry
from .contactmap import ContactMap, ModeCentroid, binding_mode, registered_map
from .geometry import LoopInternalCoords, cdr3_internal_coords
from .contacts import all_atom_contacts, ca_distance_matrix
from .ddg import alanine_scan, ddg, neoantigen_ddg, reference_delta
from .mechanics import coupling_residues, interface_springs, rupture, stiffness_tensor
from .oracle import summarize_structure
from .orient import substitute_tcr
from .pipeline import PipelineResult
from .pipeline import run as run_pipeline
from .potential import Potential, derive_tcren, derive_tcren_loo
from .refine import check_register, fix_register, refine_peptide, substitute_peptide
from .refine.interface import interface_energy
from .scoring import RecognitionMatrix, recognition_matrix, score_peptides, score_structures
from .scoring_rank import background_peptides, percentile_rank
from .structure import Structure, import_structure, parse_structure

__version__ = "2.2.3"

__all__ = [
    "annotate_batch", "native_peptide",
    "q_score", "phi_bind", "strain_z", "zscore",
    "potential",
    "Potential",
    "derive_tcren",
    "derive_tcren_loo",
    "parse_structure",
    "import_structure",
    "Structure",
    "all_atom_contacts",
    "ca_distance_matrix",
    "ContactMap", "ModeCentroid", "binding_mode", "registered_map", "geometry",
    "cdr3_internal_coords", "LoopInternalCoords",
    "score_peptides", "recognition_matrix", "RecognitionMatrix",
    "score_structures",
    "percentile_rank",
    "background_peptides",
    "ddg",
    "alanine_scan",
    "neoantigen_ddg",
    "reference_delta",
    "binder_score",
    "BINDER_MODEL",
    "is_real_interface",
    "interface_springs",
    "stiffness_tensor",
    "rupture",
    "coupling_residues",
    "interface_energy",
    "interface_clashes",
    "has_clashes",
    "ClashReport",
    "summarize_structure",
    "run_pipeline",
    "PipelineResult",
    "substitute_peptide",
    "substitute_tcr",
    "refine_peptide",
    "check_register",
    "fix_register",
    "__version__",
]
