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
from .cohort import (F_TERMS, Q_FEATURES_GEOM, coupling, f_invert_by_iptm, f_score, phi_bind, q_f,
                     q_f_iptm, q_coupled, q_iptm, q_score, strain_z, zscore)
from .clashes import ClashReport, has_clashes, interface_clashes
from .stability import StabilityReport, contact_stability
from . import geometry, torsions
from .contactmap import ContactMap, ModeCentroid, binding_mode, registered_map
from .geometry import LoopInternalCoords, cdr3_internal_coords
from .torsions import cdr3_torsions, chain_torsions, residue_torsions
from .contacts import all_atom_contacts, ca_distance_matrix, peptide_internal_contacts
from .cpl import (ResponseMatrix, equimolar_effect, mutation_effect, position_scan,
                  response_matrix)
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

# One source of truth: pyproject.toml. Hard-coding it here as well meant a release could ship with
# `tcren info` reporting the previous version, which is what happened at 2.3.2. The fallback covers
# a source tree that was never installed (running from a checkout without `pip install -e .`).
try:
    from importlib.metadata import PackageNotFoundError, version as _pkg_version
    __version__ = _pkg_version("tcren")
except PackageNotFoundError:  # pragma: no cover - only an uninstalled source tree
    __version__ = "0.0.0+unknown"

__all__ = [
    "annotate_batch", "native_peptide",
    "q_score", "q_iptm", "f_score", "q_f", "q_coupled", "coupling", "Q_FEATURES_GEOM", "F_TERMS",
    "q_f_iptm", "f_invert_by_iptm",
    "phi_bind", "strain_z", "zscore",
    "potential",
    "Potential",
    "derive_tcren",
    "derive_tcren_loo",
    "parse_structure",
    "import_structure",
    "Structure",
    "all_atom_contacts",
    "peptide_internal_contacts",
    "ca_distance_matrix",
    "ContactMap", "ModeCentroid", "binding_mode", "registered_map", "geometry", "torsions",
    "cdr3_internal_coords", "LoopInternalCoords",
    "score_peptides", "recognition_matrix", "RecognitionMatrix",
    "score_structures",
    "percentile_rank",
    "background_peptides",
    "ddg",
    "alanine_scan",
    "neoantigen_ddg",
    "reference_delta",
    "response_matrix", "ResponseMatrix", "mutation_effect", "position_scan", "equimolar_effect",
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
    "contact_stability",
    "cdr3_torsions", "chain_torsions", "residue_torsions",
    "StabilityReport",
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
