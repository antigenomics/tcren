"""tcren — structure-based prediction of TCR recognition of epitopes.

A Python re-implementation of the TCRen method (Karnaukhov et al. 2024), extended to
the full TCR-pMHC picture: structure parsing, TCR/MHC annotation, canonical orientation,
contacts, and a configurable per-interface statistical potential, plus percentile rank,
fast ΔΔG, and a one-call oracle facade composing these for the paper notebooks.
"""

from . import potential
from .binder import is_real_interface
from .paper.helpers import annotate_batch
from .refine.anchors import native_peptide
from .cohort import (PHI_TERMS, Q_FEATURES_GEOM, coupling, phi_score, q_coupled, q_score,
                     strain_z, zscore)
from .clashes import ClashReport, has_clashes, interface_clashes
from .topology.pose import (
    POSE_FEATURES,
    POSE_FEATURES_CONTACT,
    POSE_FEATURES_DEGREE,
    POSE_FEATURES_SHELL,
    pose_consistency,
)
from .mechanics.stability import StabilityReport, contact_stability
from . import geometry, torsions
from .contactmap import ContactMap, ModeCentroid, binding_mode, registered_map
from .geometry import LoopInternalCoords, cdr3_internal_coords
from .torsions import cdr3_torsions, chain_torsions, residue_torsions
from .contacts import all_atom_contacts, ca_distance_matrix, peptide_internal_contacts
from .mechanics.dynamics import Stability, peptide_stability, stability_table
from .energetics.rotamers import contact_probabilities, repack, soft_energy
from .stacking import RING_ATOMS, ring_stacking
from .topology.surface import (SurfaceMap, surface_distance, surface_map, surface_stats,
                      surface_table, surface_tree)
from .cpl import (ResponseMatrix, equimolar_effect, mutation_effect, position_scan,
                  response_matrix)
from .energetics.mutation import alanine_scan, ddg, neoantigen_ddg, reference_delta
from .mechanics.springs import coupling_residues, interface_springs, rupture, stiffness_tensor
from .oracle import summarize_structure
from .docking import substitute_tcr
from .pipeline import PipelineResult
from .pipeline import run as run_pipeline
from .potential import Potential, derive_tcren, derive_tcren_loo
from .refine import check_register, fix_register, refine_peptide, substitute_peptide
from .refine.interface import interface_energy
from .energetics.scoring import (RecognitionMatrix, intra_peptide_energy, recognition_matrix, score_peptides,
                      score_structures)
from .scoring_rank import background_peptides, percentile_rank
from .structure import Structure, import_structure, parse_structure
from .structure.io import mean_bfactor

# One source of truth: pyproject.toml. Hard-coding it here as well meant a release could ship with
# `tcren info` reporting the previous version, which is what happened at 2.3.2.
#
# An EDITABLE install adds a second way to be wrong: the dist-info is written once, at install time,
# so after a version bump `importlib.metadata` keeps reporting the version that was current when
# `pip install -e .` last ran. It read 2.13.0 against a source tree at 2.25.0 -- twelve releases
# stale, and every provenance stamp written in between recorded the wrong version. When the package
# is imported from a source checkout, pyproject.toml is therefore the authority.
def _resolve_version() -> str:
    from pathlib import Path as _P
    src = _P(__file__).resolve().parents[2] / "pyproject.toml"
    if src.exists():
        for line in src.read_text().splitlines():
            if line.startswith("version"):
                return line.split("=", 1)[1].strip().strip('"\'')
    try:
        from importlib.metadata import PackageNotFoundError, version as _pkg_version
        return _pkg_version("tcren")
    except PackageNotFoundError:  # pragma: no cover - only an uninstalled source tree
        return "0.0.0+unknown"


__version__ = _resolve_version()

__all__ = [
    "annotate_batch", "native_peptide",
    "q_score", "phi_score", "q_coupled", "coupling", "Q_FEATURES_GEOM", "PHI_TERMS",
    "strain_z", "zscore",
    "potential",
    "Potential",
    "derive_tcren",
    "derive_tcren_loo",
    "parse_structure",
    "import_structure",
    "mean_bfactor",
    "Structure",
    "all_atom_contacts",
    "peptide_internal_contacts",
    "ring_stacking", "contact_probabilities", "soft_energy", "repack",
    "Stability", "peptide_stability", "stability_table",
    "SurfaceMap", "surface_map", "surface_stats", "surface_table",
    "surface_distance", "surface_tree",
    "RING_ATOMS",
    "ca_distance_matrix",
    "ContactMap", "ModeCentroid", "binding_mode", "registered_map", "geometry", "torsions",
    "cdr3_internal_coords", "LoopInternalCoords",
    "score_peptides", "recognition_matrix", "RecognitionMatrix", "intra_peptide_energy",
    "score_structures",
    "percentile_rank",
    "background_peptides",
    "ddg",
    "alanine_scan",
    "neoantigen_ddg",
    "reference_delta",
    "response_matrix", "ResponseMatrix", "mutation_effect", "position_scan", "equimolar_effect",
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
    "pose_consistency",
    "POSE_FEATURES",
    "POSE_FEATURES_CONTACT",
    "POSE_FEATURES_SHELL",
    "POSE_FEATURES_DEGREE",
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
