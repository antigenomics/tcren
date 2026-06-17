"""Canonical TCR-pMHC orientation: MHC-frame superposition + chain renaming."""

from .chains import CHAIN_RENAME, rename_chains, select_primary_complex
from .docking import DockingAngles, crossing_incident_from_vector, docking_angles
from .exceptions import detect_reverse_dock
from .frame import CanonResult, build_canonical_frame, canonical_frame
from .pipeline import (
    align_to_canonical,
    canonicalize_structure,
    check_oriented_complex,
    run_folder,
    run_superimpose,
)
from .superimpose import superimpose

__all__ = [
    "CanonResult", "canonical_frame", "build_canonical_frame",
    "detect_reverse_dock", "CHAIN_RENAME", "select_primary_complex", "rename_chains",
    "canonicalize_structure", "align_to_canonical", "check_oriented_complex", "run_folder",
    "superimpose", "run_superimpose",
    "DockingAngles", "docking_angles", "crossing_incident_from_vector",
]
