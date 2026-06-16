"""MHC mapping and partitioning."""

from .linker import (
    MhcAlignmentCheck,
    check_against_mhc,
    detect_linked_peptide,
    split_linked_peptides,
)
from .mapper import MhcCall, apply_mhc_calls, map_mhc
from .reference import build, reference_fasta
from .regions import annotate_mhc, annotate_mhc_batch, partition_chain, partition_mhc

__all__ = [
    "MhcCall", "map_mhc", "apply_mhc_calls", "build", "reference_fasta",
    "annotate_mhc", "annotate_mhc_batch", "partition_chain", "partition_mhc",
    "check_against_mhc", "detect_linked_peptide", "split_linked_peptides",
    "MhcAlignmentCheck",
]
