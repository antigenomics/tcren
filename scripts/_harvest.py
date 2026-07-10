"""Shared junction harvest for the gap-placement scripts.

Pulls every ``C ... [FW]GXG`` omega-loop junction out of a structure set, labels each with its
chain type (TRA/TRB) and species by matching the recovered IMGT CDR3 against the curated
markup, and collapses crystal redundancy.

Redundancy is not bookkeeping. 374 structures yield ~377 junctions but only ~205 unique beta
sequences: the same clonotype is crystallised several times. Pairing without collapsing first
inflates every confidence interval, because the unit of independence is the junction, not the
pair.

2026-07-10
"""
from __future__ import annotations

import csv
import glob
import gzip
import os
import warnings
from dataclasses import dataclass

import numpy as np

warnings.filterwarnings("ignore")

AA = "ACDEFGHIKLMNPQRSTVWY"


@dataclass(frozen=True)
class Loop:
    """One junction, with the provenance needed to resample by junction rather than by pair."""
    pdb: str
    chain: str
    chain_type: str      # 'TRA' | 'TRB'
    species: str         # 'Human' | 'Mouse'
    seq: str             # AIRR junction, anchors included
    ca: np.ndarray

    @property
    def cdr3(self) -> str:
        return self.seq[1:-1]


def _markup(path: str) -> dict[str, dict[str, str]]:
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            out[row["pdb.id"].lower()] = row
    return out


def harvest(pattern: str, markup_csv: str, relax_length: bool = True) -> list[Loop]:
    """Every omega-loop junction in ``pattern``, typed and speciated from ``markup_csv``."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import index_to_one, three_to_index

    from tcren.loops import find_junctions, is_omega_loop

    mk = _markup(markup_csv)
    parser = PDBParser(QUIET=True)
    out: list[Loop] = []
    for path in sorted(glob.glob(pattern)):
        pdb = os.path.basename(path).split(".")[0].lower()
        row = mk.get(pdb)
        if row is None:
            continue
        try:
            with gzip.open(path, "rt") as fh:
                model = parser.get_structure(pdb, fh)[0]
        except Exception:
            continue
        for chain in model:
            res = [r for r in chain if r.id[0] == " " and "CA" in r]
            if not 90 <= len(res) <= 130:      # a TCR/Ig variable domain
                continue
            try:
                seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in res)
            except Exception:
                continue
            ca = np.array([r["CA"].get_coord() for r in res], dtype=float)
            for j in find_junctions(seq, ca):
                if not (set(j.seq) <= set(AA) and is_omega_loop(j.ca, relax_length=relax_length)):
                    continue
                # Type the chain by which curated CDR3 the recovered loop reproduces, rather
                # than by chain id -- chain letters are not consistent across depositions.
                if j.cdr3 == row.get("cdr3a"):
                    ct = "TRA"
                elif j.cdr3 == row.get("cdr3b"):
                    ct = "TRB"
                else:
                    continue
                out.append(Loop(pdb, chain.id, ct, row.get("species", "?"), j.seq, j.ca))
    return out


def collapse(loops: list[Loop]) -> tuple[list[Loop], dict[tuple[str, str], int]]:
    """One representative per (chain_type, sequence); also the crystal multiplicity of each.

    Returns ``(unique_loops, multiplicity)``. The multiplicity is worth reporting: if the same
    sequence's backbone varies as much across crystals as the effects being measured, the
    signal is inside crystallographic noise.
    """
    by_key: dict[tuple[str, str], list[Loop]] = {}
    for lp in loops:
        by_key.setdefault((lp.chain_type, lp.seq), []).append(lp)
    reps = [sorted(v, key=lambda x: (x.pdb, x.chain))[0] for _, v in sorted(by_key.items())]
    return reps, {k: len(v) for k, v in by_key.items()}


def crystal_noise_floor(loops: list[Loop]) -> list[float]:
    """CA-RMSD between crystals of the *same* junction sequence: the resolution limit."""
    from tcren.loops import kabsch_rmsd

    by_key: dict[tuple[str, str], list[Loop]] = {}
    for lp in loops:
        by_key.setdefault((lp.chain_type, lp.seq), []).append(lp)
    out = []
    for group in by_key.values():
        for a, b in ((group[i], group[j]) for i in range(len(group)) for j in range(i + 1, len(group))):
            if len(a.ca) == len(b.ca):
                out.append(kabsch_rmsd(a.ca, b.ca))
    return out
