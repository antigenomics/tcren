"""Shuffled decoys: wrong-TCR-on-real-pMHC negatives for TCR-recognition models.

Given a set of co-framed (oriented) TCR-pMHC complexes, keep each complex's **pMHC** (peptide + MHC) intact and
graft on a **different** complex's **TCR** — a within-MHC-class derangement, so no complex keeps its own TCR.
The result is a physically-implausible recognition complex: a real, correctly-presented pMHC with the wrong
TCR docked over it. Both the TCR:peptide and TCR:MHC interfaces are therefore mismatched, while the peptide:MHC
interface is untouched (a clean internal control: peptide:MHC energy is invariant under the graft).

These decoys are the true negatives that a one-class "plausible complex" density lacks. A classifier trained on
real (label 1) vs shuffled (label 0) learns TCR-recognition compatibility from structures alone, with **no
binding-assay labels** — useful as a general, label-free recognition prior and as a supplementary benchmark.

Inputs MUST be **canonically oriented** (all superposed into the common MHC frame): run
:func:`tcren.orient.run_folder` (``tcren orient``) or :func:`tcren.orient.superimpose` (``tcren superimpose``)
first. The graft is then a **direct chain replacement with no per-pair alignment** — deliberately unlike
:func:`tcren.orient.graft.substitute_tcr`, which superposes the donor MHC onto the *host* MHC pairwise. Because
every complex already sits in the one canonical frame, dropping in the donor TCR as-is lets it keep its own
native docking angle relative to the canonical MHC, so the decoy set spans the **real MHC–TCR docking-angle
variance** across the whole database rather than forcing every TCR onto one host's pose. Chains are typed by
:func:`tcren.annotation.classify_chains` + :func:`tcren.mhc.annotate_mhc`, so the graft is by chain *type*.

CLI: ``tcren shuffle -s oriented/ -o shuffled/ --n 10``.
"""
from __future__ import annotations

import random
from collections.abc import Iterable, Iterator
from dataclasses import replace
from pathlib import Path

from .structure.model import Structure

_TCR_TYPES = frozenset({"TRA", "TRB", "TRG", "TRD"})
_ID_POOL = "ABYZWVUTSRQP"


def _tcr_chains(s: Structure):
    return [c for c in s.chains if c.chain_type in _TCR_TYPES]


def _pmhc_chains(s: Structure):
    return [c for c in s.chains if c.chain_type not in _TCR_TYPES]


def mhc_class(s: Structure) -> str | None:
    """MHC class (``"MHCI"``/``"MHCII"``) from the annotated MHC-chain supertype, or ``None`` if unknown."""
    for c in s.chains:
        if c.chain_supertype in ("MHCI", "MHCII"):
            return c.chain_supertype
    return None


def graft_tcr(pmhc_source: Structure, tcr_source: Structure, pdb_id: str | None = None) -> Structure:
    """Build a decoy: the pMHC chains of ``pmhc_source`` + the TCR chains of ``tcr_source``.

    Both structures must already share a coordinate frame (be oriented into the canonical MHC frame). **No
    coordinate transform is applied** — chains are copied as-is, so the grafted TCR keeps its own native
    docking angle relative to the canonical MHC (see the module docstring for why this beats pairwise
    alignment). TCR chain ids that collide with a pMHC chain id are reassigned.

    Args:
        pmhc_source: the complex whose peptide + MHC (and their coordinates) are kept.
        tcr_source: the complex whose TCR is grafted on.
        pdb_id: id for the decoy (default ``"<pmhc>__tcr_<tcr>"``).

    Returns:
        The decoy :class:`~tcren.structure.model.Structure`.

    Raises:
        ValueError: if ``pmhc_source`` has no pMHC or ``tcr_source`` has no TCR.
    """
    pmhc = _pmhc_chains(pmhc_source)
    tcr = _tcr_chains(tcr_source)
    if not pmhc:
        raise ValueError(f"no pMHC chains in {pmhc_source.pdb_id!r}")
    if not tcr:
        raise ValueError(f"no TCR chains in {tcr_source.pdb_id!r}")
    used = {c.chain_id for c in pmhc}
    grafted = []
    for c in tcr:
        cid = c.chain_id
        if cid in used:
            cid = next((x for x in _ID_POOL if x not in used), cid)
            c = replace(c, chain_id=cid)
        used.add(cid)
        grafted.append(c)
    pid = pdb_id or f"{pmhc_source.pdb_id}__tcr_{tcr_source.pdb_id}"
    return Structure(pid, pmhc + grafted, pmhc_source.complex_species, tcr_source.cell_type)


def make_decoys(structures: Iterable[Structure], n_per: int = 10, within_class: bool = True,
                seed: int = 0) -> Iterator[Structure]:
    """Yield decoy structures: each input pMHC paired with ``n_per`` distinct *other* TCRs.

    Within each MHC class (if ``within_class``) the TCR sources are a random selection excluding the pMHC's own
    complex, so no decoy reproduces a real pairing. Reproducible for a given ``seed``.

    Args:
        structures: co-framed, chain-typed + MHC-annotated TCR-pMHC structures.
        n_per: decoys generated per input pMHC (capped at the pool size).
        within_class: only graft TCRs from complexes of the same MHC class.
        seed: RNG seed.

    Yields:
        Decoy :class:`~tcren.structure.model.Structure` objects.
    """
    structures = list(structures)
    rng = random.Random(seed)
    cls = [mhc_class(s) if within_class else "_all" for s in structures]
    by_cls: dict[str | None, list[int]] = {}
    for i, c in enumerate(cls):
        by_cls.setdefault(c, []).append(i)
    for i, s in enumerate(structures):
        if not _pmhc_chains(s):
            continue
        pool = [j for j in by_cls.get(cls[i], []) if j != i and _tcr_chains(structures[j])]
        if not pool:
            continue
        for j in rng.sample(pool, min(n_per, len(pool))):
            yield graft_tcr(s, structures[j])


def _load_annotated(struct_dir: str | Path, organism: str = "human") -> list[Structure]:
    """Parse + chain-type + MHC-annotate every structure in a folder (one mmseqs pass per organism)."""
    from .annotation import classify_chains
    from .annotation.arda_adapter import _import_arda
    from .mhc import annotate_mhc_batch
    from .paper.helpers import _batch_annotate
    from .structure import parse_structure, structure_id_from_path, structure_paths

    structs: list[Structure] = []
    for path in structure_paths(Path(struct_dir)):
        try:
            structs.append(parse_structure(path, pdb_id=structure_id_from_path(path)))
        except Exception:
            pass
    if not structs:
        return []
    records = _batch_annotate(structs, _import_arda())
    for i, s in enumerate(structs):
        classify_chains(s, organism=organism, autodetect_species=True, precomputed_records=records[i])
    annotate_mhc_batch(structs)
    return structs


def run_shuffle(struct_dir: str | Path, out: str | Path, n: int = 10, seed: int = 0,
                within_class: bool = True, organism: str = "human", compress: bool = False) -> int:
    """Load a folder of oriented complexes, generate ``n`` decoys per pMHC, write them, return the count."""
    from .structure.io import write_structure

    structs = _load_annotated(struct_dir, organism=organism)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    written = 0
    for decoy in make_decoys(structs, n_per=n, within_class=within_class, seed=seed):
        write_structure(decoy, out / f"{decoy.pdb_id}.pdb{'.gz' if compress else ''}")
        written += 1
    return written
