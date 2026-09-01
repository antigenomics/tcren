"""Batched featurisation: one structure set -> one row per structure.

The dispatch layer. It owns the two things a whole-set run needs and a single-structure call does
not: the **single** arda call per organism plus the single mmseqs MHC search that annotate the set,
and the process pool that featurises it. Which columns each family contributes is
:mod:`tcren.descriptors.catalogue`; how the interface block is computed is
:mod:`tcren.descriptors.compute`.
"""
from __future__ import annotations

import math
from collections.abc import Sequence

from .catalogue import DESCRIPTORS, FAMILIES, PEPTIDE_INTERNAL_FEATURES
from .compute import (
    _peptide_internal_columns,
    _placement_columns,
    _stability_clash_columns,
    _symmetry_columns,
    recognition_features,
)

def recognition_table(items, *, organism: str = "human", full: bool = False,
                      threads: int = 1, chunk: int = 64,
                      autodetect_species: bool = True, mechanics: bool = False,
                      include: Sequence[str] | None = None, radii: Sequence[float] = (7.0, 8.0),
                      _mmseqs_threads: int = 0) -> list[dict]:
    """Batched feature (+score) extraction for a whole set of TCR–pMHC structures.

    ``items`` is an iterable of ``(id, structure-or-path)``. The set is annotated with a **single**
    arda call per organism (:func:`tcren.paper.helpers._batch_annotate`) and a **single** mmseqs MHC
    search (:func:`tcren.mhc.annotate_mhc_batch`) — the dataset-scale path that avoids the per-structure
    annotation cost — then :func:`recognition_features` (``full=``) is extracted for each. This emits
    **descriptors only**: the fitted composites and cohort-relative scores that used to ride along
    here were removed in 2.26.0, and scoring is :func:`tcren.reliability.s_score` on the table.
    ``full`` also appends the
    intra-peptide columns :data:`PEPTIDE_INTERNAL_FEATURES` (``Phi_pep_int``, ``n_pep_int``) — the
    peptide's contact energy with itself, which the interface energies omit. Returns one row dict per
    structure (``complex.id`` + features); a structure that fails yields
    ``{"complex.id": id, "error": ...}`` so the batch stays resilient.

    The two stages run **in sequence** and never compete for the machine.

    *Search* is one arda call per organism plus one mmseqs MHC search, each given every core, over
    the whole set. *Featurisation* is where the time actually goes — a 100-pose probe spends 96 s
    there against 2.4 s of arda and 0.9 s of MHC search — and it is pure Python/numpy, so
    ``threads`` > 1 runs it in that many **worker processes**. The flag keeps its name for
    compatibility; it has always meant "how much of this machine may I use".

    It used to mean concurrent *threads* over ``chunk``-sized batches, which was the wrong shape
    twice over: the GIL serialised the 94 % of the work that dominates, and each batch spawned its
    own mmseqs, so N batches asked for N x cores. Sharding the same work across independent
    subprocesses was measured 8x faster, which is what this now does directly.

    ``chunk`` is retained for signature compatibility and is no longer used.

    ``autodetect_species`` searches ``organism`` **and** mouse so a mis-declared cohort is still
    typed correctly. That doubles the annotation cost, so pass ``False`` when the organism is known
    — it halves the mmseqs work and changes nothing else.

    ``mechanics`` appends the :mod:`tcren.mechanics` koff proxies (``n_spring``, ``S_tot``,
    ``K_tens``, ``K_shear``, ``aniso``, ``rupture_force``, ``rupture_work``, ``couple_*``) to the
    same rows. They need the same annotated structure the descriptors do, so computing them here
    costs only their own arithmetic — running ``tcren mechanics`` separately repeats the whole
    parse and both mmseqs searches, and returns a second table keyed differently.
    """
    import os as _os

    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..mhc import annotate_mhc_batch
    from ..annotation.batch import _batch_annotate
    from ..structure import Structure, import_structure

    items = list(items)
    ids, structs, rows = [], [], []
    for id_, src in items:
        try:
            structs.append(src if isinstance(src, Structure) else import_structure(src))
            ids.append(id_)
        except Exception as exc:  # noqa: BLE001
            rows.append({"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"})

    if structs:                       # stage 1: one arda call per organism + one MHC search, all cores
        cores = _mmseqs_threads or (_os.cpu_count() or 1)
        orgs = (organism, "mouse") if autodetect_species else (organism,)
        recs = _batch_annotate(structs, _import_arda(), organisms=orgs, threads=cores)
        for i, s in enumerate(structs):
            try:
                classify_chains(s, organism=organism, autodetect_species=autodetect_species,
                                precomputed_records=recs[i])
            except Exception:  # noqa: BLE001 - MHC-only / unannotatable chains stay unset
                pass
        annotate_mhc_batch(structs, threads=cores)

    # stage 2: featurisation, the part that actually costs (94 % of wall time on a 100-pose probe:
    # 96 s against 2.4 s of arda and 0.9 s of MHC search). It is pure Python/numpy, so processes.
    work = [(id_, s, organism, full, mechanics, include, tuple(radii))
            for id_, s in zip(ids, structs)]
    if threads > 1 and len(work) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=min(threads, len(work))) as ex:
            rows.extend(ex.map(_featurise_one, work, chunksize=max(1, len(work) // (threads * 4))))
    else:
        rows.extend(_featurise_one(w) for w in work)

    return rows


def _featurise_one(args) -> dict:
    """One structure -> one row. Module-level and self-contained so it pickles to a worker process.

    The structure arrives already annotated: chain typing and the MHC call are batch operations and
    belong to the single search in :func:`recognition_table`, not to a per-structure worker.
    """
    id_, s, organism, full, mechanics, include, radii = args
    if include is not None:
        return _featurise_families(id_, s, organism, include, radii)
    try:
        feats = recognition_features(s, organism=organism, full=full, annotate=False)
        row = {"complex.id": id_, **feats, **_stability_clash_columns(s), **_symmetry_columns(s)}
        if full:                              # the intra-peptide term costs a second contact map
            row.update(_peptide_internal_columns(s))
        if mechanics:
            from ..mechanics.springs import interface_mechanics
            row.update(interface_mechanics(s))
        return row
    except Exception as exc:  # noqa: BLE001
        return {"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"}


def _featurise_families(id_, s, organism: str, include, radii) -> dict:
    """One structure -> one row holding exactly the catalogued descriptors of the requested families.

    Only what is asked for is computed: ``tcren features -i topology`` never builds the energies, and
    ``-i placement`` never runs the spring network. The returned row is filtered against
    :data:`DESCRIPTORS`, so a column exists in the output if and only if the catalogue names it —
    which is what makes the families a partition of the feature table rather than a label on it.
    """
    want = set(include)
    unknown = want - set(FAMILIES)
    if unknown:
        raise ValueError(f"unknown feature families {sorted(unknown)}; expected {FAMILIES}")
    row: dict[str, float] = {}
    try:
        if want & {"placement", "interface", "energetics"}:
            row.update(recognition_features(s, organism=organism, full=True, annotate=False))
            row.update(_symmetry_columns(s), **_peptide_internal_columns(s))
        if want & {"interface", "kinetics"}:                 # clash + contact fragility share a pass
            row.update(_stability_clash_columns(s))
        if "placement" in want:
            row.update(_placement_columns(s))
        if "topology" in want:
            row.update(_footprint_columns(s, radii))
        if "potts" in want:
            from ..potts import score_structure
            row.update({k: v for k, v in score_structure(s).items() if k != "pdb.id"})
        if "kinetics" in want:
            from ..mechanics.springs import interface_mechanics
            row.update(interface_mechanics(s))
    except Exception as exc:  # noqa: BLE001 - keep the batch alive, one bad structure is one bad row
        return {"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"}
    keep = {n for n, (fam, _) in DESCRIPTORS.items() if fam in want}
    keep |= {f"fp_{k}_r{r:g}" for r in radii for k in ("b0", "b1", "chi", "b0_frac")} if "topology" in want else set()
    return {"complex.id": id_, **{k: v for k, v in row.items() if k in keep}}
