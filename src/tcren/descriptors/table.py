"""Batched featurisation: one structure set -> one row per structure.

The dispatch layer. It owns the two things a whole-set run needs and a single-structure call does
not: the **single** arda call per organism plus the single mmseqs MHC search that annotate the set,
and the process pool that featurises it. Which columns each family contributes is
:mod:`tcren.descriptors.catalogue`; how the interface block is computed is
:mod:`tcren.descriptors.compute`.
"""
from __future__ import annotations

from collections.abc import Sequence

from .catalogue import DESCRIPTORS, FAMILIES
from ..runtime import check_parallel_entrypoint, parallel_map, resolve_threads, single_threaded
from .compute import (
    _footprint_columns,
    _peptide_internal_columns,
    _placement_columns,
    _stability_clash_columns,
    _symmetry_columns,
    recognition_features,
)

@single_threaded
def recognition_table(items, *, organism: str = "human", full: bool = False,
                      threads: int = 1, chunk: int = 64,
                      autodetect_species: bool = True, mechanics: bool = False,
                      include: Sequence[str] | None = None, radii: Sequence[float] = (7.0, 8.0),
                      _mmseqs_threads: int = 0, annotate: bool = True) -> list[dict]:
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

    Search and featurisation run in sequence with one fixed ``threads`` CPU budget.
    Default 1 is serial; 0 resolves all allocated cores. Search receives that budget;
    descriptors use at most that many spawn workers over coarse contiguous chunks.
    Each worker uses one native compute thread. No numerical or annotation pool is
    nested inside the descriptor workers. Workers preserve parsed-input order;
    parsing failures remain explicit rows.

    ``annotate=False`` reuses chain and MHC annotations already attached to input structures.
    Parallel calls require an importable script with an ``if __name__ == "__main__"`` guard;
    use the CLI for batches or ``threads=1`` in notebooks and stdin sessions.

    ``chunk`` is retained for signature compatibility and is no longer used.

    ``autodetect_species`` searches human **and** mouse for either supported ``organism`` so a mis-declared cohort is still
    typed correctly. That doubles the annotation cost, so pass ``False`` when the organism is known
    — it halves the mmseqs work and changes nothing else.

    ``mechanics`` appends the :mod:`tcren.mechanics` koff proxies (``n_spring``, ``S_tot``,
    ``K_tens``, ``K_shear``, ``aniso``, ``rupture_force``, ``rupture_work``, ``couple_*``) to the
    same rows. They need the same annotated structure the descriptors do, so computing them here
    costs only their own arithmetic — running ``tcren mechanics`` separately repeats the whole
    parse and both mmseqs searches, and returns a second table keyed differently.
    """
    from ..annotation import classify_chains
    from ..annotation.arda_adapter import _import_arda
    from ..mhc import annotate_mhc_batch
    from ..annotation.batch import _batch_annotate
    from ..structure import Structure, import_structure

    threads = resolve_threads(threads)
    items = list(items)
    check_parallel_entrypoint(threads, len(items))
    cores = resolve_threads(_mmseqs_threads) if _mmseqs_threads else threads
    if cores > threads:
        raise ValueError("annotation threads cannot exceed the total CPU budget")
    ids, structs, rows = [], [], []
    for id_, src in items:
        try:
            structs.append(src if isinstance(src, Structure) else import_structure(src))
            ids.append(id_)
        except Exception as exc:  # noqa: BLE001
            rows.append({"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"})

    if structs and annotate:          # stage 1: one query batch per organism/MHC, within budget
        other = "mouse" if organism == "human" else "human"
        orgs = ((organism, other) if autodetect_species and organism in ("human", "mouse")
                else (organism,))
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
    rows.extend(parallel_map(_featurise_one, work, threads))

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
            # the published descriptors ride on the topology already built: the gap is the
            # difference of the two surface height fields, and the graph pair reuses the contact
            # map. Both fail soft to NaN, so an unmappable groove is a thin row, not a dead batch.
            from ..topology.literature import literature_features
            row.update(literature_features(s))
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
