Performance
===========

Per-stage wall time (best of *n*) on one TCR-pMHC complex (1ao7), Apple M-series, single thread.
Reproduce the core stages with::

   $ RUN_BENCHMARK=1 pytest -k benchmark -s

.. list-table::
   :header-rows: 1
   :widths: 46 14 46

   * - stage
     - time
     - notes
   * - parse a gzipped structure
     - ~17 ms
     - ``.pdb.gz`` / ``.cif.gz``
   * - contact map (5 Å, cKDTree)
     - ~9 ms
     - per structure
   * - score 1000 candidate peptides
     - ~11 ms
     - ~10 µs/peptide (vectorised)
   * - ΔΔG alanine scan (9-mer)
     - ~11 ms
     - virtual-matrix path; no atoms move
   * - binder P(bind) (features + model)
     - ~49 ms
     - native geometry, no external tool
   * - peptide refine (2000-step DOPE MC)
     - ~320 ms
     - knowledge-based rigid-body refinement
   * - annotate — MHC map (1 structure)
     - ~670 ms
     - one mmseqs2 search
   * - annotate — TCR (arda, 1 structure)
     - ~1.5 s
     - mmseqs2 startup-bound; ~0.2 s/str batched
   * - full pipeline (no superimpose)
     - ~2.2 s
     - parse → annotate → contacts → score
   * - superimpose onto the canonical DB (per query)
     - ~2.8 s
     - aligns to every same-class DB structure

Peak resident memory
--------------------

.. list-table::
   :header-rows: 1
   :widths: 42 12 46

   * - workload
     - peak RSS
     - notes
   * - single-structure pipeline (no orient)
     - ~200 MB
     - parse → annotate → contacts → score → refine
   * - + ``superimpose`` (loads canonical DB)
     - ~780 MB
     - holds the Canonical2026 reference set in RAM

Well under the 32 GB workstation ceiling either way. If you are only scoring or computing ΔΔG,
skip ``superimpose`` (``tcren pipeline --no-superimpose``) to stay at the ~200 MB working set.

Threading model
---------------

Annotation (TCR chain typing + MHC mapping) is the only compute-heavy step. It is always run as a
**single batched mmseqs2 search** over every chain in the input set — mmseqs2 parallelises internally,
so it is never called per structure and never wrapped in Python threads (a fork-based pool would also
deadlock after mmseqs2/BLAS spawn their own threads). Batching amortises the fixed ~1.5 s mmseqs2
startup: ~0.2 s/structure across a set, versus ~1.5 s/structure one at a time.

Threads (:func:`tcren.orient.run_folder`'s ``threads`` / ``tcren orient -t N``) are used **only** for
the embarrassingly-parallel, mmseqs-free stages — structure parsing, the Kabsch/SVD alignment, and
writing oriented files — and, by extension, any PyMOL/Rosetta/FlexPepDock rendering and relaxation.

Scaling a screen
----------------

Placing pMHCs against a TCR panel, refining, and scoring each complex is embarrassingly parallel:
references are annotated/oriented **once** (not per complex), so the hot loop is only
``refine + contacts + score`` (~0.3–0.9 s/complex here). Distribute complexes across cores/machines —
throughput scales linearly with no shared state.
