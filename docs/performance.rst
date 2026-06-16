Performance
===========

Per-stage timings on a TCR-pMHC complex (1ao7), Apple M3, single thread. Reproduce with::

   $ RUN_BENCHMARK=1 pytest -k benchmark -s

========================================  ===========  ====================================
stage                                     time         notes
========================================  ===========  ====================================
parse a gzipped structure                 ~19 ms       ``.pdb.gz`` / ``.cif.gz``
contact map (5 Å, cKDTree)                ~9 ms        per structure
score 1000 candidate peptides             ~8 ms        ~8 µs/peptide (vectorised)
annotate (TCR + MHC), **batched**         ~213 ms/str  one mmseqs2 call for the whole set
peak RSS, single-structure pipeline       ~195 MB
========================================  ===========  ====================================

Threading model
---------------

Annotation (TCR chain typing + MHC mapping) is the only compute-heavy step. It is always run as a
**single batched mmseqs2 search** over every chain in the input set — mmseqs2 parallelises internally,
so it is never called per structure and never wrapped in Python threads (a fork-based pool would also
deadlock after mmseqs2/BLAS spawn their own threads). Batching amortises the fixed ~1.5 s mmseqs2
startup: ~213 ms/structure across a set, versus ~1.5 s/structure one at a time.

Threads (:func:`tcren.orient.run_folder`'s ``threads`` / ``tcren orient -t N``) are used **only** for
the embarrassingly-parallel, mmseqs-free stages — structure parsing, the Kabsch/SVD alignment, and
writing oriented files — and, by extension, any PyMOL/Rosetta/FlexPepDock rendering and relaxation.
