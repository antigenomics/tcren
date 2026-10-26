Batching and performance
========================

Pass a whole structure set to one command. TCR annotation searches a batch per
requested organism; MHC annotation is a separate batch. Never run one search
per structure or wrap the internally threaded search engine in a Python pool.
``screen`` reuses its annotated structures for descriptor calculation.

.. code-block:: console

   $ tcren screen -s models.tar.gz --confidence confidence.tsv -t 0 -o scores.tsv

For ``screen``, ``-t 0`` is the default and selects the allocated CPU budget;
``-t 1`` is explicit serial execution and a positive count bounds the budget.
Other commands expose their own defaults in ``--help``. Use one parallel layer:
search and coarse descriptor workers run in separate phases, with native worker
thread limits. Do not launch several all-core commands concurrently.

What limits scaling
-------------------

Annotation startup, file parsing and serial reductions can dominate small inputs.
Recognition/topology and contact counting use bounded coarse work; the structural
energy/reference loop retains serial work. Thus additional CPUs do not imply
linear speedup, particularly for a single structure. Canonical orientation,
shape calculations and raw structural exports have different costs; compare
like-for-like commands and output quantities.

Memory depends on input atom counts, references and active workers. There is no
fixed workstation-size assumption or universal per-model memory claim. Select
inputs before loading, keep archives compressed, and avoid retaining duplicate
expanded coordinates or calculated feature caches. For large sets, use the
command's documented chunk/batch options and a CPU budget that fits measured
memory. Do not trade a failed annotation for a silently omitted input.

Measure the actual workload
---------------------------

Record the installed wheel/version, reference identities, exact command, input
count, CPU allocation, wall time and peak process-tree memory. Sum-of-process
RSS can double-count shared pages; label it accordingly. Measure live thread
counts as well as active compute threads to detect nested oversubscription.
Compare serial and bounded execution on the same raw input and check exact
membership plus numerical agreement before scaling to the full cohort.

Reference acquisition/index construction is setup work and should be reported
separately from steady reference-backed scoring. Fresh scientific reproduction
still rebuilds calculated descriptors from raw sources. A faster frozen-feature
replay is a different operation and cannot stand in for end-to-end timing.

Use :doc:`reproduction` for complete calculation/provenance boundaries and
:doc:`getting-started` for reference setup. Timing claims should name the actual
hardware and workload; measurements from an older release are not current
performance guarantees.
