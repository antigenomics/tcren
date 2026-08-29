Structure summary (oracle)
==========================

:func:`tcren.summarize_structure` is the one-call entry point the paper notebooks use:
it takes a single TCR–peptide–MHC structure and returns a bundle of ready-to-tabulate
frames by composing the pipeline (S1+S2), the percentile rank (S3) and the alanine scan
(S4). Nothing is re-derived — the facade only orchestrates those functions, so its
``scores`` frame is byte-identical to :func:`tcren.pipeline.run`'s scores.

API
---

The full signature and argument reference live with the module autodoc:
:func:`tcren.oracle.summarize_structure` (re-exported as ``tcren.summarize_structure``).

The five returned frames:

.. list-table::
   :header-rows: 1
   :widths: 14 16 70

   * - Key
     - Source
     - Contents
   * - ``scores``
     - S1+S2 (``run``)
     - One row of per-interface energies (``tcr_peptide``, ``tcr_mhc``, ``peptide_mhc``,
       ``total``) and ``rmsd`` when ``superimpose=True``.
   * - ``rank``
     - S3 (``percentile_rank``)
     - One row: the native peptide's energy and its ``rank_pct`` against a random
       pMHC background.
   * - ``ddg``
     - S4 (``alanine_scan``)
     - Per-position alanine scan (``pos``/``wt_aa``/``ddG``); empty unless ``alanine=True``.
   * - ``markup``
     - S1+S2 (``run``)
     - The per-residue region-markup table.
   * - ``contacts``
     - S1+S2 (``run``)
     - The annotated residue-contact table.

Example
-------

The script below turns a PDB into the five summary CSVs (default: the bundled ``1ao7``
fixture):

.. literalinclude:: ../scripts/summarize_structure_example.py
   :language: python
   :caption: scripts/summarize_structure_example.py

Run it with the activated ``tcren`` environment:

.. code-block:: console

   $ python scripts/summarize_structure_example.py complex.pdb summary/
   scores       1 x 6  -> summary/scores.csv
   rank         1 x 5  -> summary/rank.csv
   ddg          9 x 3  -> summary/ddg.csv
   markup     605 x 7  -> summary/markup.csv
   contacts   512 x 18 -> summary/contacts.csv

Command line
------------

The ``rank`` and ``ddg`` frames are also available as standalone CLI subcommands.

``tcren rank`` — percentile-rank a peptide's TCRen energy against a random pMHC
background. With no ``-c/--candidates`` it ranks each structure's own native peptide:

.. code-block:: console

   $ tcren rank -s complex.pdb -o rank.csv
   $ tcren rank -s complex.pdb -c candidates.txt --background 5000 --seed 1 -o rank.csv

The output carries ``complex.id``, ``peptide``, ``score`` (native energy), ``rank_pct``
(fraction of background scoring at least as well — lower energy is a better binder, so a
small ``rank_pct`` flags a strong binder) and ``n_background``. ``--background-source``
points at a FASTA/text file of epitopes to sample the background from instead of drawing
it uniformly at random.

``tcren ddg`` — fast ΔΔG of peptide mutations (virtual-matrix path; no atoms move).
``ddG = E(native) - E(mutant)``, and lower energy binds better, so a positive value is
stabilising -- the mutant improves on the native:

.. code-block:: console

   $ tcren ddg -s complex.pdb --native LLFGYPVYV --alanine-scan -o ddg.csv
   $ tcren ddg -s complex.pdb --native LLFGYPVYV --mutant LLFGYPVYA --mutant LLFAYPVYV -o ddg.csv

Pass exactly one of ``--alanine-scan`` (one row per position mutated to alanine) or one
or more ``--mutant`` (neoantigen mode). Both subcommands share the ``--interface``,
``--regions``, ``-p/--potential`` and ``--cutoff`` options with ``tcren score``.

Which interface, and why it matters for a library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--interface tcr_peptide`` (the default) answers a **recognition** question and is structurally
blind to presentation. On 1ao7 the C-terminal anchor substitution reads exactly zero, because the
TCR never touches that position:

.. code-block:: console

   $ tcren ddg -s 1ao7.pdb --native LLFGYPVYV --mutant LLFGYPVYA --interface tcr_peptide
   # ddG = 0.0000

``--interface complex`` sums both peptide-bearing interfaces — the TCR:peptide potential plus
``--mhc-potential`` (Miyazawa-Jernigan by default) over peptide:MHC — which is the convention
:func:`tcren.cpl.response_matrix` has always used for a response-matrix cell:

.. code-block:: console

   $ tcren ddg -s 1ao7.pdb --native LLFGYPVYV --mutant LLFGYPVYA --interface complex
   # ddG = -0.9740

Use ``complex`` to **rank** whole peptides from a combinatorial or activation library: the assay
fires only if the peptide is presented *and* the receptor engages, so a peptide whose anchors are
destroyed must score badly, and under ``tcr_peptide`` it scores like any other. The two effects are
**not separable** in a library that varies every position — report the per-interface terms beside
the complex rather than reading either alone as the mechanism.
