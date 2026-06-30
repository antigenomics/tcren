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
``ddG = E(native) - E(mutant)``, so a positive value is destabilising:

.. code-block:: console

   $ tcren ddg -s complex.pdb --native LLFGYPVYV --alanine-scan -o ddg.csv
   $ tcren ddg -s complex.pdb --native LLFGYPVYV --mutant LLFGYPVYA --mutant LLFAYPVYV -o ddg.csv

Pass exactly one of ``--alanine-scan`` (one row per position mutated to alanine) or one
or more ``--mutant`` (neoantigen mode). Both subcommands share the ``--interface``,
``--regions``, ``-p/--potential`` and ``--cutoff`` options with ``tcren score``.
