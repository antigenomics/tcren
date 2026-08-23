Potentials
==========

A TCRen potential is a directed table of residue–residue contact energies: the value at
``(from, to)`` is the energy of a TCR-side residue ``from`` contacting a peptide-side
residue ``to``. Values are in dimensionless statistical-potential units, **not** in
:math:`k_\mathrm{B}T`, and the matrix is genuinely asymmetric — transposing it is not a
no-op.

Every matrix shipped in ``tcren`` is listed below, and every one derived by us records the
exact command that reproduces it. That contract is machine-readable in
``src/tcren/data/potentials.json`` and enforced by
``tests/regression/test_shipped_potentials.py``: a potential-shaped CSV in the data folder
with no manifest entry fails the suite, and so does a recipe that no longer reproduces its
file to within ``1e-9``.

Shipped matrices
----------------

.. list-table::
   :header-rows: 1
   :widths: 16 26 14 44

   * - Key
     - File
     - Cells
     - Provenance
   * - ``classic``
     - ``TCRen_potential.csv``
     - 380
     - The *Nat. Comput. Sci.* 2022 derivation, and the historical default returned by
       :func:`tcren.potential.tcren`. Derived from the committed oracle contact maps
       restricted to non-redundant structures.
   * - ``tcren2``
     - ``TCRen2_potential.csv``
     - 380
     - **TCRen2 — the default TCR:peptide potential since 2.11.0.** Redundancy-balanced
       derivation over the **362 fully annotated αβ** ``Native2026`` crystals (``--ab-only``),
       weighting each structure on both the epitope and the receptor axis (``--balance
       both``). This is the matrix the TCRen2 manuscript reports.
   * - ``mj1996``
     - ``MJ1996_contact_energies.csv``
     - 400
     - Miyazawa & Jernigan 1996 contact energies, published table. Used on the two
       presentation interfaces (TCR:MHC and peptide:MHC), where a selection-derived
       potential is inappropriate.
   * - ``mj_keskin``
     - ``MJ_Keskin_potentials.csv``
     - 400
     - Keskin *et al.* residue contact potentials, published table.

``classic`` and ``tcren2`` differ substantially — Pearson :math:`r` = 0.867 over the 380
shared cells, with a maximum absolute difference of 0.943 on a TCRen2 range of
2.95 (-1.25 to +1.70).
They are not interchangeable, and a score computed under one cannot be compared with a
score computed under the other. ``tcren2`` also changed at 2.11.0 (374 → 362 structures),
so scores are not comparable across that boundary either.

Cysteine is dropped from the ``from`` axis as a data convention, which is why the TCRen
matrices carry 380 cells (19 × 20) rather than 400.

Reproducing them
----------------

``classic``, from the committed oracle contact maps:

.. code-block:: console

   $ tcren derive-potential \
       -i tests/assets/oracle/data/contact_maps_PDB.csv \
       --summary tests/assets/oracle/data/summary_PDB_structures.csv \
       --nonred \
       -o TCRen_potential.csv

``tcren2``, from the reference crystals — fetch them once with ``tcren fetch-data``, then:

.. code-block:: console

   $ tcren derive-potential \
       --structure-dir "$TCREN_DATA_DIR/Native2026" \
       --balance both \
       --ab-only \
       -o TCRen2_potential.csv

``--ab-only`` keeps the 362 complexes that carry both CDR3s and a peptide. Without it the 12
single-chain, γδ and pMHC-only files are kept, and because ``--balance`` skips a structure with a
null on any axis, each of those would enter at weight 1.0 — the maximum.

The second command runs the whole pipeline — parse, annotate, contacts, derive — over all
374 structures in roughly 20 seconds, so there is no reason to cache a derived matrix
rather than rebuild it.

Weighting
---------

The PDB is redundant on more than one axis, and the two are comparable in size. Over the
374 ``Native2026`` crystals there are 230 distinct epitopes and 226 distinct receptors;
212 structures share their epitope with at least one other and 223 share their receptor,
with largest groups of 9 and 10 respectively. Correcting only for epitope leaves receptor
bias untouched, and the two corrections pull the matrix in measurably different directions
— an epitope-balanced and a receptor-balanced matrix agree at only Pearson *r* = 0.86,
which is *less* than either agrees with the unweighted matrix.

``--balance epitope|tcr|both``
    For structure :math:`i` with :math:`n_a(i)` structures sharing its value on axis
    :math:`a`, the weight is the mean of the per-axis inverse counts

    .. math::
        w_i = \\frac{1}{|A|} \\sum_{a \\in A} \\frac{1}{n_a(i)}

    The **mean**, not the product, is the point. A previously unseen receptor against a
    nine-times-crystallized epitope is a genuinely new recognition event, and scores
    :math:`(1/9 + 1)/2 = 0.556` rather than the :math:`1/9` a product rule would give it.
    A structure unique on every axis gets 1.0; a true re-solve, duplicated on all of them,
    gets :math:`1/n`. With one axis this is plain inverse frequency. Overall scale cancels
    in the log-odds, so no normalization is applied. The epitope axis keys on the peptide
    sequence and the receptor axis on the CDR3α/CDR3β pair jointly. Helpers:
    :func:`tcren.potential.balanced_weights` and its single-axis alias
    :func:`tcren.potential.epitope_weights`.

``--redundancy-t``
    A different operation: clusters αβ complexes by CDR3α/CDR3β/peptide *distance* and
    keeps one representative per cluster, so unlike exact-identity balancing it also
    catches near-duplicates such as point mutants — at the cost of a threshold, and of
    conflating the two axes into one. :func:`tcren.potential.cluster_weights` down-weights
    instead of excluding, keeping every structure's data.

Both feed ``derive_tcren``'s ``weights`` argument, which multiplies each structure's
contribution to the amino-acid pair counts.

Choosing one in code
--------------------

.. code-block:: python

   from tcren.potential import Potential, tcren

   classic = tcren()                                   # the 2022 default
   tcren2 = Potential.from_csv("TCRen2_potential.csv")  # the manuscript's matrix

Adding a matrix
---------------

Ship the CSV in ``src/tcren/data/`` and add an entry to
``src/tcren/data/potentials.json`` in the same commit, giving ``file``, ``description`` and
a ``source`` of ``contact-maps``, ``structure-dir``, ``published`` or ``unknown``. The first
two must carry the flags that reproduce the file; the regression test will run them. A new
entry with source ``unknown`` fails the suite unless it is also listed under
``known_unresolved``, which exists to record historical files whose provenance cannot be
recovered — not as a place to put new ones.
