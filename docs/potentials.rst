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
       derivation over the **359 fully annotated αβ** ``Native2026`` crystals the RCSB still
       holds — the only
       kind ``derive-potential`` accepts — weighting each structure on both the epitope and
       the receptor axis (``--balance both``). This is the matrix the TCRen2 manuscript
       reports.
   * - ``tcren2_em``
     - ``TCRen2_sc_potential.csv``
     - 380
     - **TCRen2-EM — the self-consistent matrix, shipped beside TCRen2 and not in place of
       it.** The same estimator, reference state, pseudocount and crystals as ``tcren2``;
       what changes is what is counted. Every *available* residue pair contributes, weighted
       by the probability the coupled contact model gives it, and the matrix is re-derived and
       refitted until it stops moving
       (:func:`tcren.potts.fit_selfconsistent`, :func:`tcren.potential.tcren2_em`). It moves
       TCRen2 less than the 2022 → 2026 re-derivation did — Pearson *r* = 0.944 against
       TCRen2, against *r* = 0.867 between TCRen2 and ``classic`` — and the movement lands in
       the cells the crystal tally left thin.
   * - ``mj1996``
     - ``MJ1996_contact_energies.csv``
     - 400
     - Miyazawa & Jernigan 1996 contact energies, published table. Five of its 210 unique
       pairs differ from AAindex3 ``MIYS960101`` by 0.04–0.28 (correlating at 0.99978);
       left untouched and pinned by the tests, with ``aaindex("MIYS960101")`` as the
       curated alternative.
   * - ``mj_keskin``
     - ``MJ_Keskin_potentials.csv``
     - 800
     - Two matrices in one file, both **identified against AAindex3 on 2026-08-29 at 400
       of 400 cells exactly**: ``mj()`` is ``MIYS990106``, Miyazawa & Jernigan **1999**,
       and ``keskin()`` is ``KESO980101``, the solvent-mediated interfacial form. The "MJ"
       one had carried an unknown-provenance warning since 2026-08-11.
   * - ``bt1999``
     - ``BT1999_contact_energies.csv``
     - 400
     - Betancourt–Thirumalai, AAindex3 ``BETM990101``: Miyazawa–Jernigan re-referenced with
       **Thr as the reference solvent**, so every Thr entry is exactly ``0.00``. Parsed from
       the record rather than retyped.
   * - *(resource)*
     - ``aaindex3.txt``
     - 47 records
     - The whole of AAindex3, verbatim. See :ref:`aaindex-resource`.

``classic`` and ``tcren2`` differ substantially — Pearson :math:`r` = 0.867 over the 380
shared cells, with a maximum absolute difference of 0.943 on a TCRen2 range of
2.95 (-1.25 to +1.70).
They are not interchangeable, and a score computed under one cannot be compared with a
score computed under the other. ``tcren2`` also changed at 2.11.0 (374 → 362 structures) and
again at 3.3.0 (362 → 359, when three withdrawn PDB entries were dropped), so scores are not
strictly comparable across either boundary. The second change is small — the two matrices agree
at *r* = 0.9996 over all 380 cells with a largest shift of 0.104 — but it is still a different
matrix.

Cysteine is dropped from the ``from`` axis as a data convention, which is why the TCRen
matrices carry 380 cells (19 × 20) rather than 400.

Reading one out
---------------

Any table the scoring commands resolve can be written to a CSV under the same name:

.. code-block:: console

   $ tcren potential-table -p tcren2_em -o TCRen2_sc_potential.csv
   TCRen2-EM: 380 pairs over 20 residues -> TCRen2_sc_potential.csv

``-p`` takes a bundled name or a path, exactly as it does on ``score``, ``ddg``, ``cpl`` and
``rank``, and nothing is computed — the output is the shipped matrix in long form.

Reproducing them
----------------

``classic``, from the committed oracle contact maps:

.. code-block:: console

   $ tcren derive-potential \
       -i tests/assets/oracle/data/contact_maps_PDB.csv \
       --summary tests/assets/oracle/data/summary_PDB_structures.csv \
       --nonred \
       -o TCRen_potential.csv

For an explicit historical recipe, a CSV contact table may carry a numeric
per-contact weight column. The CLI exposes the existing derivation API:

.. code-block:: console

   $ tcren derive-potential -i weighted_contacts.csv --weight-column vote --keep-cys -o historical_potential.csv

Weights must be finite and nonnegative, with at least one positive value.
They multiply ``--balance`` weights if both are supplied. ``--keep-cys`` retains
receptor cysteine entries; ``--drop-cys`` removes them. Omitting both retains
the variant default (classic drops receptor cysteine; am keeps it).
These explicit overrides are unavailable with ``--loo``.

``tcren2``, from the reference crystals — fetch them once with ``tcren fetch-data``, then:

.. code-block:: console

   $ tcren derive-potential \
       --structure-dir "$TCREN_DATA_DIR/Native2026" \
       --balance both \
       -o TCRen2_potential.csv

That is the whole command, and it is the one the released matrix was derived with. There is no
flag controlling the α:β restriction because it is unconditional: a structure missing either CDR3
or the peptide is out of scope for this potential, and the command says so as it drops them. On
``Native2026`` it keeps **362 of 374** — the 12 it drops are 3 pMHC-only, 8 single-chain and one
γδ. It then drops the PDB entries the RCSB has withdrawn, listed with their reasons in
``data/pdb_withdrawn.csv``, which removes three more and leaves **359**. Two of those three,
``1jtr`` and ``4zez``, are superseded by ``1mwa`` and ``5jzi``, both already in the set, so keeping
them counted one complex twice in a set built to be non-redundant, and the redundancy weighting
cannot see it because the two carry different ids. The third, ``4euq``, has no successor: the RCSB
records that its TCR α chain was modelled from the wrong gene at positions 93–98, which is the
CDR3α region this potential reads. ``--keep-withdrawn`` restores all three, for auditing the
difference and nothing else. Making it optional would be a trap rather than a convenience, because ``--balance`` skips a
structure carrying a null on any weighting axis and the derivation then defaults it to weight 1.0,
the maximum, so the least usable files would have counted for the most.

The second command runs the whole pipeline — parse, annotate, contacts, derive — over all
374 structures in roughly 20 seconds, so there is no reason to cache a derived matrix
rather than rebuild it.

``tcren2_em``, from the same crystals, through the self-consistent loop:

.. code-block:: console

   $ tcren potts selfcon \
       -s "$TCREN_DATA_DIR/Native2026" \
       --balance both --margin 0.5 --damping 0.5 \
       --steps 25 --tol 5e-3 --draws 400 \
       -o TCRen2_sc_potential.csv \
       --model-out potts_tcr_peptide_sc.json \
       --trace-out TCRen2_sc_trace.tsv

This one is **minutes rather than seconds**, and it is the only shipped matrix that does not
reproduce bit for bit. Its E-step is a block-Gibbs sample, so the derived cells carry
Monte-Carlo error that no seed removes across a re-run at a different draw count; the manifest
records a ``tolerance`` of 5e-3 for it and ``reproduce: benchmark``, so
:mod:`tests.regression.test_shipped_potentials` checks it only under ``RUN_BENCHMARK=1``. Every
other entry still reproduces at 1e-9 and must keep doing so — a loose tolerance is right for a
sampled estimator and wrong for a closed-form one, where it would hide exactly the drift that
test exists to catch. :doc:`potts` derives the loop and states the reduction gate:
at ``--margin 0`` the crystal pins every site and the probability-weighted tally reproduces the
plain one exactly.

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
    gets :math:`1/n`. With one axis this is plain inverse frequency. No normalization is
    applied: a fixed pseudocount makes the absolute weight scale part of the estimator. The epitope axis keys on the peptide
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

.. _aaindex-resource:

Every published contact matrix: the AAindex3 resource
-----------------------------------------------------

``src/tcren/data/aaindex3.txt`` is the upstream AAindex3 flat file, bundled verbatim: 47
records over the 20 amino acids, each transcribed by GenomeNet's curators from a published
table. Bundling the file rather than a converted subset means the provenance *is* the
record, and adding a matrix to a comparison costs a string rather than a transcription.

.. code-block:: python

   from tcren.potential import aaindex, catalogue, entry, identify

   catalogue()                       # all 47, with kind / symmetry / mean / citation fields
   aaindex("MOOG990101")             # one as a Potential, ready to score with
   entry("MIYS960103").description   # the record itself, including the non-energy tables

Three things the catalogue is for:

**Kind.** Not every entry is an energy. 42 are contact energies; two are contact *counts*
(``TANS760102``, ``MIYS960103``) and three are side-chain centre *distances*
(``BONM030104``–``BONM030106``). :func:`~tcren.potential.aaindex` refuses the last two
kinds, because scoring a contact map with a count table is a silent category error;
:func:`~tcren.potential.entry` still returns them when that is what you want.

**Reference state.** Read the ``mean`` column. A matrix with mean near zero is a
*pair-contact* form with the one-body transfer term removed (``MIYS990106`` −0.079,
``BETM990101`` −0.057); one with a large negative mean is a *raw contact energy* that still
carries it (``KESO980101`` −3.547, ``MIYS960101`` −3.166). Comparing across the two groups
compares reference states as well as derivations.

**Symmetry.** The six ``ZHAC*`` entries are environment-dependent — row secondary structure
against column secondary structure — and three of them are asymmetric by construction, so
they are directed potentials and :meth:`~tcren.potential.Potential.decompose` refuses them.

One caveat is carried in the code as well as here: **AAindex's PMID field sometimes cites the
paper that tabulated a matrix rather than the one that derived it** — ``MIYS850102`` carries
Bastolla 2001. Check the entry's own author, title and journal fields before citing.

Identifying an unlabelled matrix
--------------------------------

:func:`~tcren.potential.identify` compares a potential cell by cell against every AAindex3
entry and returns the accessions ordered by maximum absolute difference. It exists because
two of the matrices shipped here arrived with no recorded upstream table:

.. code-block:: python

   >>> from tcren.potential import identify, mj, keskin
   >>> identify(mj())[:2]
   [('MIYS990106', 0.0), ('BASU010101', 0.6525)]
   >>> identify(keskin())[:2]
   [('KESO980101', 0.0), ('LIWA970101', 2.77)]

An identification needs both halves: an exact match *and* a distant runner-up. Here the
runner-up is off by 0.65 and 2.77 respectively, so neither is a coincidence.

Splitting a potential into what it actually measures
-----------------------------------------------------

An interface score is a sum over contacts, so the exact split
:math:`e(a,b) = \mathrm{mean} + H(a) + H(b) + J(a,b)` of
:meth:`~tcren.potential.Potential.decompose` carries straight through to the score.
:meth:`~tcren.potential.Potential.components` returns the three parts as scorable
potentials:

.. list-table::
   :header-rows: 1
   :widths: 12 34 54

   * - Part
     - Matrix
     - What its interface sum equals
   * - ``size``
     - the grand mean everywhere
     - ``mean × (number of contacts)`` — an interface-area term
   * - ``comp``
     - :math:`H(a) + H(b)`
     - a degree-weighted composition term
   * - ``pair``
     - :math:`J(a, b)`
     - the interaction proper, one-body parts removed

.. code-block:: python

   from tcren.potential import mj1996

   parts = mj1996().components()          # {"size": ..., "comp": ..., "pair": ...}

Scoring a structure with each part in turn says which of three very different things a
potential is reading on an interface: how *big* it is, what it is *made of*, or which
residue *faces which*. That matters because a matrix with no positive entries has a large
negative mean, so its interface sum is dominated by the contact count — and a result
obtained with one can be an interface-area effect wearing a chemical name.

Adding a matrix
---------------

Ship the CSV in ``src/tcren/data/`` and add an entry to
``src/tcren/data/potentials.json`` in the same commit, giving ``file``, ``description`` and
a ``source`` of ``contact-maps``, ``structure-dir``, ``published`` or ``unknown``. The first
two must carry the flags that reproduce the file; the regression test will run them. A new
entry with source ``unknown`` fails the suite unless it is also listed under
``known_unresolved``, which exists to record historical files whose provenance cannot be
recovered — not as a place to put new ones.

Three optional keys exist for recipes that are not ``derive-potential`` over a fixed table, and
``tcren2_em`` is the first entry to use any of them:

``command``
    The CLI subcommand, default ``derive-potential``.
``structure_flag``
    The flag that subcommand takes its structure directory on, default ``--structure-dir``.
    ``potts selfcon`` takes ``-s``.
``tolerance`` and ``reproduce``
    How close the rebuild has to land, default ``1e-9``, and whether the check is deferred to
    ``RUN_BENCHMARK=1``. **Set a tolerance above 1e-9 only when the estimator itself is
    sampled.** A closed-form derivation over a fixed contact table is deterministic and must
    reproduce bit for bit; loosening it there hides drift instead of tolerating noise. State the
    measured floor in ``tolerance_reason`` so the number is auditable rather than chosen.
