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
     - **TCRen2.** Redundancy-balanced derivation over the 374 ``Native2026`` crystals,
       weighting each structure on both the epitope and the receptor axis (``--balance
       both``). This is the matrix the TCRen2 manuscript reports.
   * - ``dfire2``
     - ``DFIRE2_potential.csv``
     - 400
     - Residue-level DFIRE2 over every inter-chain residue pair of the 374 ``Native2026``
       crystals. A *physics*-reference contact energy, against TCRen's selection-reference
       log-odds, and so an independent baseline rather than a variant.
   * - ``tcren2_dfire``
     - ``TCRen2_dfire_potential.csv``
     - 380
     - TCRen2 with the DFIRE distance and DFIRE2 rotation corrections added.
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

``classic`` and ``tcren2`` differ substantially — Pearson :math:`r` = 0.875 over the 380
shared cells, with a maximum absolute difference of 0.846 on a TCRen2 range of 2.88
(-1.14 to +1.74).
They are not interchangeable, and a score computed under one cannot be compared with a
score computed under the other.

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
       -o TCRen2_potential.csv

The second command runs the whole pipeline — parse, annotate, contacts, derive — over all
374 structures in roughly 20 seconds, so there is no reason to cache a derived matrix
rather than rebuild it.

DFIRE reference states
----------------------

A contact count discards two things the coordinates already hold: how far apart the two
residues are, and how they are turned relative to each other. DFIRE (Zhou & Zhou 2002)
supplies the radial reference for the first — in a *finite* globular system the number of
pairs at separation :math:`r` grows as :math:`r^{1.61}`, not as the ideal-gas :math:`r^2` —
and DFIRE2 (Yang & Zhou 2008) adds the orientation coordinate. ``tcren derive-dfire`` builds
both from the same structures and returns, per amino-acid pair,

``E0``
    the orientation-free DFIRE energy of a contact, the direct analogue of one TCRen cell;
``C_dist``
    the change in TCRen's own log-odds when each contact is weighted by the DFIRE volume
    element :math:`(r/r_c)^{-1.61}` rather than counted once;
``C_rot``
    :math:`-\mathrm{KL}(P(\cos\theta_a, \cos\theta_b \mid a, b) \,\|\, \mathrm{uniform})`,
    the orientational free energy a contact-only count cannot see. It is :math:`\le 0` by
    construction and its magnitude is the pair's orientational information in nats.

``DFIRE2 = E0 + C_rot``, and the corrected TCRen is ``TCRen2 + C_dist + C_rot``.

**Why the corrections transfer and the resolved potential does not.** A distance- and
orientation-resolved 20×20 potential needs an occupancy per pair *per distance bin per
orientation cell*. The 374 crystals hold about 8,000 TCR:peptide contacts, and on that
interface alone the median amino-acid pair has **11** orientable contacts — so with a count
floor set where the estimator is trustworthy, **not one of the 400 cells qualifies**. Pooled
over every inter-chain pair of every interface the same floor admits 202 of 400. The
corrections are one number per pair and are properties of packing geometry rather than of TCR
biology, so they are estimated on the wide sample (``--scope all``, the default) and added to
the sparse TCR:peptide derivation.

.. code-block:: console

   $ tcren derive-dfire --structure-dir "$TCREN_DATA_DIR/Native2026" \
       --scope all --emit dfire2 -o DFIRE2_potential.csv

   $ tcren derive-dfire --structure-dir "$TCREN_DATA_DIR/Native2026" \
       --scope all --emit corrected --correct tcren2 --terms dist,rot \
       -o TCRen2_dfire_potential.csv

``--emit corrections`` writes the three columns themselves, and ``--emit radial`` the
distance-resolved :math:`u(a, b, r)` behind them. Glycine has no Cα→Cβ direction, so its
cells carry no rotation term; a pair with fewer than
:data:`~tcren.potential.dfire.MIN_ORIENTED` orientable contacts is given zero rather than the
value its Miller–Madow-corrected divergence would suggest, because at those counts the
estimator's residual bias is the size of the effect.

Smoothing a sparse matrix
-------------------------

Split-half derivations of the 380-cell TCRen2 matrix -- two disjoint halves of the reference
crystals, correlated over the 231 cells with at least ten observations -- agree at Pearson
*r* = 0.45. The cells are **undersampled, not overfitted**, and the tail is worse than the median:
tryptophan, cysteine and methionine columns hold a handful of observations each.

``--smooth-beta`` applies substitution-matrix pseudocounts (:mod:`tcren.potential.smoothing`,
the Henikoff scheme used by PSI-BLAST): a cell's prior is the BLOSUM62-weighted average of the
cells that *were* observed, blended in with weight :math:`\beta / (n + \beta)`. If Ile:Leu
contacts are common, Val:Leu is not really unknown. The BLOSUM background is recovered from the
published scores by a linear solve rather than transcribed, and agrees with the usually quoted
values to within the rounding of those scores (Ala 0.082 against 0.074, Trp 0.012 against 0.013).

.. warning::

   **Do not tune** :math:`\beta` **on split-half reproducibility.** It rises monotonically with
   smoothing -- 0.45 at :math:`\beta = 0` to 0.56 at :math:`\beta = 100` -- for the trivial reason
   that both halves are being pulled toward the same prior, and it would reach 1 in the limit where
   the matrix contains no data at all. On the Yang/Garcia B*27:05 potency series the same smoothing
   takes Spearman(:math:`\Phi`, log EC50) from +0.60 to +0.11 at :math:`\beta = 10`. Choose
   :math:`\beta` on a held-out endpoint, and note that on this data set the answer was zero.

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
