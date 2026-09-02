Beyond the contact sum
======================

Every score in :doc:`assess` and :doc:`potentials` is, at bottom, a sum of a residue-pair potential
over a contact list. Five things a TCR:pMHC interface does are invisible to such a sum, and ``tcren``
measures each of them with its own instrument. This page is what those instruments are for; the
signatures are in :doc:`tcren`.

What a contact potential can and cannot express
-----------------------------------------------

A contact energy is not purely an interaction: burying a residue against *any* partner costs
something that depends on that residue alone. :meth:`~tcren.potential.model.Potential.decompose`
separates the two exactly, and only the pair part ``J`` is beyond what a per-position model can
already write down.

.. code-block:: python

   from tcren.potential import mj, mj1996, mj_partition_energy

   d = mj1996().decompose()          # e(a,b) = mean + H(a) + H(b) + J(a,b), J double-centred
   d.h("F"), d.j("F", "W")           # one-body term; the genuinely pairwise remainder
   d.energy("F", "W")                # reassembles the original value

   f = mj1996().hydrophobicity_fit()  # C0 + C1(q_a + q_b) + C2 q_a q_b
   f.r2, f.eigenvalue_share           # 0.98 on MJ1996, 0.84 on the bundled mj
   mj_partition_energy()["F"]         # 4.37 -- MJ's own one-body scale (larger = more hydrophobic)

Where a potential has that shape the interaction term is only ``C2 q_a q_b``, so it **cannot prefer
one pair of side chains over another of equal hydrophobicity**. Both calls refuse a directed
potential: TCRen is TCR-to-peptide and must not be split this way.

Peptide conformational stability: what a contact model cannot see
------------------------------------------------------------------

A contact potential scores whichever conformation it is handed. It cannot tell a peptide whose own
side chains **hold** it in the TCR-facing conformation from one that merely happens to have been
modelled there -- both present the same contact list. :mod:`tcren.mechanics.dynamics` puts the
backbone in motion: it samples peptide phi/psi by Metropolis Monte Carlo against DOPE and reports how
far the peptide wanders, not a better pose.

.. code-block:: python

   from tcren import peptide_stability, stability_table

   peptide_stability(structure).rmsf        # ensemble spread, A -- larger = floppier
   stability_table([s1, s2])["delta_rmsf"]  # intra-peptide term ON vs OFF, paired

**The hypothesis it was built to test** (Sewell, 2026-08): intra-peptide interactions stabilise the
productive bulge a TCR reads, so a poor binder could make many contacts and still fail to stabilise
the productive peptide conformation -- which would explain why an additive contact model describes
some systems well and others badly.

Tested on the combinatorial-peptide-library set: about 160 best-binder and 160 worst-binder modelled
complexes for each of seven clones, 2,102 structures. AUC is best-against-worst discrimination.

.. list-table::
   :header-rows: 1
   :widths: 20 25 25

   * - clone
     - contact energy
     - stability
   * - ila1
     - 0.348
     - **0.862**
   * - 868
     - 0.537
     - **0.677**
   * - sb27
     - 0.570
     - **0.934**
   * - mel8
     - 0.690
     - **0.876**
   * - 4c6
     - **0.955**
     - 0.519
   * - 1e6
     - **0.973**
     - 0.707
   * - mel5
     - **0.974**
     - 0.859

**Stability beats the contact energy in 4 of the 4 clones where the contact model fails, and 0 of
the 3 where it works.** Mean AUC over the failing clones goes 0.536 to 0.837; over the working ones
the contact energy stays ahead, 0.967 against 0.695. Combining the two as a within-clone z-sum lifts
the mean AUC from 0.721 to 0.826, improved in 5 of 7 -- though with seven clones that paired test is
underpowered (Wilcoxon p = 0.22).

**The intra-peptide term is a switch, and flipping it does what the hypothesis says.** Removing the
peptide's contacts with itself lets the *best* binders' backbones wander further (delta rmsf
= +0.021 A, s.e. 0.005, 4.4 sigma) and leaves the *worst* binders unchanged (+0.002 A, s.e. 0.007);
best against worst p = 0.042. The same term sharpens the stability discrimination itself, by +0.024
AUC on average and in 5 of 7 clones (Wilcoxon p = 0.078).

So the **mechanism** is supported while the **system** originally guessed is not: 4c6 is one of the
clones the contact model handles well here (0.955), and the ones it fails on are ila1, 868, sb27 and
mel8. Caveats worth carrying: these are modelled structures, the Monte Carlo is knowledge-based
rather than molecular dynamics (no solvent, no force field, no time), delta rmsf is a mechanistic
signal and not a useful classifier on its own (AUC 0.526), and every clone-level test has n = 7.

Side-chain repack: what a local minimiser cannot do
-----------------------------------------------------

:func:`tcren.energetics.rotamers.repack` (native ``_relax.repack``) places every side chain in the
chi rotamer the DOPE potential prefers. The rigid-body refiner moves the peptide and leaves every chi
where it found it, so a full-atom model whose side chains a predictor placed keeps them -- which is
most of why a pairwise contact energy stops discriminating on generated poses.

.. code-block:: python

   from tcren import repack

   fixed, report = repack(structure)   # report: n_conformers, energy, p_best per residue

Like for like -- same wrong-rotamer input (chi1 rotated 120 degrees), same 33 to 42 side-chain atoms,
same crystal reference:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * -
     - peptide side-chain RMSD (A)
     - time
   * - input (wrong chi1)
     - 4.131
     - --
   * - ``repack``
     - **2.364**
     - **6 ms**
   * - OpenMM (anchor-restrained minimisation)
     - 4.133
     - 3,103 ms

OpenMM leaves them where they are, and that is not a defect in OpenMM: a local minimiser cannot cross
the torsional barrier between two rotamer basins, so relaxing clashes and re-sampling rotamers are
different operations and only a discrete packer does the second. Over eight structures the packer
recovers side-chain RMSD 3.93 to 1.66 A, 8 of 8 improved, median 6 ms.

It rotates the side chains a model **has**. It cannot rebuild ones
:func:`tcren.refine.substitute.substitute_peptide` stripped; that is side-chain *construction*, and
it is still open.

Footprint shape: what the contacts say before they are scored
---------------------------------------------------------------

Every other scorer here sums over contacts. :mod:`tcren.topology.footprint` reads the same contact
map as a **shape** -- which of the six CDR loops touched what, and whether the resulting footprint is
one connected patch. No potential, no reference structure, no fitted parameter, and **no canonical
orientation**: every descriptor is invariant under rigid motion, so unaligned inputs are fine. Only
chain typing and CDR markup are needed, which the command line does in one batched annotation pass.

Coverage is the composition over cells -- the 6 CDR loops times {peptide, MHC}, optionally splitting
the peptide into thirds -- summarised by the normalised Shannon entropy and by the Hill numbers
(`Hill 1973 <https://doi.org/10.2307/1934352>`_,
`Jost 2006 <https://doi.org/10.1111/j.2006.0030-1299.14714.x>`_), where ``D2`` is the *effective
number of engaged cells*. Topology joins the contacted pMHC residues at a C-alpha threshold and
builds the flag complex: ``fp_b0_*`` counts footprint patches and ``fp_b1_*`` its holes. Coverage and
topology are only weakly related, which is why they belong in one channel, and why that channel is
read as ``T``, a directional score against the native crystals rather than a hand-written z-sum.

.. code-block:: console

   $ tcren features -s structures/ -i topology -o shape.tsv

.. code-block:: python

   from tcren.topology.footprint import footprint_batch, footprint_features
   from tcren.reliability import t_score

   row = footprint_features(structure)   # one dict, the shape descriptors at the default two radii
   row["D2_pep24"], row["fp_b0_r7"], row["L_canon"]

   table = footprint_batch("structures/")  # polars frame, one row per structure
   T = t_score(table)                      # the shape score, fit-free, one row is enough

The cyclomatic number of the bipartite contact graph (``E - V + C``) is deliberately not offered:
with of order thirty contacts among of order thirty residues it is dominated by ``E`` and simply
tracks interface size. The patch count is scale-free instead.

Surface topology: what a TCR meets before it binds
----------------------------------------------------

A contact potential scores an interface that already exists. :mod:`tcren.topology.surface` describes
the pMHC *beforehand*: the peptide sits in a groove between two helices, and a TCR coming down meets
one surface, so the descriptor is a height field ``h(x, y)`` over that groove with hydropathy and
charge painted on. The method follows
`SURFMAP <https://doi.org/10.1021/acs.jcim.1c01269>`_ (surface shell, per-cell feature, 8-neighbour
smoothing, Manhattan map distance, hierarchical tree) and
`Protein Surface Topography <https://doi.org/10.1074/jbc.RA119.010494>`_ (centre the chart on the
functional site). A flat raster rather than SURFMAP's equal-area spherical chart, because the
TCR-facing surface is an open, near-planar patch that a plane does not distort.

.. code-block:: python

   from tcren import surface_map, surface_stats, surface_distance, surface_tree

   smap = surface_map(structure)            # channels: h, phobic, charge; source: peptide/helix/floor
   surface_stats(smap)["frac_above_ridge"]  # how much peptide surface clears the MHC helix crests
   ids, d = surface_distance([m1, m2, m3])  # pairwise map distance -> epitopes cluster

Two things worth knowing, because both were defects first.

**The frame is refit from every structure** -- z from the groove-floor plane normal, **y from the
peptide**, origin on the peptide centroid. The floor's own principal axis is *not* the groove axis
(its beta-strands run across the groove), which put the two helices diagonally across the map.
Because the frame is intrinsic, maps compare without prealigning the inputs, which is SURFMAP's
standing caveat.

**Heights come from ray casting in the groove frame**, not from Shrake-Rupley surface points. Sphere
sampling is fixed in global axes, so the same structure rotated gave a different map (median cell
moved 1.35 A, ``relief`` by 19 %). Ray casting is exactly equivariant and needs no probe test: the
highest surface in a column is by definition the one nothing is above.

**"Featureless" becomes a number.** Over the 374 Canonical2026 complexes (230 distinct epitopes), the
epitopes the literature *names* as featureless and as bulged separate completely.

.. list-table::
   :header-rows: 1
   :widths: 20 34 16 15 15

   * - epitope
     - source
     - rank by ``frac_above_ridge``
     - ``frac_above_ridge``
     - ``relief`` (A)
   * - LPEPLPQGQLTAY
     - EBV BZLF1 13-mer, HLA-B\*35 -- **bulged**
     - **2 / 230**
     - 0.749
     - 2.81
   * - HPVGEADYFEY
     - HCMV pp65 11-mer, HLA-B\*35:08 -- **bulged**
     - 5 / 230
     - 0.562
     - 3.59
   * - EPLPQGQLTAY
     - EBV BZLF1 11-mer, HLA-B\*35 -- **bulged**
     - 8 / 230
     - 0.416
     - 2.54
   * - LLFGYPVYV
     - HTLV-1 Tax, HLA-A\*02:01 -- prominent P5-Tyr
     - 46 / 230
     - 0.145
     - 2.15
   * - GILGFVFTL
     - influenza M1, HLA-A\*02:01 -- **featureless**
     - 139 / 230
     - **0.000**
     - 1.16
   * - TAFTIPSI
     - HIV RT 8-mer, HLA-B\*51:01 -- **featureless**
     - 205 / 230
     - **0.000**
     - 0.95

Five of the eight most-protruding epitopes are literature-named bulged HLA-B\*35 epitopes; both named
featureless ones have *no* peptide surface clearing the helix crest at all. Structure-level AUC is
1.000 on ``relief``, ``peak_to_valley`` and ``frac_above_ridge`` (p <= 0.001, 9 featureless against 5
bulged structures) -- though with two distinct epitopes per group that is a 2-against-2 comparison,
so the properly powered evidence is the trend over all 279 class-I structures: ``frac_above_ridge``
rises from 0.054 (8-mers) to 0.569 (13-mers), Spearman on ``relief`` +0.414, p = 5.5e-13.

``notebooks/surface_topology.py`` draws the elevation, charge and hydropathy maps and reproduces this
comparison.

Ring stacking: the geometry an identity cannot carry
------------------------------------------------------

A contact potential scores a pair by identity, so two rings face to face at 3.5 A score exactly like
the same two residues brushing past edge on. :func:`tcren.stacking.ring_stacking` measures the
difference and returns **no energy**:

.. code-block:: python

   from tcren import ring_stacking

   ring_stacking(structure, cutoff=7.5)  # centroid_distance, interplanar_angle, vertical, lateral

``interplanar_angle`` near 0 is face to face, near 90 edge to face; a parallel-displaced stack shows
a small ``vertical`` with a few angstrom of ``lateral``. Proline is included -- its pyrrolidine ring
packs face on against aromatics through CH-pi contacts.
