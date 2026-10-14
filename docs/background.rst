The amino-acid background, in 3D
================================

:func:`tcren.ddg.reference_delta` references a **specific** substitution: replace the chain with
poly-alanine, recompute the contact map, subtract. :func:`tcren.energetics.mutation.smoothed_reference`
references the **distribution** over all twenty instead, but does it *virtually* -- it holds the
contact map fixed and re-indexes the potential, so a contact a long arginine reaches across is still
charged to the alanine that notionally replaces it.

:mod:`tcren.energetics.background` takes that background structurally. Every candidate residue is
built at every rotamer of a chi grid, weighted by DOPE against the environment, and each
``(position, candidate, partner)`` resolves to one of **three** outcomes rather than one number.

Why alanine was the right specific reference
--------------------------------------------

:func:`tcren.refine.substitute.substitute_residues` truncates every target at Cbeta. Alanine's heavy
atoms are exactly ``N``, ``CA``, ``C``, ``O`` and ``CB``, so for an alanine target the truncation
*is* the mutant -- no rotamer, no relaxation, no choice. Every other target is a Cbeta stub carrying
the right label and the wrong reach.

That is why the shipped scores reference alanine, and the size of what the reference removes is
measured, not assumed. Over 370 Native2026 crystals and 8,002 TCR:peptide contacts at 5 A, only
**37.77 %** survive truncation to backbone plus Cbeta and **45.03 %** are made by side-chain atoms
on both sides at once.

Three outcomes, not one
------------------------

no contact
   No surviving rotamer of that candidate comes within ``cutoff`` of the partner. ``p = 0``.

contact
   Some surviving rotamer does. ``p`` is the Boltzmann-weighted fraction that do.

collision
   **Every** rotamer of that candidate overlaps the environment beyond ``clash_tol``. The candidate
   is inadmissible and leaves the background sum entirely, its prior mass renormalised over the
   survivors.

The third is the one the virtual scan cannot express, and it is not cosmetic. A plain
``w / w.sum()`` over all rotamers reports an impossible substitution at its *least-bad* rotamer and
lets it set the reference free energy. Here the free energy is taken over the admissible set
:math:`A_i`,

.. math::

   \varphi_i(a) = \sum_j p_{ij}(a)\, e(a, y_j), \qquad
   F_i = -\beta^{-1}\ln \sum_{a \in A_i} \tilde p_a\, e^{-\beta \varphi_i(a)}

and :math:`|A_i|` is itself a read-out: a position where seventeen of twenty residues collide is a
different statement about the interface from one where all twenty fit.

.. code-block:: python

   from tcren.energetics import background as bg
   from tcren.potential import tcren2

   scan, positions, partners = bg.scan_structure(structure, "tcr_peptide", "to")
   out = bg.reference_free_energy(scan, positions, partners, tcren2(), side="to")

   out["dPhi"]          # per position: the native residue's energy, referenced
   out["varPhi"]        # its spread over the residues that could actually be there
   out["n_admissible"]  # |A_i|

What it looks like on a crystal
--------------------------------

Over 60 Native2026 crystals and 596 peptide positions -- 585,720 ``(pair, candidate)`` cells:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - outcome
     - cells
     - share
   * - collision
     - 84,826
     - 14.48 %
   * - no contact
     - 477,259
     - 81.48 %
   * - contact
     - 23,635
     - 4.04 %

:math:`|A_i|` reads mean **17.17** of 20 (p05 12, p95 20); 24.3 % of positions admit all twenty and
4.2 % admit ten or fewer. The anchor pockets are where it bites: on 1AO7 the HLA-A\*02:01 P2 pocket
admits 3 of 20 and the solvent-exposed P8 admits 19.

Building the side chain
------------------------

:mod:`tcren.refine.sidechain` places every heavy atom of any of the nineteen non-glycine types from
a backbone and a chi vector. **Chemistry is declared, geometry is measured**:
:data:`~tcren.refine.sidechain.TOPOLOGY` names the atoms, their frames and which chi drives each
torsion and carries no number; the bond lengths, angles and fixed torsion offsets are fitted from
374 Canonical2026 crystals over 176,475 residues and ship as ``data/sidechain_geometry.json``.

Round-tripping a deposited residue through its own chi reads **0.0784 A** mean heavy-atom RMSD
(median 0.0466, p95 0.2323) over 54,954 residues -- worst ARG at 0.2210 A, best SER at 0.0207 A.

.. code-block:: python

   from tcren.refine.sidechain import build_sidechain, chi_angles, rebuild_residue

   chis = chi_angles(residue)                    # what the crystal adopts
   atoms = build_sidechain(n, ca, c, cb, "W", [0.0, 0.0])
   trp = rebuild_residue(residue, "W", [0.0, 0.0])

Cost, and how it is calibrated
-------------------------------

The scan is positions x 20 candidates x ``n_grid ** n_chi`` rotamers, so it runs as **one batched
call per structure** into ``_relax.substitution_scan``, which releases the GIL and threads over
positions itself. Do not call it per residue and do not wrap it in a process pool.

Against the :func:`~tcren.energetics.background.scan_numpy` reference it agrees to **2.2e-16** on
the contact probability and exactly on the integer outputs, at a **100x** speed-up; **27.5 ms**
per structure over 60 crystals (median 23.4).

The threads run over the ``20 P`` ``(position, candidate)`` cells rather than over positions, from
an atomic cursor: the peptide side has nine positions, so a per-position split would cap the pool
at nine however many cores are free, and the cells are not equal work -- a tryptophan is nine atoms
and an alanine none. A clashed rotamer is abandoned at the first overlap, since its energy is never
weighted and its contacts never counted. Together: peptide **505 ms single-threaded to 23.1 ms on
16 cores**, receptor **6,372 to 234 ms**, both still exact against ``scan_numpy``.

The chi grid and the clash tolerance are picked against one invariant -- **a residue the crystal
deposited must be admissible at its own position**, because the crystal is the proof that it fits.
Over 400 deposited residues in 40 crystals:

.. list-table::
   :header-rows: 1
   :widths: 20 20 30 30

   * - step
     - clash_tol
     - native admissible
     - ms / structure
   * - 120 deg
     - 0.8 A
     - 0.9300
     - 4.0
   * - **60 deg**
     - **0.8 A**
     - **0.9725**
     - **23.1**
   * - 30 deg
     - 0.6 A
     - 0.9700
     - 350

Halving the step again does not help, so the residual 2.75 % is strained crystal geometry rather
than a coarse grid. ``scripts/calibrate_background_scan.py`` is that sweep.
