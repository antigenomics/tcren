Getting started
===============

Installation
------------

From PyPI (binary wheels ship the C++ extensions and pull in the TCR-annotation backend):

.. code-block:: console

   $ pip install tcren
   $ tcren build-mhc-ref   # once: the MHC allele reference is built on demand, not bundled

For development — a repo-local ``.venv`` via `uv <https://docs.astral.sh/uv/>`_, an editable
install, and the reference data fetched into ``data/`` (no conda; needs only ``uv`` and a C++
compiler):

.. code-block:: console

   $ bash setup.sh
   $ source .venv/bin/activate

The TCR-annotation backend ``arda`` (mmseqs2-based) is a normal dependency, published to PyPI as
`arda-mapper <https://pypi.org/project/arda-mapper/>`_ (it imports as ``arda``); ``uv``/``setup.sh``
pull it in automatically. From ``arda-mapper >= 2.5.7`` it auto-fetches its own reference **and a
static mmseqs2 binary** on first use — no conda/bioconda and no ``ARDA_HOME`` to set. ``tcren`` also
builds five small pybind11/C++ kernels on install:
``tcren._align`` (MHC-pseudosequence alignment), ``tcren._refine`` (DOPE Monte-Carlo refinement),
``tcren._relax`` (DOPE interface energy), ``tcren._fold`` (CCD loop closure) and ``tcren._geom``
(binder interface geometry).

.. _mhc-reference:

The MHC allele reference — build it once
----------------------------------------

MHC chains are mapped against a curated allele reference that is **built from IMGT on demand, not
bundled in the wheel**. Run it once after installing; every command that annotates a structure
needs it:

.. code-block:: console

   $ tcren build-mhc-ref

.. _data-locations:

Where the data lives
--------------------

:func:`tcren.paths.tcren_home` is the root of tcren's on-disk reference data, and everything else
hangs off it: the allele reference in ``database/mhc/``, its mmseqs index in ``data/mhc_cache/``,
and the structure sets ``fetch-data`` populates in ``data/`` (``Native2026``, ``Canonical2026``,
``PDB_date.tsv``). It resolves in three steps:

#. ``$TCREN_HOME``, when set — the way to put the data on a shared volume;
#. otherwise the source checkout, recognised by its ``pyproject.toml``, so a development install
   reads and writes the repo's own ``data/``;
#. otherwise ``$XDG_CACHE_HOME/tcren`` (in practice ``~/.cache/tcren``), which an installed wheel
   can write and which survives an upgrade.

``$TCREN_DATA_DIR`` overrides the ``data/`` subdirectory on its own
(:func:`tcren.paths.data_dir`).

Command line
------------

End-to-end candidate-epitope scoring from a structure:

.. code-block:: console

   $ tcren score -s complex.pdb -c candidates.txt -o ranked.csv

Scoring structures
------------------

``tcren scoring`` reads structures and writes numbers: the three interface contact energies
(TCRen for TCR↔peptide, Miyazawa--Jernigan for TCR↔MHC and peptide↔MHC) and their total
:math:`\Phi`, one row per structure. It is *scoring only* — the preparation steps
(canonicalisation, region mapping, Cα / contact / atom-distance matrices) are the separate
``tcren annotate``, ``tcren superimpose`` and ``tcren contacts`` commands. In the library it is
``tcren.run_pipeline(structure)``.

.. code-block:: console

   $ tcren scoring -s complex.pdb.gz -o scores.csv

``-s`` takes a file, a directory, a ``.tar.gz``, a quoted glob, a ``.txt`` manifest with one path
per line, a comma-separated list, or a repeated flag — mix them freely:

.. code-block:: console

   $ tcren scoring -s a.pdb.gz -s b.pdb.gz -o scores.csv
   $ tcren scoring -s 'models/*.pdb.gz' -o scores.csv
   $ tcren scoring -s models.txt -o scores.csv

Two options change what is reported:

``--delta``
   adds the poly-alanine reference :math:`\Delta\Phi_I=\Phi_I(\text{peptide})-\Phi_I(\text{poly-Ala})`
   per interface, plus ``dPhi_total``. :math:`\Delta\Phi_{\mathrm{TCR:MHC}}` is identically zero — the
   peptide is not in that interface. On a *fixed* contact map :math:`\Delta\Phi` is :math:`\Phi`
   minus a constant and reorders nothing; use it when each candidate carries its **own** generated
   pose, where raw :math:`\Phi` partly reads the pose the predictor chose.

``--geometry``
   adds the interface descriptors (``burial``, ``n_pep_contacted``, ``chain_balance``, ``n_hbond``,
   ``pitch``, ``crossing``) and ``Q``, the directional decorrelated interface-quality score
   (:func:`tcren.q_score`), standardised against the native-crystal reference so it is defined for a
   single structure. For the complete descriptor catalogue plus ``P(real)``, use
   ``tcren recognize``.

Columns are named as in ``tcren recognize`` (``Phi_tcr_pep``, ``dPhi_pep_mhc``, …), but the key is
not: ``tcren scoring`` emits ``pdb.id`` and ``tcren recognize`` emits ``complex.id``, so rename
one of them before joining the two tables.

Inputs accept ``.pdb``/``.cif``/``.pdb.gz``/``.cif.gz``, a directory, or a ``.tar.gz`` batch;
identifiers are resolved from the file names:

.. code-block:: console

   $ tcren contacts -s batch.tar.gz -o contacts.csv --interface tcr_peptide
   $ tcren annotate -s complex.cif.gz -o markup.csv --regions mhc --pseudo

``tcren annotate`` emits one per-residue markup table covering TCR (CDR/FR), MHC groove
(helices/floor) and peptide; ``--regions all|tcr|mhc|peptide`` filters it to one chain class and
``--pseudo`` additionally marks the NetMHCpan MHC pseudosequence residues (region ``MPS``). It
replaces the old separate ``tcren mhc`` command.

There are two orientation commands (chains are renamed ``A``\=Vα, ``B``\=Vβ, ``C``\=peptide,
``D``\=MHCα, ``E``\=MHCβ/β2m):

* ``tcren superimpose`` brings a **new** structure into the canonical frame by superposing its
  conserved MHC groove Cα onto a canonical *database*. It detects the input's MHC class and
  species, selects every database structure of the same class and species, superposes against
  each (sequence alignment fixes the residue correspondence), and **averages** the rigid
  transforms — translations by mean, rotations by the chordal (SVD-orthonormalised) mean — into
  one consensus placement. The database defaults to ``data/Canonical2026`` (populated at install).

* ``tcren orient`` **builds** a canonical database from native complexes by deriving the
  per-class canonical frame (this is how ``Canonical2026`` itself is produced). Annotation runs
  as a single batched mmseqs2 call; ``-t`` threads only the structural alignment and write.

.. code-block:: console

   $ tcren superimpose -s complex.pdb -o oriented/
   $ tcren orient -s data/Native2026 -o data/Canonical2026 -t 8

Both need the reference sets in ``data/``; ``setup.sh`` runs ``tcren fetch-data`` at install to
populate them (the shipped database's ``orient_metadata.json`` comes with the package, so the
fetch only has to bring down structures). ``orient`` writes its own metadata to
``<out>/orient_metadata.json``, which is what ``superimpose`` reads back.
Structure outputs are plain ``.pdb`` by default — add ``--mmCIF`` for ``.cif`` and
``--compress`` for a trailing ``.gz`` (these flags apply to every command that writes a structure).

Fetch recent TCR-pMHC structures from the RCSB into ``data/pdb_recent`` (mmCIF ``.cif.gz``,
validated to have all five required chains):

.. code-block:: console

   $ tcren fetch-recent --discover --after 2024-01-01

What tcren can answer
---------------------

From a single TCR–peptide–MHC structure (crystal or model), each task is one command — or, where the
task has no command of its own, one call:

=================================================  ==========================================================
question                                           command
=================================================  ==========================================================
Which candidate epitopes does this TCR recognise?  ``tcren score -s c.pdb -c candidates.txt -o ranked.csv``
Is this peptide a strong binder for this TCR?       ``tcren rank -s c.pdb -o rank.csv``
How does a mutation change recognition (ΔΔG)?       ``tcren ddg -s c.pdb --native EPI --alanine-scan``
Which of these modelled TCRs bind?                  ``tcren features -s models/ -o f.tsv``, then ``tcren recognize --features f.tsv -o s.tsv``
Three-interface energy Φ (and ΔΦ, and Q)?           ``tcren scoring -s c.pdb --delta --geometry``
Substitute a peptide and relax its pose?            ``tcren refine -s c.pdb --substitute KQWLVWLFL -o out/``
What does a TCR meet on this pMHC surface?         ``tcren surface -s c.pdb -o surf.csv --compare``
Does the peptide hold its own conformation?        ``peptide_stability(structure).rmsf``
=================================================  ==========================================================

Case studies
------------

* **Screen candidate epitopes.** ``tcren score`` ranks a candidate list by TCRen energy on the
  native contact map (no re-docking) — the drop-in for the original ``run_TCRen.R``. Add
  ``tcren rank`` to place the top hit's energy in a random-background percentile.

* **Charge the candidate for its own conformation.** Every interface energy sums over contacts
  between two *different* chains, so a candidate held in the template's peptide conformation by its
  own side chains costs the same as one that is not. ``tcren score --intra-weight w`` adds that
  omitted term. It is sparse by design — an extended class-I 9-mer makes zero to two internal
  contacts at 5 Å with sequence separation ≥ 3 — so it separates candidates only where the peptide
  is genuinely bulged or self-packed. See :func:`tcren.intra_peptide_energy`.

* **Neoantigen / alanine ΔΔG.** ``tcren ddg`` re-scores mutants on the native contacts:
  ``--alanine-scan`` for a per-position sensitivity profile, or ``--mutant`` (repeatable) for
  specific neoantigen substitutions. Positive ΔΔG = stabilising (the mutant scores lower).

* **Rank candidate TCRs against a fixed pMHC.** ``tcren features`` then ``tcren recognize
  --features`` gives ``S``, the recommended score — three fit-free blocks over interface
  geometry, footprint topology and contact energetics, each a directional score against the
  Native2026 crystals, so it is defined for a **single structure** and its value does not depend on
  what else was scored alongside it. See :func:`tcren.reliability.s_score` and :doc:`kit`.

* **Substitute + refine a pose.** ``tcren refine --substitute`` threads a new equal-length peptide
  onto the backbone and runs a knowledge-based Monte-Carlo refinement scored by the DOPE atom-level
  potential — deliberately *independent* of the TCRen/MJ scoring potentials so the pose is not
  optimised against the quantity it is later scored with. This is not physics relaxation: it moves
  the peptide as a rigid body and leaves crystal strain, missing hydrogens and off-rotamer side
  chains elsewhere in the complex exactly where the input file put them. For that,
  ``scripts/relax_openmm.py`` minimizes every atom of the complex in amber14 with GBn2 implicit
  solvent -- the energy family an all-atom MD run evaluates, and 10-30x faster than a Rosetta
  FastRelax. It relieves the interface strain of an AlphaFold forced pose without moving the model
  off its pose, and it puts a structure in the state MD scores it in rather than the state it was
  deposited in, which matters whenever a residue-level contact score is regressed against an
  all-atom energy. Minimization, not equilibration. Needs ``openmm`` and ``pdbfixer``, which are not
  tcren dependencies; the script's docstring gives the environment and the sharding command.

* **Map the surface a TCR meets, before any TCR is there.** A contact potential scores an interface
  that already exists; ``tcren surface`` describes the pMHC beforehand, as a height field over the
  groove with hydropathy and charge painted on. The per-structure scalars turn "featureless" into a
  number — ``relief``, ``peak_to_valley`` and ``frac_above_ridge``, the fraction of peptide surface
  clearing the MHC helix crests — and ``--compare`` writes the pairwise map distance, under which
  structures of the same epitope cluster. Over the 374 Canonical2026 complexes the literature-named
  bulged epitopes rank 2nd, 5th and 8th of 230 while both named featureless ones sit at exactly
  0.000. ``notebooks/surface_topology.py`` draws all three channels; see
  :doc:`notebooks/surface_topology` and :mod:`tcren.topology.surface`.

  .. code-block:: python

     from tcren import surface_map, surface_stats, surface_distance
     smap = surface_map(structure)              # channels h / phobic / charge on a 64×32 raster
     surface_stats(smap)["frac_above_ridge"]    # how much peptide clears the helix crests
     ids, d = surface_distance([m1, m2, m3])    # pairwise map distance -> epitopes cluster

* **Ask whether the peptide holds its own conformation (backbone dynamics).** The same contact list
  comes back whether the peptide's own side chains *hold* it in the TCR-facing conformation or it
  merely happens to have been modelled there, so a static score cannot separate the two.
  :func:`tcren.peptide_stability` samples peptide φ/ψ by Metropolis Monte Carlo against DOPE and
  reports how far it wanders — ``rmsf`` (ensemble spread) and ``drift`` — rather than a better pose.
  ``intra_weight`` switches the peptide's contacts with itself on and off, and
  :func:`tcren.stability_table` runs both at one seed to give the paired ``delta_rmsf``. On the CPL
  set (2102 modelled complexes, seven clones) stability separates best from worst binders in 4/4
  clones where the additive contact energy fails and 0/3 where it works. Not MD: no solvent, no force
  field, no time, so ``rmsf`` compares between structures run at the same settings, never against an
  MD RMSF in Å. See :mod:`tcren.mechanics.dynamics`.

  .. code-block:: python

     from tcren import peptide_stability, stability_table
     peptide_stability(structure).rmsf          # ensemble spread, Å -- larger = floppier
     stability_table([s1, s2])["delta_rmsf"]    # intra-peptide term ON vs OFF, paired

* **Graft a TCR onto another pMHC (build a chimera).** ``tcren substitute-tcr`` takes a *host* and a
  *donor* TCR:pMHC complex and returns a new complex with the **host peptide + MHC** and the **donor
  TCR**. ``--by mhc`` superposes the two MHC grooves, so the donor TCR keeps its native docking
  geometry; ``--by tcr`` superposes the two TCRs, so the donor TCR inherits the host's docking pose.
  Useful for cross-docking and for building TCR:pMHC models to score. As a library call:

  .. code-block:: python

     from tcren import substitute_tcr
     from tcren.annotation import classify_chains
     from tcren.mhc import annotate_mhc
     from tcren.structure import import_structure

     host = import_structure("hostA.pdb"); classify_chains(host); annotate_mhc(host)
     donor = import_structure("donorB.pdb"); classify_chains(donor); annotate_mhc(donor)
     chimera = substitute_tcr(host, donor, by="mhc")   # host pMHC + donor TCR

* **Interface energy and koff mechanics.** ``tcren energy`` reports the DOPE interaction energy across
  the peptide↔partner interface (``e_native``; add ``--relax`` for the post-refinement energy and the
  gap — the ΔΔG scorer). ``tcren mechanics`` treats the contact map as a network of breakable springs
  and reports the stiffness tensor (``K_tens``/``K_shear``), a steered-rupture force/work, and coupling
  residues. On ATLAS these mechanical measures track the dissociation off-rate (``koff``, a Bell–Evans
  rupture quantity) better than the equilibrium ΔG/Kd — the physically apt axis for the TCR
  mechanosensor. Both are also library calls:

  .. code-block:: python

     from tcren import interface_energy, stiffness_tensor, rupture
     e = interface_energy(structure)               # DOPE interface energy (lower = more favourable)
     k = stiffness_tensor(structure)               # {"K_tens": ..., "K_shear": ..., "n_spring": ...}
     r = rupture(structure, direction="tensile")   # {"rupture_force": ..., "rupture_work": ...}

Library
-------

Score candidate epitopes against a structure:

.. code-block:: python

   from tcren import parse_structure, ContactMap, score_peptides
   from tcren.annotation import classify_chains
   from tcren.potential import tcren

   structure = parse_structure("complex.pdb.gz")     # .pdb/.cif/.pdb.gz/.cif.gz
   classify_chains(structure, organism="human")      # TRA/TRB via arda, peptide, MHC
   contact_map = ContactMap.from_structure(structure)
   ranked = score_peptides(contact_map, ["KQWLVWLFL", "RLLHPHHPL"], tcren())

Charge each candidate for the contacts it makes with itself as well (off by default):

.. code-block:: python

   from tcren import intra_peptide_energy
   from tcren.potential import mj

   contact_map = ContactMap.from_structure(structure, peptide_internal=True)
   ranked = score_peptides(contact_map, ["KQWLVWLFL", "RLLHPHHPL"], tcren(),
                           intra_weight=0.5, intra_potential=mj())
   intra_peptide_energy(contact_map, mj())          # the term alone, for the native peptide

Iterate over a batch (file, directory, or ``.tar.gz``):

.. code-block:: python

   from tcren.structure import iter_structures

   for pdb_id, structure in iter_structures("batch.tar.gz"):
       classify_chains(structure, organism="human")
       ...

Orient into the canonical frame, layer contacts, and read the docking geometry:

.. code-block:: python

   from tcren.mhc import annotate_mhc
   from tcren.docking import canonicalize_structure, superimpose, docking_angles
   from tcren.contacts import multi_contacts, ContactDefinition

   annotate_mhc(structure)
   oriented, info = canonicalize_structure(structure)   # z=MHC->TCR, y=peptide, x=thin
   oriented, info = superimpose(structure)              # onto data/Canonical2026 (class+species ensemble)
   layers = multi_contacts(structure, ContactDefinition(d1=5, d2=8, d3=12))
   angles = docking_angles(structure)                   # crossing + incident angle

Build a 2D complementarity map and summarise contacts by region pair:

.. code-block:: python

   from tcren.project2d import (project_structure, residue_markup_table,
                                contacts_table, region_pair_summary)
   from tcren.viz import render_complementarity_map

   proj = project_structure(structure)
   svg = render_complementarity_map(residue_markup_table(structure, proj),
                                    contacts=contacts_table(structure, threshold=5.0))
   summary = region_pair_summary(structure, kind="closest")   # also "cb" (8 A) / "ca" (12 A)

Data
----

Structures come from the Hugging Face dataset
`isalgo/tcren_structures <https://huggingface.co/datasets/isalgo/tcren_structures>`_ (all gzipped):
``Native2022`` (the 2022 paper set), ``Native2026`` (the 2026 set the current potential is derived
from), and ``Canonical2026`` (``Native2026`` re-oriented). When orienting a new complex an installed
library lazily fetches the canonical reference structures (1ao7/1fyt) from the Hub, so no local
dataset is required.
