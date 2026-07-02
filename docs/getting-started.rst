Getting started
===============

Installation
------------

From PyPI (binary wheels ship the C++ extensions and pull in the TCR-annotation backend):

.. code-block:: console

   $ pip install tcren

For development — the conda toolchain, an editable install, and the reference data fetched
into ``data/``:

.. code-block:: console

   $ bash setup.sh
   $ conda activate tcren

The TCR-annotation backend ``arda`` (mmseqs2-based) is a normal dependency, published to PyPI as
`arda-mapper <https://pypi.org/project/arda-mapper/>`_ (it imports as ``arda``); ``pip``/``setup.sh``
pull it in automatically. From ``arda-mapper >= 2.0.3`` it auto-fetches its own reference on first
use — no ``ARDA_HOME`` to set. ``tcren`` also builds three small pybind11/C++ kernels on install
(``tcren._align`` MHC-pseudosequence alignment, ``tcren._refine`` DOPE refinement, ``tcren._fold``
CCD loop closure).

Command line
------------

End-to-end candidate-epitope scoring from a structure:

.. code-block:: console

   $ tcren score -s complex.pdb -c candidates.txt -o ranked.csv

The full pipeline — annotate → superimpose → resmarkup / canonical Cα / contacts → per-interface
energies (TCRen for TCR↔peptide, MJ for TCR↔MHC and peptide↔MHC) plus the total — is one command
(``tcren.run_pipeline(structure)`` in the library):

.. code-block:: console

   $ tcren pipeline -s complex.pdb -o scores.csv

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
populate them. Structure outputs are plain ``.pdb`` by default — add ``--mmCIF`` for ``.cif`` and
``--compress`` for a trailing ``.gz`` (these flags apply to every command that writes a structure).

Fetch recent TCR-pMHC structures from the RCSB into ``data/pdb_recent`` (mmCIF ``.cif.gz``,
validated to have all five required chains):

.. code-block:: console

   $ tcren fetch-recent --discover --after 2024-01-01

What tcren can answer
---------------------

From a single TCR–peptide–MHC structure (crystal or model), each task is one command:

=================================================  ==========================================================
question                                           command
=================================================  ==========================================================
Which candidate epitopes does this TCR recognise?  ``tcren score -s c.pdb -c candidates.txt -o ranked.csv``
Is this peptide a strong binder for this TCR?       ``tcren rank -s c.pdb -o rank.csv``
How does a mutation change recognition (ΔΔG)?       ``tcren ddg -s c.pdb --native EPI --alanine-scan``
Is this modelled TCR a binder or a non-binder?      ``tcren binder -s model.pdb -o binder.csv``
Full three-interface energy breakdown?              ``tcren pipeline -s c.pdb -o scores.csv``
Substitute a peptide and relax its pose?            ``tcren refine -s c.pdb --substitute KQWLVWLFL -o out/``
=================================================  ==========================================================

Case studies
------------

* **Screen candidate epitopes.** ``tcren score`` ranks a candidate list by TCRen energy on the
  native contact map (no re-docking) — the drop-in for the original ``run_TCRen.R``. Add
  ``tcren rank`` to place the top hit's energy in a random-background percentile.

* **Neoantigen / alanine ΔΔG.** ``tcren ddg`` re-scores mutants on the native contacts:
  ``--alanine-scan`` for a per-position sensitivity profile, or ``--mutant`` (repeatable) for
  specific neoantigen substitutions. Positive ΔΔG = destabilising.

* **Rank candidate TCRs against a fixed pMHC.** ``tcren binder`` scores AlphaFold/TCRmodel2 models
  from interface geometry (size, dual-chain balance, H-bonds, buried ΔSASA) plus a CDR1/2-vs-CDR3α
  TCRen term — AlphaFold-orthogonal signal that ranks binders above non-binders using no external
  tool. See :func:`tcren.binder.binder_score`.

* **Substitute + refine a pose.** ``tcren refine --substitute`` threads a new equal-length peptide
  onto the backbone and runs a knowledge-based Monte-Carlo refinement scored by the DOPE atom-level
  potential — deliberately *independent* of the TCRen/MJ scoring potentials so the pose is not
  optimised against the quantity it is later scored with. This is not physics relaxation; use Rosetta
  FlexPepDock for that.

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

Iterate over a batch (file, directory, or ``.tar.gz``):

.. code-block:: python

   from tcren.structure import iter_structures

   for pdb_id, structure in iter_structures("batch.tar.gz"):
       classify_chains(structure, organism="human")
       ...

Orient into the canonical frame, layer contacts, and read the docking geometry:

.. code-block:: python

   from tcren.mhc import annotate_mhc
   from tcren.orient import canonicalize_structure, superimpose, docking_angles
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
