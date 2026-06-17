Getting started
===============

Installation
------------

.. code-block:: console

   $ bash setup.sh
   $ conda activate tcren

``setup.sh`` creates the ``tcren`` conda environment and installs ``tcren`` in editable mode.
The TCR-annotation backend ``arda`` (mmseqs2-based) is pulled in automatically as a pinned
git dependency (tag ``2.0.1``); its C++ extension builds against the conda toolchain.

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
