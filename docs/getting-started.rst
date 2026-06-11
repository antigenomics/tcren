Getting started
===============

Installation
------------

.. code-block:: console

   $ bash setup.sh
   $ conda activate tcren

``setup.sh`` creates the ``tcren`` conda environment, installs the sibling ``arda``
package (TCR annotation backend), and installs ``tcren`` in editable mode.

Command line
------------

The ``tcren score`` command is an end-to-end replacement for the legacy
``run_TCRen.R``:

.. code-block:: console

   $ tcren score -s example/input_structures -c example/candidate_epitopes.txt -o out.csv

Library
-------

.. code-block:: python

   from tcren import parse_structure, ContactMap, score_peptides
   from tcren.annotation import classify_chains
   from tcren.potential import tcren

   structure = parse_structure("complex.pdb")
   classify_chains(structure, organism="human")
   contact_map = ContactMap.from_structure(structure)
   ranked = score_peptides(contact_map, ["KQWLVWLFL"], tcren())
