tcren package
=============

The package is laid out in layers: what a structure *is* (parsing, annotation, contacts), what can
be *measured* on it (topology, energetics, mechanics, docking geometry), the *catalogue* that names
every measurement, and the *scores* built on top. Each layer only reaches downwards.

Three layers are documented in prose elsewhere and are not repeated here: :doc:`assess` for
:mod:`tcren.score`, :doc:`reliability` for :mod:`tcren.reliability`, and :doc:`potts` for
:mod:`tcren.potts`.

.. note::

   Nine top-level modules are **deprecated locations**, kept so existing imports keep working:
   ``tcren.ddg``, ``tcren.dynamics``, ``tcren.footprint``, ``tcren.interface_graph``,
   ``tcren.pose``, ``tcren.rotamers``, ``tcren.scoring``, ``tcren.stability`` and ``tcren.surface``,
   as is the whole ``tcren.orient`` package. Each re-exports its new home, and the new home is what
   is documented below. Import the canonical name in new code.

Structure I/O
-------------

tcren.structure.model module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.structure.model
   :members:
   :undoc-members:
   :show-inheritance:

tcren.structure.io module
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.structure.io
   :members:
   :undoc-members:
   :show-inheritance:

Annotation
----------

tcren.annotation.arda_adapter module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.annotation.arda_adapter
   :members:
   :undoc-members:
   :show-inheritance:

tcren.annotation.chains module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.annotation.chains
   :members:
   :undoc-members:
   :show-inheritance:

tcren.annotation.cgene module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.annotation.cgene
   :members:
   :undoc-members:
   :show-inheritance:

tcren.annotation.batch module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.annotation.batch
   :members:
   :undoc-members:
   :show-inheritance:

MHC mapping
-----------

tcren.mhc.imgt module
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mhc.imgt
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mhc.reference module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mhc.reference
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mhc.mapper module
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mhc.mapper
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mhc.domains module
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mhc.domains
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mhc.regions module
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mhc.regions
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mhc.linker module
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mhc.linker
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mhc.pseudo module
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mhc.pseudo
   :members:
   :undoc-members:
   :show-inheritance:

Contacts and per-residue geometry
---------------------------------

tcren.contacts.geometry module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.contacts.geometry
   :members:
   :undoc-members:
   :show-inheritance:

tcren.contacts.definitions module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.contacts.definitions
   :members:
   :undoc-members:
   :show-inheritance:

tcren.contacts.table module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.contacts.table
   :members:
   :undoc-members:
   :show-inheritance:

tcren.contactmap module
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.contactmap
   :members:
   :undoc-members:
   :show-inheritance:

tcren.contact_types module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.contact_types
   :members:
   :undoc-members:
   :show-inheritance:

tcren.stacking module
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.stacking
   :members:
   :undoc-members:
   :show-inheritance:

tcren.torsions module
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.torsions
   :members:
   :undoc-members:
   :show-inheritance:

tcren.geometry module
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.geometry
   :members:
   :undoc-members:
   :show-inheritance:

tcren.clashes module
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.clashes
   :members:
   :undoc-members:
   :show-inheritance:

The descriptor catalogue
------------------------

tcren.descriptors package
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.descriptors
   :members:
   :undoc-members:
   :show-inheritance:

tcren.recognition module
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.recognition
   :members:
   :undoc-members:
   :show-inheritance:

tcren.descriptors.catalogue module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.descriptors.catalogue
   :members:
   :undoc-members:
   :show-inheritance:

tcren.descriptors.compute module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.descriptors.compute
   :members:
   :undoc-members:
   :show-inheritance:

tcren.descriptors.table module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.descriptors.table
   :members:
   :undoc-members:
   :show-inheritance:

Topology: the shape of the contact set
--------------------------------------

tcren.topology package
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.topology
   :members:
   :undoc-members:
   :show-inheritance:

tcren.topology.footprint module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.topology.footprint
   :members:
   :undoc-members:
   :show-inheritance:

tcren.topology.graph module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.topology.graph
   :members:
   :undoc-members:
   :show-inheritance:

tcren.topology.surface module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.topology.surface
   :members:
   :undoc-members:
   :show-inheritance:

tcren.topology.literature module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.topology.literature
   :members:
   :undoc-members:
   :show-inheritance:

tcren.topology.pose module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.topology.pose
   :members:
   :undoc-members:
   :show-inheritance:

Energetics: sums of a pair potential, and differences of them
-------------------------------------------------------------

tcren.energetics package
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.energetics
   :members:
   :undoc-members:
   :show-inheritance:

tcren.energetics.scoring module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.energetics.scoring
   :members:
   :undoc-members:
   :show-inheritance:

tcren.energetics.mutation module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.energetics.mutation
   :members:
   :undoc-members:
   :show-inheritance:

tcren.energetics.rotamers module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.energetics.rotamers
   :members:
   :undoc-members:
   :show-inheritance:

Mechanics: the contact map as something that can break
------------------------------------------------------

tcren.mechanics package
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mechanics
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mechanics.springs module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mechanics.springs
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mechanics.stability module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mechanics.stability
   :members:
   :undoc-members:
   :show-inheritance:

tcren.mechanics.dynamics module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.mechanics.dynamics
   :members:
   :undoc-members:
   :show-inheritance:

Docking geometry and canonical orientation
------------------------------------------

tcren.docking package
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking
   :no-index:
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.angles module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.angles
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.tcrdock_geometry module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.tcrdock_geometry
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.frame module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.frame
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.align module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.align
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.superimpose module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.superimpose
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.pipeline module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.pipeline
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.chains module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.chains
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.graft module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.graft
   :members:
   :undoc-members:
   :show-inheritance:

tcren.docking.exceptions module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.docking.exceptions
   :members:
   :undoc-members:
   :show-inheritance:

Potentials
----------

tcren.potential.model module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.potential.model
   :members:
   :undoc-members:
   :show-inheritance:

tcren.potential.derive module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.potential.derive
   :members:
   :undoc-members:
   :show-inheritance:

tcren.potential.redundancy module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.potential.redundancy
   :members:
   :undoc-members:
   :show-inheritance:

tcren.potential.aaindex module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.potential.aaindex
   :members:
   :undoc-members:
   :show-inheritance:

The score set
-------------

The five read-outs and what each answers are in :doc:`assess`; the machinery
behind them is here.

tcren.score.transform module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.score.transform
   :members:
   :undoc-members:
   :show-inheritance:

tcren.score.model module
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.score.model
   :members:
   :undoc-members:
   :show-inheritance:

tcren.score.fit module
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.score.fit
   :members:
   :undoc-members:
   :show-inheritance:

Fit-free cohort scores
----------------------

tcren.cohort module
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.cohort
   :members:
   :undoc-members:
   :show-inheritance:

tcren.binder package
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.binder
   :members:
   :undoc-members:
   :show-inheritance:

tcren.binder.noise module
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.binder.noise
   :members:
   :undoc-members:
   :show-inheritance:

Epitope scoring and ranking
---------------------------

tcren.scoring_rank module
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.scoring_rank
   :members:
   :undoc-members:
   :show-inheritance:

tcren.cpl module
~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.cpl
   :members:
   :undoc-members:
   :show-inheritance:

tcren.shuffle module
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.shuffle
   :members:
   :undoc-members:
   :show-inheritance:

tcren.pipeline module
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.pipeline
   :members:
   :undoc-members:
   :show-inheritance:

tcren.oracle module
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.oracle
   :members:
   :undoc-members:
   :show-inheritance:

Peptide substitution and refinement
-----------------------------------

tcren.refine package
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.substitute module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.substitute
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.register module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.register
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.anchors module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.anchors
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.rmsd module
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.rmsd
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.interface module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.interface
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.model module
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.model
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.oracle_flexpep module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.oracle_flexpep
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.engines.base module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.engines.base
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.engines.dope module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.engines.dope
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.engines.ccd module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.engines.ccd
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.engines.openmm_engine module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.engines.openmm_engine
   :members:
   :undoc-members:
   :show-inheritance:

tcren.refine.engines.promod3_engine module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.refine.engines.promod3_engine
   :members:
   :undoc-members:
   :show-inheritance:

Data paths, metadata and provenance
-----------------------------------

tcren.paths module
~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.paths
   :members:
   :undoc-members:
   :show-inheritance:

tcren.metadata module
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.metadata
   :members:
   :undoc-members:
   :show-inheritance:

tcren.provenance module
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.provenance
   :members:
   :undoc-members:
   :show-inheritance:

Analysis
--------

tcren.analysis module
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.analysis
   :members:
   :undoc-members:
   :show-inheritance:

2D projection and visualization
-------------------------------

tcren.project2d.frame module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.project2d.frame
   :members:
   :undoc-members:
   :show-inheritance:

tcren.project2d.tables module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.project2d.tables
   :members:
   :undoc-members:
   :show-inheritance:

tcren.project2d.pockets module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.project2d.pockets
   :members:
   :undoc-members:
   :show-inheritance:

tcren.viz.svg2d module
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.viz.svg2d
   :members:
   :undoc-members:
   :show-inheritance:

tcren.viz.surface2d module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.viz.surface2d
   :members:
   :undoc-members:
   :show-inheritance:

tcren.viz.pocket3d module
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.viz.pocket3d
   :members:
   :undoc-members:
   :show-inheritance:

tcren.viz.pymol module
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.viz.pymol
   :members:
   :undoc-members:
   :show-inheritance:

tcren.viz.palette module
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.viz.palette
   :members:
   :undoc-members:
   :show-inheritance:

Reference data and reproduction
-------------------------------

tcren.recent module
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.recent
   :members:
   :undoc-members:
   :show-inheritance:

tcren.paper.bootstrap module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.paper.bootstrap
   :members:
   :undoc-members:
   :show-inheritance:

tcren.paper.helpers module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.paper.helpers
   :members:
   :undoc-members:
   :show-inheritance:

Command line
------------

tcren.cli module
~~~~~~~~~~~~~~~~~~

.. automodule:: tcren.cli
   :members:
   :undoc-members:
   :show-inheritance:
