.. _reproducing-the-manuscript:

Reproducing the TCRen2 manuscript
=================================

Every number the TCRen2 paper prints comes out of three repositories, and the boundary between
them is deliberate:

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - repository
     - what it holds
     - located by
   * - ``tcren`` (this one)
     - the library: every descriptor, potential, score and fit
     - installed from PyPI, or the checkout you are reading
   * - ``2026-tcren2-code``
     - the benchmark: datasets, task definitions, evaluation, the ledger
     - ``$TCREN2_BENCH``, default ``~/vcs/projects/2026-tcren2-code``
   * - ``2026-tcren2-ms``
     - the LaTeX and the generated publication layer
     - ``$TCREN2_MS``, default ``~/vcs/manuscripts/2026-tcren2-ms``

**No AUC is computed here and no featurisation happens there.** This page walks the four commands
this library contributes — ``fetch-data``, ``features``, ``assess``, and the two refits that
regenerate the frozen objects — and then hands off to the benchmark's ``recompute.sh`` and the
manuscript's ``make``.

The released version
--------------------

The manuscript reports **tcren 3.3.0**. Pin it, because two boundaries in the project's history
change scores rather than only code:

.. code-block:: console

   $ pip install "tcren==3.3.0"
   $ tcren info
   tcren 3.3.0
   arda: available (.../site-packages/arda)
   default potential TCRen2: 380 pairs (karnaukhov2022: 380)

The other stamp to record is the **descriptor-catalogue digest**, which ``info`` does not print:

.. code-block:: console

   $ python -c "from tcren.provenance import registry_digest; print(registry_digest())"

That digest is what makes a feature table readable. A table is only valid under the catalogue that
wrote it, and both ``tcren recognize --features`` and ``tcren assess`` **refuse** a table stamped
with a different one, naming the command that regenerates it. The digest moved at 2.29.0 and again
at 2.30.0.

Two things not to cross without noticing:

* **The TCR:peptide potential moved twice.** ``tcren2`` replaced ``karnaukhov2022`` as the default
  at 2.11.0 (374 → 362 derivation structures) and was re-derived at 3.3.0 (362 → 359, when three
  withdrawn PDB entries were dropped). The 3.3.0 matrix agrees with its predecessor at Pearson
  *r* = 0.9996 over all 380 cells, but a Φ computed under one is still not a Φ computed under the
  other. See :doc:`potentials`.
* **The score layer was rebuilt at 3.0.0.** Everything fitted out of fold was removed at 2.28.0 and
  the five read-outs of :doc:`assess` replaced it. Nothing before 3.0.0 emits them.

Step 0 — a virtualenv and the reference data
--------------------------------------------

``bootstrap.sh`` at the repository root does this step end to end, and it is **outside any time
budget**: it downloads about 620 MB from Hugging Face and builds the IMGT/HLA MHC allele reference,
which is minutes of network and minutes of ``mmseqs``, on a machine and a connection this project
cannot predict.

.. code-block:: console

   $ ./bootstrap.sh                 # venv + editable install + Native2026/Canonical2026 + MHC ref
   $ ./bootstrap.sh --sets all      # every structure set the benchmark reads (about 620 MB)

Under it are two library commands and one environment variable.

**Where the data goes.** Since 2.12.1 every on-disk root resolves through
:func:`tcren.paths.tcren_home` — ``$TCREN_HOME``, else the source checkout recognised by its
``pyproject.toml``, else ``$XDG_CACHE_HOME/tcren``. The structure sets resolve under
``$TCREN_DATA_DIR``. Set them, or accept the cache default; do not derive a path from the source
layout, because an installed wheel then cannot find the reference at all.

**The structure sets.** ``tcren fetch-data`` reads its listing from the Hub rather than a
hard-coded tuple, so ``--list`` shows what actually exists:

.. code-block:: console

   $ tcren fetch-data                        # Native2026 + Canonical2026, the orientation refs
   $ tcren fetch-data --list                 # every set on isalgo/tcren_structures
   $ tcren fetch-data --sets all             # all of them
   $ tcren fetch-data --sets cpl,tcrvdb      # two of them

The sets the five benchmark blocks read, and which block each one serves. ``--list`` is the
authority on what exists; ``as_case`` and ``immrep23_negatives`` are also deposited and are
not in this table.

.. list-table::
   :header-rows: 1
   :widths: 26 30 12 32

   * - ``--sets`` name
     - benchmark dataset
     - structures
     - block
   * - ``Native2026``
     - Native2026
     - 374
     - the set the potential is derived from
   * - ``Canonical2026``
     - Native2026, re-oriented
     - 374
     - the default ``superimpose`` database
   * - ``cpl``
     - CPL
     - 2,103
     - 1 — peptide ranking for a fixed receptor
   * - ``tcrvdb``
     - TCRvdb
     - 618
     - 2 — receptor ranking, functionally validated
   * - ``vdjdb_binder_benchmark``
     - VDJdb-binder-benchmark
     - 1,089
     - 3 — receptor ranking under template scarcity
   * - ``atlas_ala``
     - ATLAS-affinity
     - 199
     - 4 — kinetics and the ergodicity bridge
   * - ``garcia_b27``
     - Garcia-B27
     - 7
     - 4 — peptide potency and ergodicity
   * - ``vdjdb_mut2023_benchmark``
     - VDJdb-mut2023-benchmark
     - 81
     - 3 — the mutational-scan cohort
   * - ``Bobisse`` / ``Bigot``
     - Bobisse, Bigot
     - 7 / 14
     - the two neoantigen rankings
   * - ``Native2022`` / ``PolyV2022``
     - the 2022 reproductions
     - 242 / 163
     - the Nat Comput Sci 2022 oracle

**The MHC allele reference is built, not downloaded.** It comes from IMGT/HLA plus the UniProt
mouse H-2 records, on demand:

.. code-block:: console

   $ tcren build-mhc-ref

Without it every structure fails to annotate. ``bootstrap.sh`` runs it and fails loudly if it does
not produce a reference.

Step 1 — descriptors
--------------------

One flat table, one row per structure, columns catalogued in ``tcren.recognition.DESCRIPTORS``:

.. code-block:: console

   $ tcren features -s "$TCREN_DATA_DIR/vdjdb_binder_benchmark" \
       -i placement,interface,topology,energetics \
       -t 0 -o feats.tsv

Three properties of this command matter for a reproduction and are easy to lose:

* **The whole set is annotated in one call.** One ``arda`` call per organism and one ``mmseqs`` MHC
  search for the set, never one per structure. Wrapping ``tcren features`` in a process pool
  reverses that: every worker spawns its own ``mmseqs``, which threads internally, and the run gets
  slower. Use ``-t`` to set the thread count of the two batched searches.
* **Chunk for memory, not for speed.** A parsed complex is about 1 MB, so a corpus of tens of
  thousands does not fit in one pass. ``--chunk N --resume`` splits it and keeps the batched
  annotation per chunk:

  .. code-block:: console

     $ tcren features -s big_corpus/ --chunk 2000 --resume -o feats.tsv

  The chunked table equals the whole-set table exactly; a test asserts it.
* **Every table is stamped.** ``feats.tsv`` gets a sibling ``feats.provenance.json`` carrying the
  command, the released version and a SHA-256 digest of the descriptor catalogue. That stamp is
  what lets a reader tell a table written under 3.3.0 from one written under 2.30.0.

The families are described in :doc:`features`; ``potts`` and ``kinetics`` are off by default because
they are the expensive two, and ``--all`` adds them.

Step 2 — the score set
----------------------

.. code-block:: console

   $ tcren assess --features feats.tsv -o assessed.tsv

One row per structure, and **every read-out in it is defined for a single structure** — the
transform, the class means and the covariance are frozen on a hold-out that ships in the wheel, so
nothing is estimated from the rows being scored:

``peptide_score``, ``pose_score``, ``binder_score``, the ``channel_*`` columns and the fit-free
predecessor ``S``. :doc:`assess` is the reference.

**Which columns you get depends on which families you featurised**, because a channel is a marginal
of one family's sub-block. The four default families give ``channel_placement``,
``channel_interface``, ``channel_shape`` and ``channel_energetics``; ``channel_mechanics`` needs
``-i kinetics`` at the ``features`` step. Over the 374 ``Native2026`` crystals with the four default
families, ``assess`` writes 374 rows by 11 columns.

Two things decide whether the table is right:

* ``--peptide`` when **the peptide is what varies** across the structures being compared — a
  combinatorial library or a mutational scan. It restores the five descriptors computed without the
  receptor. On a receptor benchmark those are constant across every structure of one epitope on one
  allele, so a model reading them reaches the cohort's name without reading an interface; that is
  why they are marginalized out by default. Block 1 (CPL) uses ``--peptide``; blocks 2 and 3 do not.
* **An ipTM column, if the generator wrote one.** ``assess`` picks up a column named ``iptm`` on its
  own (``--iptm-column`` renames it) and then adds ``confidence_residual``, ``binder_iptm`` and the
  AlphaFold-band diagnostic; without one it says so and skips all three. ``binder_iptm`` is
  ``binder_score`` plus ``logit(ipTM)`` — two log-odds added with no coefficient to fit — and it is
  the recommended read when a confidence is available.

Step 3 — the two frozen objects, regenerated
--------------------------------------------

Neither is taken on trust, and both are commands rather than prose.

**The frozen score model.** Its hold-out manifest ships inside the wheel;
``score.holdout_manifest()`` returns the 8,292 structures it names, with their dataset, epitope,
binder label and ipTM. The descriptors are not shipped — 8,292 rows by 147 columns is 19 MB — so
the refit is three commands:

.. code-block:: console

   $ tcren fetch-data --sets all
   $ tcren features -s <the structures the manifest names> -o hold.tsv
   $ tcren fit-holdout --features hold.tsv -o refit.npz

``refit.npz`` matches the shipped arrays to a relative 1e-5, which is the bound the test suite
asserts and the bound the optimiser earns: each of the 55 Yeo-Johnson lambdas comes from a Brent
search that stops at its own tolerance of about 1.5e-8, and every array fitted through the
transform inherits that. It is not bit-identical across platforms and the docstring says so.

**The TCR:peptide potential.**

.. code-block:: console

   $ tcren derive-potential \
       --structure-dir "$TCREN_DATA_DIR/Native2026" \
       --balance both \
       -o TCRen2_check.csv \
       --contacts-out derivation_contacts.csv

That is the whole recipe, and it is the one the released matrix was derived with;
``data/potentials.json`` records it verbatim and ``tests/regression/test_shipped_potentials.py``
reproduces the file from it. ``--contacts-out`` writes the contact table **after** every filter, so
a downstream analysis counting the structures behind the shipped matrix reads the derivation's own
set rather than assembling one: 7,816 TCR:peptide contacts over 358 structures, one fewer than the
359 the derivation keeps, because one kept structure makes no TCR:peptide contact at all. See
:doc:`potentials` for what the two filters drop and why.

Step 4 — the figures and the tables
-----------------------------------

The evaluation and the LaTeX live in the other two repositories. ``reproduce.sh`` at this
repository's root runs the library half and then hands off to them, locating each by an environment
variable:

.. code-block:: console

   $ ./reproduce.sh                       # the library half, then the benchmark, then the LaTeX
   $ TCREN2_BENCH=~/src/bench ./reproduce.sh
   $ ./reproduce.sh --library-only        # stop after this repository's own checks

What it runs there:

.. code-block:: console

   $ cd "$TCREN2_BENCH" && ./recompute.sh    # every fixture, table and .dat the manuscript reads
   $ cd "$TCREN2_MS"    && make              # figures, then main.pdf + supplementary.pdf

``recompute.sh --list`` prints its stages, and a single stage can be re-run by name. It reads the
same ``$TCREN_DATA_DIR`` this page sets, and it calls the **installed** ``tcren`` — so the version
that produced a number is whatever ``tcren info`` reported when it ran, not whatever is checked out
here.

.. note::

   **The corpus-scale featurisations do not fit in a 30-minute budget.** A single pass over the
   descriptor-harmonization corpus is more than twenty thousand structures and runs for hours. Their
   outputs are published as a derived-table deposit that ``tcren fetch-data --tables`` pulls, each
   file carrying the verbatim command, the released version and the catalogue digest that made it.
   ``reproduce.sh`` re-runs those commands live on a small deterministic slice and asserts the two
   agree, so what is fetched is precomputed rather than unreproducible. **That deposit is not yet
   published**: ``tcren fetch-data --tables`` currently fails against
   ``isalgo/tcren2_reproduce``, which does not exist on the Hub. ``reproduce.sh`` therefore does not
   call it, and no stage it runs depends on it.

What to check when a number does not match
------------------------------------------

In this order, because each one explains the ones below it:

#. ``tcren info`` — the version and the catalogue digest. A digest mismatch is the whole answer.
#. The ``.provenance.json`` beside the feature table — which command wrote it, under which version.
#. Which potential the Φ columns were computed under. ``tcren info`` names the default; a score
   computed under ``karnaukhov2022`` is not comparable with one computed under ``tcren2``, and the
   3.3.0 and 2.11.0 ``tcren2`` matrices are not comparable with each other either.
#. Whether ``--peptide`` was passed. It changes which descriptors reach the model, not how they
   are combined.
#. ``$TCREN_DATA_DIR``. If ``Canonical2026`` is absent every structure fails to orient, and a
   command that used to exit 0 with a table that simply had no energy columns now exits non-zero
   instead — that guard exists because the silent version cost a debugging session.
