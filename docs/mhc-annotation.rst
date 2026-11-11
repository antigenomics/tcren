MHC annotation and supplied alleles
===================================

Build the curated sequence reference once with ``tcren build-mhc-ref`` before
annotating structures. The search index is then built as needed. The FASTA
reference itself is not automatically downloaded by scoring commands.

The standard structure score needs MHC chain roles, class and groove annotation.
Groove regions are projected from canonical class/role sequences. Exact allele
identity and the optional NetMHCpan pseudosequence annotation are not inputs to
the default cohort score. Allele identity is context for the explicit allele
confidence proxy.

A cropped MHC chain may match several reference alleles identically. The sequence
search selects a deterministic representative and retains co-best alleles among
its returned hits. This is not a genotype call or an exhaustive equivalence
class: the search returns at most 50 hits per chain. ``query_coverage`` and
``exact_observed_match`` describe the representative alignment, not full-length
allele identity. ``search_hit_limit_reached`` signals a possibly truncated list.

Supply known context in ``metadata.tsv`` beside the structures, keyed by unique
``id`` (or ``complex.id``). Use ``mhc`` for an allele such as ``HLA-A*02:01`` or
``A*02``. Class II can use ``mhc.a`` and ``mhc.b``, or a slash-separated
``DRA*01:01/DRB1*01:01`` pair in ``mhc``. Supplied field resolution is retained;
``A*02`` is never expanded into an arbitrary two-field allele.

The mapped call records ``allele_source=supplied``, the separate
``representative_allele``, and ``supplied_support``. Exact observed-sequence
support means the observed chain is contained in an existing reference sequence
consistent with the declaration. It cannot resolve differences outside that
chain. Non-exact sequence support or absence from returned co-best hits emits a
warning and remains visible in ``supplied_support``; absence from a limited hit
list alone does not establish a genotype conflict. Contradictory class/chain
roles or conflicting metadata fail explicitly. Supplied metadata does not
change the mapped role or groove projection. The energy-score table exposes MHCa evidence with the ``mhc.`` prefix.
For both chains of a class-II complex, Python callers can inspect
``structure.mhc_calls``; the MHCa table fields do not describe the MHCb allele.
