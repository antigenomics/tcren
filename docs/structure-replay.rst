Sequence-verified structure replay
==================================

``tcren.validation.replay.compare_replay`` compares Cα positions after explicit
sequence-based chain assignment. Supply literal sequences for roles A (TRA),
B (TRB), C (peptide) and D (MHC). It reports fits on all shared positions, TCR,
and MHC, measuring each component in each fit frame. MHC-fitted TCR RMSD retains
docking displacement; TCR-fitted TCR RMSD measures internal receptor agreement.
The all-position fit includes only declared shared query residues, not omitted
constant domains, beta-2 microglobulin or unresolved reference atoms.

Default model mapping rejects longer constructs. ``model_extended=True`` and
``reference_extended=True`` explicitly allow exact query spans in longer chains.
The terminal-missing budgets are stated separately. For a known substituted
reference, ``reference_spans`` can supply multiple exact segments per role with
``output_chain``, ``input_start``, ``output_start`` and ``length``. Starts are
zero-based sequence indices, not author residue numbers. The library checks
literal amino acids, overlap and unique role chains; it never chooses a
substitution alignment by structural RMSD. Omitted positions stay visible in
``spans[].input_indices`` and coverage counts.

``inspect_replay`` adds existing tcren geometry/contact diagnostics and explicit
user-declared RMSD alert thresholds. Findings require investigation rather than
silently repairing, filtering or ranking a model. An empty TCR–peptide contact
set is a finding even if global confidence is high. Contact counts cover whole
mapped output chains; exact RMSD coverage is a separate scope. Confidence must
come from the native predictor and must not be inferred from RMSD or native
experimental B factors.

The repository command ``scripts/compare_structure_replays.py`` reads a JSON
manifest with unique ``comparisons`` IDs and an ``alerts`` mapping, and writes a
fresh report with input hashes. Each comparison supplies ``id``, ``reference``,
``model``, ``sequences`` and optional ``mapping`` keyword arguments. Example alert
key ``MHC/TCR`` means fit on MHC and measure TCR. All values are in Å. Thresholds
must be predeclared for the intended diagnostic scope; there is no universal
unlike-method RMSD accuracy gate. Preserve every producer rank and distinguish
confidence-selected rank from post-hoc best agreement.

Prediction-to-prediction agreement measures reproducibility, not experimental
accuracy or affinity. Native comparisons require documented construct identity,
covered residues and template overlap. Different covered spans must not be
pooled into an apparently uniform accuracy ranking.

Geometry and ensemble review
----------------------------

The extreme interchain heavy-atom overlap guard defaults to 1 Å and records its
threshold separately. This conservative alarm catches impossible overlaps even
for a coordinate-identical replay. It is not comprehensive stereochemical
validation; ordinary steric clashes still require reviewing the raw distances.
The optional threshold must be positive and at most 2 Å, the diagnostic search
radius. Previously written reports retain their original diagnostic scope.

``replay_ensemble_consistency`` compares a small producer ensemble using pairwise
MHC-fitted TCR Cα RMSD and a mandatory declared threshold. Complete-linkage
clusters bound every within-cluster pair. It reports the unique largest cluster,
its fraction and whether the confidence-selected model belongs to it. Ties or a
selected model outside that cluster trigger review. A unique largest cluster need
not be a majority. No old or native coordinates enter clustering; agreement can
still be inaccurate. This diagnostic never replaces the producer-selected rank.

``validation.structures.same_atom_coordinates`` checks exact ordered heavy-atom
identity and unaligned coordinates while ignoring producer timestamps and other
CIF metadata. Use native ranking metadata to identify the selected sample, then
verify its atoms against the producer's copied top-model file. File-byte identity
alone can fail because those two writes may cross a timestamp boundary. Source
file hashes should still be retained for provenance.
