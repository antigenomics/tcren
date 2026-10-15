# Retired material and pre-audit records

Archived on 2026-10-14 at the author’s request. MANIFEST.tsv records original paths,
SHA-256 values and whether a file was moved or snapshotted. Source/resources/tests of
the geometry experiment moved together; no scientific results were regenerated.
Snapshots retain the pre-reconciliation instructions and mixed-source CLI/reproduction
context. They are historical records, not current instructions or supported entry points.

To reproduce the retired experiment, restore its complete source/resource/test closure
in an isolated checkout at the matching recorded revision; do not add this directory to
PYTHONPATH. The active contact model, geometry helpers, self-consistent fit and MD
models remain in src. The attic is excluded from distributions and ordinary tests.
