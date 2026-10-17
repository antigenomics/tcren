# Separate contact graphs and intrinsic backbone geometry — 17 October 2026

Author-requested Figure 1/descriptor clarification reuses ContactMap, SciPy's sparse
connected-components kernel and existing C-alpha angle/dihedral primitives. No new
potential, fitted combination, threshold, binding label or MD operation.

- `export --quantity contact-topology`: CDR:peptide, CDR:MHC and peptide:MHC,
  unique physical residue edges, contacted vertices, graph beta0/beta1 and explicit
  empty-interface availability. Existing pooled graph and spectral readouts unchanged.
- `export --quantity backbone`: peptide and each TRA/TRB CDR1/2/3, C-alpha bonds,
  angles/torsions with physical residue windows. Apo alpha-beta receptors supported.
- Effective rank is the normalized order-2 Hill rank of a Gaussian proximity matrix,
  sigma=5 Angstrom, q the smaller resolved residue count. Graph cycles, spatial
  flag-complex holes, docking angles, loop apex and peptide bulge are distinct.

Verification: exact installed noneditable wheel; all nine export tests pass,
including a known cycle/disconnected/empty graph and backbone rigid rotation/translation.
Canonical-frame test passes with the local runtime reference root: ten focused checks.
Ruff `src tests` and Sphinx `-W` pass. The broad fast run was stopped at 298 seconds:
710 passed, four skipped, one missing-default-reference failure; that failure passes
when TCREN_HOME points at the checkout's populated references. This is not a full-suite pass.
Initial source-only testing could not load native extensions; final checks use the wheel.

Fresh benchmark audit: `results/cpl-shape-audit-20261017-c` in the owning benchmark.
All seven native CPL matrices are computed with both potential generations, and
seven canonical contact graphs plus thirteen backbone inputs use batched CLI exports.
Earlier partial diagnostic runs are rejected, not resumed or used as feature inputs.
Native sparse graph/table operations suffice; no worker pool or new dependency.
CPU budget is one per scientific CLI, compiler budget two. Exact bound 868 construct
matching and independent nonbinding/deformation truth remain gates, not closed claims.
