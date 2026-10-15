The manuscript benchmark calculates virtual receptor poly-alanine referencing directly in bench/scripts/tcr_reference.py: Φ(real receptor)−Φ(poly-Ala receptor) on the same contact map, separately for all/CDR/CDR12/CDR3 and TCR:peptide/TCR:MHC with explicit potentials.

Existing CLI routes are not equivalent. score --delta references peptides; receptor-side alanine scans rebuild individual residues structurally, and the virtual mode does not cover this exact receptor calculation. Important manuscript calculations therefore bypass the CLI.

The adapter's score_one catches every exception→None; callers discard those rows without exact input/output coverage. Actual inspected batches retained all1,089 structures; silent loss is a source-level defect, not an observed loss in those batches.

Follow-up in tcren: one batched CLI route with precise virtual/structural semantics, physical chain/region identity, potential identities/units, complete key accounting, retained errors and nonzero all-failed status. Audit the equation before adopting any final fit-free score; do not add a fitted replacement.
