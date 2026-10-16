# Cognate saved-pose validation

The follow-up reuses existing `features` and `export pmhc-cells` CLIs; no new
library module, fit, score direction, threshold or physical model is needed.
The benchmark owns raw saved predictions, exact-reference DockQ and matched R
summaries; the manuscript owns the formatted publication layer.

Fifteen saved predictions from three sequence-verified cognate complexes extend
three single-model checks. Each target supplies five correlated models, all class
I and covered templates; exact-reference template inclusion remains separate.
Overall and TCR-pMHC-specific ipTM are different saved fields. Both select an
inferior structural match for the clearest reversal, while another target has
strong confidence ordering. Offset, coverage and graph diagnostics provide
within-target structural associations, with every weaker/undefined comparison
retained. Better DockQ can coexist with greater peptide clash burden; that
quantity includes peptide:MHC as well as peptide:TCR overlaps.

These checks establish individual structural ranking failures, not a universal
forced-pose detector or the causal source of a prediction error. Historical
reference aliases, missing per-interface lineage and independent target/class/
template validation remain evidence gates. No assay population or original score
has changed. No MD files or jobs were accessed. Source details stay in the private
benchmark and manuscript records.
