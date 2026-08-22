# TCRen2_potential.csv — do not try to improve it

The matrix is a plain Boltzmann inversion of contact counts over the 374 `Native2026` crystals: a
5 Å heavy-atom cutoff, `--balance both`, a flat pseudocount of 1, and nothing else. Between
2026-08-21 and 2026-08-22 we tried five ways to estimate it better — a 6 Å and a per-structure
resolution-aware cutoff, DFIRE distance and orientation corrections, BLOSUM62 substitution
pseudocounts, BLOSUM62 nearest-donor imputation of under-observed cells, and composition-weighted
pseudocounts drawn from the actual CDR3 and peptide residue frequencies — and measured each against
all three endpoints that read the energy: CPL activation ranking (7 clones, 2,103 models, paired per
clone), TCRvdb receptor ranking against DESeq2 log fold change (n = 423, paired bootstrap), and the
Yang/Garcia B\*27:05 potency series (16 peptides).

**Nothing improved all three.** The only contrast that cleared its confidence interval was
nearest-donor imputation on TCRvdb (ΔSpearman +0.063 [+0.022, +0.109]), and it cost the potency
series its whole signal (+0.706 → +0.137); on CPL every variant was inside noise once paired
(Wilcoxon p ≥ 0.22). Removing the pseudocount entirely is the one change that clearly *breaks*
things: 17 of 380 cells go to +∞, including Ile:Ile, and 317 of 2,103 CPL models return an
undefined energy. If you are about to re-run one of these, read `../../../OBSOLETE.md` first and
then don't.
