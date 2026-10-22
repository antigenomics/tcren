# Within-residue atom-pattern potential — frozen design, 21 October 2026

This authorized candidate changes the reference observation measure, preserving
original TCRen scores and the existing physical-edge query average. It uses only
native class-I CDR→peptide and curated antibody→peptide contact tables at 5 Å.
No assay labels, responses, weights, distance/angle regression or coefficient
selection enter resource construction. No default changes.

For physical contacting residue edge e, let K_e be its number of observed physical
atom pairs, and n_e(a,i,b,j) the count with directed residue/atom labels a,i,b,j.
Saved chemical equivalences relabel atoms after physical uniqueness is checked;
they do not remove observations. Each edge supplies total mass one:

    m_abij = sum_e n_e(a,i,b,j) / K_e.
    P(i,j|a,b) = m_abij / sum_ij m_abij.
    qB(i|a) = sum_bj m_abij / sum_ibj m_abij.
    qT(j|b) = sum_ai m_abij / sum_aij m_abij.
    D_abij = -ln[P(i,j|a,b) / (qB(i|a) qT(j|b))].

All masses/probabilities and D are dimensionless. Native and antibody edges are
pooled with one unit each; no source balancing or invented effective sample size.
Resource exports raw physical counts, fractional masses, source edge totals,
conditional totals and both endpoint marginals. Correlated edges are not treated
as independent samples. Binder/target direction remains explicit.

For query edge e with supported original residue contribution E0_e:

    E_pattern = E0 + sum_e (1/K_e) sum_(i,j in e) D_abij.

Only known CDR1/2/3→peptide edges receive increments. Framework and missing
original entries retain their established policy. Unseen atom cells contribute
zero and remain in K_e; unsupported counts are exported. This is abstention, not
estimated compatibility or a normalized extension to unseen cells. Lower scores
rank binders higher, unchanged from original TCRen. The additive coefficient is
fixed at one. No fitted assay weights or scoring-time fitting.

This categorical mark distribution samples an edge uniformly, then one of its
contacting atom pairs uniformly. It is not a Bernoulli occupancy model for every
available atom pair, a joint likelihood of the whole contact graph, a direct
model of solvent exposure, or physical binding free energy. It aligns reference
and query averaging; it does not claim to eliminate all residue-composition or
model-generation bias. No discovery-derived class-exclusive blacklist is used.
