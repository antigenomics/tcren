Experimental ranking APIs
=========================

These explicit research commands are separate from the default :doc:`screening`
workflow. Their fitted parameters and evaluation populations must be recorded
when reproducing an experiment. They never replace a screen reference implicitly.

Exploratory contact-profile evaluation
--------------------------------------

``tcren contact-profile-evaluate --contacts contacts.parquet --atoms atoms.parquet
--baseline baseline.parquet -o evaluation -t 0`` evaluates residue involvement
from fresh ``export --quantity contact-atoms --interface tcr_peptide`` and
``export --quantity atom-inventory --interface all`` exports at 5 Angstrom.
The baseline requires unique ``id``, ``epitope``, binary ``y``, ``iptm``,
``prepared`` and ``matched`` columns. This experiment is separate from the
released/default screening score.

Flat atom-pair, residue-pair and zero-inclusive residue tables retain observed
TCR/CDR3 and peptide sequences. Counts are distinct observed heavy atoms,
including backbone; occupation divides contacted by observed available atoms.
It does not assert canonical completeness or polar chemical satisfaction.
Profiles use anchor-excluded CDR3 positions padded after residue 3 and peptide
positions, separately by chain. Missing positions have explicit indicators.
Separate alpha-only and beta-only CDR3 arms accompany the joint profiles.
The state-aware arm retains contacted atom count, any-contact and full-contact
indicators, and occupation restricted to partially contacted residues (zero
otherwise). This separates masses at zero and one without fitting a density.

Five outer folds and three inner folds keep transitive components sharing either
exact CDR3 chain together. Inner folds select ridge regularization, profile family
and nonnegative mixture weights using mean AP on the two priority epitopes named
in the experimental implementation. Imputation and scaling use training rows
only. Outputs retain all five epitope metrics, the no-contact mixture, fitted
parameters, folds and additive logit contributions. The convex mixture is bounded
in [0,1], not an experimentally calibrated binding probability. Descriptive
Mann--Whitney P-values are uncorrected and model-level; related receptors may be
dependent. These outputs are internal exploratory validation, not independent
external validation or automatic promotion of a fitted score.

Contact engagement additions
----------------------------

``tcren contact-engagement-evaluate`` compares engagement additions against a
supplied, fixed ``prepared`` ipTM/Betti baseline. It takes fresh ``--contacts``,
``--atoms`` and ``--baseline`` tables, plus corresponding ``--mock-contacts``,
``--mock-atoms`` and ``--mock-baseline`` inputs. Baselines additionally identify
``dataset`` and ``stratum``; ambiguous outcomes remain in the raw benchmark
inventory for reference exclusion but do not enter classification.

The mock reference excludes benchmark epitopes, either exact CDR3 chain and
coordinate duplicates (the caller checks coordinate hashes before export).
Reference residues receive equal-epitope/equal-CDR3 model weights. Inclusive
amino-acid-specific empirical tails retain zero contacts and use one tail
pseudocount. Zero-contact information is exactly zero. Contributions never
divide by residue atom capacity. Raw counts and endpoint profiles remain
separate candidates; signed categorical likelihood ratios preserve zero,
partial, full and missing-position states. These are composite ranking scores,
without independent-residue or calibrated-probability claims.

Five outer/three inner exact-chain-grouped folds choose regularization and
``S=(1-lambda)*prepared+lambda*contact``. ``--policy priority`` shares the
TTD/NLV-priority choice across epitopes. ``safeguard`` chooses by each epitope's
inner AP, retaining the fixed baseline as a candidate. ``ranking`` additionally
fits a smoothed AP objective anchored to the baseline log odds. Its aggregate
families include an explicit ipTM-minus-baseline contrast, allowing the
confidence/topology balance to change without treating the inputs as
independent evidence. Imputation, indicators and scaling use training rows.
``--policy stable`` retains those AP candidates but admits an addition for a
cohort only when its AP does not fall below the fixed baseline in any of the
three inner folds. This is a selection-stability rule, not a significance test
or a guarantee about outer-fold performance.
``--policy transport`` uses empirical quantiles from each fitted model's training
scores to map its contact output monotonically into the training baseline's
score units before mixture selection. Ties use midrank percentiles. This mapping
uses neither held-out labels nor a query-cohort distribution; its saved knots
are part of the fitted model. It changes score scale, not the physical
per-residue contact definition or residue atom capacity.
``--policy shared`` also considers a compact conditional model trained across
the TCRvdb training cohorts: absolute contacted-atom, engaged-residue and
full-contact totals for each chain/region, four separate type-tail sums and a
confidence contrast. It uses the same grouped splits and training-only mapping.
VDJdb never enters its training or selection.
``--policy bagged`` averages the three fitted member probabilities from the
existing inner training partitions instead of refitting one outer-training
model. Members retain their own imputation, scaling, score mapping and feature
contributions. All members exclude the outer test components; no new split or
bootstrap control is introduced.

Outputs retain all five TCRvdb comparisons, complete precision--recall curves,
fold parameters, feature logit contributions and reference-only VDJdb results
within each template stratum. Exact VDJdb chains are excluded from transfer
training. A positive AP delta does not establish precision dominance at every
recall. Repeated panel-guided development remains exploratory; no default
screening score is replaced by this command.

Exploratory additive contact search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``python -m tcren.contacts.additive --contacts contacts.parquet --atoms atoms.parquet
--topology topology.tsv --baseline baseline.parquet -o evaluation -t 0``
compares maximum TRA/TRB atom-pair counts, maximum distinct contacted TCR atoms,
distinct contacted peptide atoms, and both maxima divided by the peptide atom
count. Inputs are fresh ``export`` outputs at 5 Angstrom (``contact-atoms``,
``atom-inventory``, ``recognition-topology``); receptor contacts cover CDR1/2/3.
Zero-contact counts and ratios are zero. Betti means the existing recognition
``beta1/E`` on CDR:peptide + peptide:MHC + CDR:MHC, without internal edges.

The explicit exploratory score is ``iptm + a*betti + b*contact``. Shared and
separate nonnegative coefficients use zero and 1, 2, 5 times powers of ten
from 1e-5 to 50. Equal-epitope mean average precision selects coefficients;
separation p-values never select a model. Five folds keep receptors sharing
either exact chain together. Held-out scores and both baselines receive the
same training-only, within-epitope midrank transform before pooled evaluation.
Full-fit and held-out outputs remain separate. Two-sided Mann--Whitney tests
report equal-distribution separation with signed rank-biserial effects; these
are exploratory, not external validation or evidence of affinity calibration.

In-sample contact GLMs
~~~~~~~~~~~~~~~~~~~~~~

``python -m tcren.contacts.glm --contacts contacts.parquet --topology topology.tsv
--baseline baseline.parquet -o glm`` fits two pooled, unpenalized binomial GLMs
on all called models with finite ipTM and betti. Here ``betti`` is the previous
prepared score, not the isolated graph cycle fraction. The five predictors are
ipTM, betti and distinct contacted peptide, TRA and TRB heavy-atom counts.
The second fit divides each atom count by the total number of unique
TCR--peptide residue-contact edges; zero-contact densities are zero.
Standardization is numerical only; exported coefficients use input units.
These are deliberately in-sample fits, without held-out accuracy claims.
