# Contact chemistry, empirical potentials and AlphaFold confidence

Purpose: source-grounded interpretation of the contact-penalty investigation,
18 October 2026. Read-only literature research; no cohort analysis, native-contact
re-extraction, potential fitting or scientific rerun was performed. This note
does not establish the cause of the recorded GLC/YLQ changes or CPL ipTM extremes.

## Evidence and access ledger

Identifiers below were returned by the PubMed connector, rather than reconstructed
from memory. Descriptive labels deliberately avoid unverified bibliography fields.
“Complete text” means the article body was actually retrieved and read; an available
PMC link alone does not establish full-text inspection.

| Source | Verified identifiers | Access in this investigation |
|---|---|---|
| Arpeggio interaction classification | [PMID 27964945](https://pubmed.ncbi.nlm.nih.gov/27964945/); [DOI 10.1016/j.jmb.2016.12.004](https://doi.org/10.1016/j.jmb.2016.12.004); [PMC5282402](https://pmc.ncbi.nlm.nih.gov/articles/PMC5282402/) | Complete body, 10844 characters; rights-eligible connector retrieval |
| MolProbity all-atom validation | [PMID 20057044](https://pubmed.ncbi.nlm.nih.gov/20057044/); [DOI 10.1107/S0907444909042073](https://doi.org/10.1107/S0907444909042073); [PMC2803126](https://pmc.ncbi.nlm.nih.gov/articles/PMC2803126/) | Metadata and abstract only |
| Self-consistent knowledge-based protein–ligand potential analysis | [PMID 10794420](https://pubmed.ncbi.nlm.nih.gov/10794420/); [DOI 10.1110/ps.9.4.765](https://doi.org/10.1110/ps.9.4.765); [PMC2144609](https://pmc.ncbi.nlm.nih.gov/articles/PMC2144609/) | Metadata and abstract; full-text connector failed and web endpoint returned a browser challenge |
| Contact-definition comparison for fold recognition | [PMID 23140471](https://pubmed.ncbi.nlm.nih.gov/23140471/); [DOI 10.1186/1471-2105-13-292](https://doi.org/10.1186/1471-2105-13-292); [PMC3534397](https://pmc.ncbi.nlm.nih.gov/articles/PMC3534397/) | Metadata and abstract only |
| actifpTM | [PMID 40080667](https://pubmed.ncbi.nlm.nih.gov/40080667/); [DOI 10.1093/bioinformatics/btaf107](https://doi.org/10.1093/bioinformatics/btaf107); [PMC11925850](https://pmc.ncbi.nlm.nih.gov/articles/PMC11925850/) | Complete body, 18700 characters; rights-eligible connector retrieval |
| ipSAE | [PMID 39990437](https://pubmed.ncbi.nlm.nih.gov/39990437/); [DOI 10.1101/2025.02.10.637595](https://doi.org/10.1101/2025.02.10.637595); [PMC11844409](https://pmc.ncbi.nlm.nih.gov/articles/PMC11844409/) | Metadata and abstract only; preprint |
| NetTCR-struc | [PMID 40746534](https://pubmed.ncbi.nlm.nih.gov/40746534/); [DOI 10.3389/fimmu.2025.1616328](https://doi.org/10.3389/fimmu.2025.1616328); [PMC12311634](https://pmc.ncbi.nlm.nih.gov/articles/PMC12311634/) | Complete body, 50819 characters, retrieved at offsets 0 and 50000 |

## Proximity, chemical interaction and implausibility are different quantities

**Primary-source claims — Arpeggio.** The program enumerates atom pairs within
5 Å, then assigns mutually exclusive steric/covalent/van der Waals/proximal
classes and additional chemical interaction features. Proximity alone need not
represent a meaningful specific interaction. Its angle-free polar classification
differs from hydrogen-bond assignment. Atom typing and missing atoms can affect
classification. Source: [PMID 27964945](https://pubmed.ncbi.nlm.nih.gov/27964945/).

**Abstract-only corroboration — MolProbity.** Validation uses optimized hydrogen
placement, all-atom contact analysis, covalent geometry and torsions. Experimental
structures can contain local errors, including at high resolution. Source:
[PMID 20057044](https://pubmed.ncbi.nlm.nih.gov/20057044/).

**Implications/inference.** Do not call an untyped proximal contact “virtual” or
unnatural solely because it lacks a hydrogen-bond, ionic or aromatic label. It may
be ordinary packing. Conversely, satisfying a chemical feature rule does not
certify a correct docking model or experimentally detectable binding. Missing
antecedent atoms make angular geometry unassessable, rather than failed. Assay
negativity and rarity are also distinct from structural impossibility. A native
observation is not automatically an error-free geometry reference.

The existing [contact-prior dossier](CONTACT_PRIOR_SOURCES.md) separately records
verified PICCOLO rules and antibody-source access. Its PICCOLO evidence uses
donor/acceptor identity, hydrogen placement and angular criteria for hydrogen
bonds; those requirements cannot be replaced by heavy-atom distance alone.
PICCOLO: [PMID 21801404](https://pubmed.ncbi.nlm.nih.gov/21801404/),
[DOI 10.1186/1471-2105-12-313](https://doi.org/10.1186/1471-2105-12-313).
These are previously verified dossier claims, not a fresh full-text retrieval in
this investigation. The dossier also verifies that the deposited antibody survey
uses a 5 Å heavy-atom observation window; its endpoint coordinates do not supply
every antecedent needed for hydrogen-bond angles. Antibody survey:
[PMID 38161735](https://pubmed.ncbi.nlm.nih.gov/38161735/),
[DOI 10.1016/j.csbj.2023.11.056](https://doi.org/10.1016/j.csbj.2023.11.056),
[PMC10755492](https://pmc.ncbi.nlm.nih.gov/articles/PMC10755492/).
No numerical trRosetta prior or geometric-validity distribution was inspected
here, so none is inferred from its use of orientation variables.

## Frozen empirical tables and filtering

**Abstract-only primary-source claims.** The self-consistent potential study
identifies nonenergetic contact-frequency effects and chemical correlations as
important reference-state dependencies: [PMID 10794420](https://pubmed.ncbi.nlm.nih.gov/10794420/).
The fold-recognition study evaluates 45 contact definitions and explains that
some useful coarse contacts capture spatial arrangement without direct physical
contact: [PMID 23140471](https://pubmed.ncbi.nlm.nih.gov/23140471/).
Neither study validates a particular TCRen filter or its predictive effect.

**Exact conditional-probability identity, derived here.** Suppose an original
dimensionless potential has the specified log-ratio interpretation

\[
V(a,b)=-\ln\frac{p(a,b\mid C)}{q(a,b\mid C)}.
\]

Here \(a,b\) are endpoint residue types; \(C\) is the original contact event;
\(p\) is the observed distribution and \(q\) its reference distribution. For a
new filter event \(F\), the corresponding conditional log-ratio is

\[
V_F(a,b)=V(a,b)
-\ln\frac{p(F\mid a,b,C)}{q(F\mid a,b,C)}
+\ln\frac{p(F\mid C)}{q(F\mid C)}.
\]

The identity requires positive probabilities for the displayed ratios and
does not assert that every existing TCRen table uses exactly these distributions.
Its final term is constant across residue types but contributes once per retained
contact when potentials are summed; varying contact counts can therefore matter.

**Implications/inference.** Leaving the frozen matrix unchanged and filtering or
weighting its contact contributions preserves its numerical entries and baseline
provenance. It creates a new scoring operation. It does not automatically produce
a potential estimated for the filtered ensemble. A conditional-potential claim
would require matching observed/reference distributions under the new definition.
This distinction is not authorization or a requirement to refit the baseline.
Potential-entry signs are relative to the declared reference, rather than direct
measurements that each atom pair attracts or repels.

## What the penalty/sign-flip result establishes

**Mathematical inference, not a literature finding.** With lower-is-better energy,
\(E'=E+P\), \(P\ge0\), raises every individual energy but need not improve rankings.
It worsens a binder–nonbinder comparison if its relative increment moves that
comparison in the wrong direction. Literal sign flipping additionally turns
originally positive flagged entries favourable, whereas a nonnegative penalty
cannot do so. Accordingly, their different results do not identify a contradiction.

**Mechanistic hypotheses requiring the actual recorded contact evidence.** A
flag may select ordinary packing, preferentially penalize more extensive
interfaces, or alter an empirical energy component whose ordering is inverted in
the predicted ensemble. A sign-flip gain alone does not establish repulsive
chemistry, contact unnaturalness or any of these mechanisms. Relevant existing
quantities are original flagged contribution signs/magnitudes, retained physical
residue-edge counts and class-specific changes. Atom-pair multiplicity must not
silently duplicate a residue-level potential.

## AlphaFold confidence and local recognition

**Primary-source claims — actifpTM.** Ordinary ipTM incorporates interchain residue
pairs regardless of their participation in the interface. Noninteracting flexible
flanks can lower its value while the modeled interface remains similar. The paper
also motivates chain-pair confidence in multichain complexes: confidence of an
internal antibody interface need not characterize antigen engagement. Source:
[PMID 40080667](https://pubmed.ncbi.nlm.nih.gov/40080667/).

**Abstract-only primary-source claim — ipSAE preprint.** Whole-chain/accessory-region
effects motivate restricting residue pairs and changing TM-score length
normalization: [PMID 39990437](https://pubmed.ncbi.nlm.nih.gov/39990437/).
This is not a TCR-specific validation.

**Primary-source claims — NetTCR-struc.** AF-M can produce homogeneous incorrect
CDR configurations that lose native contacts. Conserved overall TCR:pMHC docking
geometry can yield favourable global docking assessment while important local
contacts remain inaccurate. The study generated negative examples by mismatching
TCRs and pMHCs; these are not the current assay-negative population. Source:
[PMID 40746534](https://pubmed.ncbi.nlm.nih.gov/40746534/).

**Implications/inference.** A model of an assay-negative pair can contain a
well-packed, chemically plausible interface. Model confidence and contact typing
alone do not establish binding specificity. For TCR:pMHC assemblies, distinguish
confidence involving TCRα:TCRβ and peptide:MHC from the local peptide:CDR recognition
interface; this is an analogy from the multichain-confidence result, not a direct
TCR experiment in the actifpTM paper. ipTM is not a simple atom-contact count.

## Limitations and unresolved attribution

Literature establishes these distinctions but does not explain the observed
GLC/YLQ score changes or identify why particular CPL models have extreme ipTM.
That attribution needs the actual matched models, chain definitions and recorded
contact-level contributions. Construct length/flexibility, alternative docking
confidence and local recognition geometry are candidate explanations, not findings
about those cohorts. External studies use different observation definitions,
targets and negative-generation procedures; their thresholds and performance
estimates have not been transferred. No new numeric geometry prior, smoothing
floor, penalty scale or classification threshold is supplied by this note.
