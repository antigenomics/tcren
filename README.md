## TCRen TCR:pMHC scoring scheme

This repo contains scripts for computing energies of TCR:pMHC complexes, scoring TCRs and antigens based on their likelihood to recognize eachother and related pre- and post-processing software.

For more details please refer to our [biorxiv preprint](https://www.biorxiv.org/content/10.1101/2022.02.15.480516v1.abstract).

See R markdown files that can be executed in [R Studio](https://www.rstudio.com/) for running TCR specificity prediction using TCRen algorithm in ``script_TCRen_peptide_screening/`` folder. The files contain R code for running predictions, examples describing most common use cases (together with sample data stored in the same folder) and a comprehensive description of the core algorithm, TCR specificity prediction pipeline and expected results & discussion. The prediction process is not memory-intensive and should run in less than 10 minutes on a conventional laptop for provided use cases.
