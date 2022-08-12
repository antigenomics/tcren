# TCRen

TCRen is a method for prediction of TCR recognition of unseen epitopes based on residue-level pairwise statistical potential

TCRen method starts from a structure of the peptide-MHC complex with the TCR of interest—either experimentally derived or based on a homology model—then extracts a TCR-peptide contact map and estimates the TCR-peptide energy of interaction for all candidate epitopes using TCRen potential, which we derived from statistical analysis of amino acid contact preferences in TCR-peptide-MHC crystal structures.

![preview](https://github.com/vadim-karnaukhov/tcren-new/blob/master/figures/Fig1.png)

## Dependencies
R (tested on R version 4.2.0)

R packages data.table, tidyverse, optparse, stringr, magrittr

Java (tested on openjdk 11.0.16)

## TCRen input 
1. A structure of TCR-peptide-MHC complex (either experimentally derived or a homology model). Several structures may be submitted at once. Structure(s) should be placed in a single folder.

* Example: ```example/input_structures```

2. A list of candidate epitopes.

* Example: ```example/candidate_epitopes.txt```

## TCRen output
A table with 4 columns: complex.id (corresponding to the name of an input structure), peptide (corresponding to the name of a candidate peptide), potential (“TCRen” if the default ```TCRen_potential.csv``` file is used) and score (TCRen estimate of energy of peptide-TCR interaction).

* Example: ```example/output_TCRen/candidate_epitopes_TCRen.csv```

## Tutorial
1. Clone the github repository for TCRen:

```$ git clone https://github.com/antigenomics/tcren-ms.git```

2. Prepare a structure (or several structures) of TCR-peptide-MHC complex. Format: “.pdb”.

* All structures for which predictions will be done should be placed in a single directory (e.g. ```example/input_structures```)

* If for the TCR of interest a crystal structure of the ternary complex (TCR-peptide-MHC) with some peptide is available (i.g. for the task of prediction of cross-reactivity of a well-known TCR), it can be downloaded directly from [PDB](https://www.rcsb.org/) and used as input. For the task of predictions for unseen TCRs, homology modeling should be used, e.g. TCRpMHCmodels tool which is implemented both as a [webserver](https://services.healthtech.dtu.dk/service.php?TCRpMHCmodels-1.0) and a [stand-alone software](https://services.healthtech.dtu.dk/cgi-bin/sw_request). Detailed instructions for TCRpMHCmodels tool use can be found on the webserver site and in README of software download. The modeling usually requires a few minutes for a single TCR-peptide-MHC complex either in the web-server or in a single CPU. 

3. Prepare a list of candidate epitopes (e.g. mutated peptides predicted as binders to host MHC for the task of prediction of neoepitope recognition). Format: “.txt” file with each candidate epitope in a separate line, with a header “peptide”. 

* Example in ```example/candidate_epitopes.csv```


4. Run TCRen pipeline. 

* It should be done from the directory containing files necessary for TCRen launching: ```run_TCRen.R```, ```TCRen_potential.csv``` and ```mir-1.0-SNAPSHOT.jar``` (java script for annotation and extraction of contacts in TCR-peptide-MHC structures)

```$ cd TCRen_pipeline```

```$ Rscript --vanilla run_TCRen.R -s ../example/input_structures/ -c ../example/candidate_epitopes.tsv -o ../example/output_TCRen/ -m 1G```

* Arguments of run_TCRen.R script:

| option                                         | Description                                      |    
|------------------------------------------------|--------------------------------------------------|
| -s, --structures                               | directory with input structures                  |
| -c, --candidates                               | file with candidate epitopes                     |
| -o, --out                                      | directory for TCRen output                       |
| -m, --memory                                   | memory allocation                                |

5. The output of TCRen can be found in the file ```candidate_epitopes_TCRen.csv``` in the folder which was set by -o flag. The content of the output file is described above in the section “TCRen output”

