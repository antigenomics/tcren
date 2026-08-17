PART1

Technical suggestions 

tcren.orient.tcr.dock_geometry and the subsequent functions are not supported for MHCII. The fix would allow the complete analysis of INSR reactive clones. 

classify_contacts rarely finds contact types besides “other”. Maybe some new approach is needed. 
Conceptual suggestions

Rotamers of the generated structure are not always in optimal positions which creates falls contacts in the modeled structure and misses the real ones. Maybe scanning through the possible rotation options and taking a weighted average of all possible potentials will result in better performance. FoldX software is actually optimizing them very fast. 


Each contact geometry and overall contact type is not taken into account while scoring. And the function to get contact types is rarely to classify them. For example stacking contacts needs a particular geometrical position of appropriate residues to be formed. 

The TCRen matrix could be recalculated with not just contacting residues but contacting residues + type of contact between them. This information may already be indirectly incorporated into the final score; nevertheless, considering it could add a degree of precision. Also - adding the contact type might lead to the increased sparsity of the data. 

At least annotation of contact types for TCRen matrix allows us to discard unexisting contacts if those were included in the initial score based entirely in proximity.


Contacts in the center of peptide might be more important when on the edges. I think for TCR to be correctly placed on pMHC the central contacts at least should not create the sterical difficulties so high TCRen potentials for this part of peptide should automatically raise the concerns about TCR:pMHC complex viability. Additionally - comparatively loose edge residues of the peptide in the complex of MHCII might not have the same impact as the center of the whole complex. Maybe we should try to include this logic to the scoring? 

structures for part 1 for whitepaper:
/projects/structures/peptide_4c6_pdb/all_pdb/
their metadata:
/projects/structures/peptide_4c6_pdb/all_pdb/meta.tsv
accessible via aldan3 client - check docs for this util to use it properly

PART2
We would also like to provide a fast way to ship pMHC with docked peptides. There is flexpepdoc, but we want our own C++ module to quickly place a peptide using tcren. We can test generated (relaxed if needed) structures to results by flexpepdoc on aldan3 and compare to native structures.

Why would we want to do it - we want to estimate paratope complementarity in general of a given pMHC - how likely would a T-cell recognize it - immunogenicity in other words. Plain comparison of RMSD doesnt account for side chain layout and properties (e.g. hydrophobicity etc). We want to compute something like https://www.pnas.org/doi/10.1073/pnas.2504783122 https://arxiv.org/html/2407.06703v3 or protein surface topology /Users/mikesh/vcs/manuscripts/2026-tcren/surface_topology and a CLI extension for tcrnet.
