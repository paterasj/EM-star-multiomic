# BRCA Subtype Differentially Expressed Gene Analysis & Pathway Enrichment 
### Conduct differentially expressed gene analysis based on EM* subtyping for BRCA Patients 


## Usage Instructions 

For now, please download brca_subtype_de_analysis.R and run using the brca_normal.rda and brca_tumor.rda files in this repository.
Please change parameters according to preference (i.e how many subtypes are present from the EM* subtyping, threshold for up/down regulated, etc.)

Further work remains on generalize this code and releasing it as a Bioconductor package.

## Output

The following ouputs are expected from this code:

1) Combined Volcano Plot with panel for each subtype differentially expressed genes. 
2) Combined Bargraph of pathway enrichment with panel for each subtype 
3) Combined GO Semantic Similarity Matrix with panel for each subtype 
