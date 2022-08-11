# GISAID and parsnp demo

This is an example of the following workflow:

- Download genomes from [GISAID](https://gisaid.org) as FASTA files
- Download a [reference genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) as a FASTA file
- Align the genomes with libmuscle (using [parsnp](https://harvest.readthedocs.io/en/latest/content/parsnp.html))
- Produce a distance matrix and phylogenetic tree (using parsnp)
- Plot the phylogenetic tree and save as a PDF file

To run this workflow, run the R script [gisaid_and_parsnp_demo.R](code/gisaid_and_parsnp_demo.R) after addressing the two prerequisites mentioned in the beginning comments of that script.
