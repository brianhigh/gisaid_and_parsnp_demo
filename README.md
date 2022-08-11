# GISAID and parsnp demo

This is an example of the following workflow:

- Download whole genomes from [GISAID](https://gisaid.org)
- Download a [reference file](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)
- Align the genomes with libMUSCLE (using [parsnp](https://harvest.readthedocs.io/en/latest/content/parsnp.html))
- Produce a distance matrix and phylogenetic tree (using parsnp)
- Plot the phylogenetic tree and save as a PDF file

To run this workflow, run the R script [scripts/gisaid_and_parsnip_demo.R](gisaid_and_parsnip_demo.R) after addressing the two prerequisites mentioned in the beginning comments of that script.
