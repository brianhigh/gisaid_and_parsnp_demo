# GISAID and parsnp demo

This is an example of the following workflow:

- Download genomes from [GISAID](https://gisaid.org) as FASTA files using [GISAIDR](https://github.com/Wytamma/GISAIDR)
- Download a [reference genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) as a FASTA file using [ape](http://ape-package.ird.fr/)
- Align the genomes with [libMUSCLE](https://bioconda.github.io/recipes/libmuscle/README.html) using [parsnp](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
- Produce a distance matrix using [snp-dists](https://github.com/tseemann/snp-dists)
- Produce a phylogenetic tree using [parsnp](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
- Plot the phylogenetic tree with [ape](http://ape-package.ird.fr/) and save as a PDF file

To run this workflow, run the R script [gisaid_and_parsnp_demo.R](code/gisaid_and_parsnp_demo.R) after addressing the two prerequisites mentioned in the beginning comments of that script.

Tested on Ubuntu 22.04.3 LTS (Intel Core i7 CPU) and macOS Monterey 12.7.3 (Intel Core i5 CPU). The conda modules do not seem to work on Apple Silicon Macs or on Windows systems (but might work with WSL).
