# Causes and consequences of clonal hematopoiesis in the UK Biobank

This is the repository, in conjunction with [https://github.com/siddhartha-kar/clonal-hematopoiesis], for the manuscript entitled [Genome-wide analyses of 200,453 individuals yields new insights into the causes and consequences of clonal hematopoiesis](https://www.medrxiv.org/content/10.1101/2022.01.06.22268846v1)

## Scripts

This repository is organized in the following directories:

- `variant calling` which includes the Nextflow pipelines used to run Mutect2 and to process the subsequent VCF files
  
- `figures` which includes a R notebook and R script to generate Fig 1, 2, and Extended Data Fig. 1-3.

NOTE: due to legal restrictions, we cannot provide data from the UK Biobank. To apply for access go to UK Biobank website [UK Biobank](https://www.ukbiobank.ac.uk/)

## Software

* `R` 4.0.2 - other versions may work, but this was the version used
* `Nextflow` version 20.07.1 - to run the pipelines. Other versions may work, but this was the version used
