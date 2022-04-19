# Hubbard-et-al-2022-Nature
Code used for analyses and visualizations for the paper "ADAR1 mutation causes ZBP1-dependent immunopathology"

## Environment Setup
We assume an instalation of R is available.
The environment is managed using `renv` if you do not have `renv` please open an R session and run

> install.packages('renv')

The R scripts will then install packages.

## RNA-seq Analysis

The Kallisto processed RNA seq data is available from [https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE200854&format=file](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE200854&format=file), it is assumed that the folder has been decompressed in the project `data` directory.

To produce the plots from the paper run `Rscript Rna_seq_analysis.R` which will produce the plots in a folder called `results`

## Nanostring Analysis

The processed data is provided in `data` and the paper plots are reproduced using `Rscript Nanostring_visualizations.R`
