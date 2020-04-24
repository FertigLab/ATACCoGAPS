# ATACCoGAPS

Package which provides tools for processing and analysis of single-cell ATAC-seq data with the Bayesian Non-Negative Matrix Factorization algorithm, CoGAPS.

Install the most recent version of CoGAPS in R with:

```
devtools::install_github("FertigLab/CoGAPS")
```

Install other needed Bioconductor packages for ATACCoGAPS using:

```
BiocManager::install(c("GenomicRanges", "projectR", "TFBSTools", "GeneOverlap", "msigdbr", "motifmatchr", "chromVAR", "GenomicFeatures", "IRanges", "fgsea", "rGREAT", "Homo.sapiens", "Mus.musculus"))
```

Install ATACCoGAPS with:

```
devtools::install_github("FertigLab/ATACCoGAPS")
```

A tutorial demonstrating the standard analysis pipeline for ATACCoGAPS using data from Schep et al, 2017 can be found here: https://rossinerbe.github.io/ATACCoGAPS_Tutorial
