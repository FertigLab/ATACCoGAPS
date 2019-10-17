# ATACCoGAPS

Package which provides tools for processing and analysis of single-cell ATAC-seq data with the Bayesian Non-Negative Matrix Factorization algorithm, CoGAPS.

Install the most recent version of CoGAPS in R via:

```
devtools::install_github("FertigLab/CoGAPS")
```

or from Bioconductor with:

```
BiocManager::install("CoGAPS")
```

Install other needed Bioconductor packages for ATACCoGAPS using:

```
BiocManager::install(c("GenomicRanges", "projectR", "TFBSTools", "GeneOverlap", "msigdbr", "motifmatchr", "chromVAR", "GenomicFeatures", "IRanges", "Homo.sapiens", "Mus.musculus"))
```

Install ATACCoGAPS via:

```
devtools::install_github("FertigLab/ATACCoGAPS")
```

A tutorial demonstrating the standard analysis pipeline for ATACCoGAPS using data from Schep et al, 2017 can be found here: https://rossinerbe.github.io/
