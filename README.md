# ATACCoGAPS

Package which provides tools for processing and analysis of single-cell ATAC-seq data with the Bayesian Non-Negative Matrix Factorization algorithm, CoGAPS, as described in our Nucleic Acids Research paper https://doi.org/10.1093/nar/gkaa349.


Install ATACCoGAPS with:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ATACCoGAPS")
```

A tutorial demonstrating the standard analysis pipeline for ATACCoGAPS using data from Schep et al, 2017 can be found here: https://rossinerbe.github.io/ATACCoGAPS_Tutorial
