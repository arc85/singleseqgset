# singleseqgset

singleseqgset is a package for gene set enrichment analysis for single-cell RNAseq data. It uses a simple underlying statistic (variance inflated Wilcoxon rank sum testing) to determine enrichment of gene sets of interest across clusters. Check out the vignette for more information on the implementation.

## Installation

You can install the development version of singleseqgset from github using devtools:

``` r
library(devtools)
install_github("arc85/singleseqgset")
library(singleseqgset)
```
