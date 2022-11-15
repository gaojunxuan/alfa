# alfa

Alignment-free sequence comparison methods for phylogenetic analysis.

## Description

alfa (Alignment-free analysis) will be an R package for carrying out comparative genomic analysis using alignment-free methods. The tool relies on k-mer based methods and different metrics to compute the similarity/dissimilarity given a set of sequences. The result can then be plotted as a heatmap or used to construct phylogenetic trees.

While there has been numerous packages for alignment-free and phylogenetic analysis on R, none of the packages aim to integrate the workflow or provide a set of tools specifically tailored for phylogenetic analysis using alignment-free methods. Our package aim to simplify the workflow for performing phylogenetic analysis and in particular creation and evaluation of phylogenetic trees using alignment-free approaches. The results obtained from running the analyses using this package can be used for downstream analyses or for boostrapping in a more conventional phylogeny building workflow.

The package is developed using R version 4.1.2 (2021-11-01) on a device running x86_64-apple-darwin17.0 (macOS 12.6).

## Installation

You can install the development version of alfa from [GitHub](https://github.com/) with:

``` r
require("devtools")
devtools::install_github("gaojunxuan/alfa", build_vignettes = TRUE)
library("alfa")
```

## Overview

The package provides the following functions for analysis.

-   `allKMers`: extract all k-mers of the given string
-   `countVectors`: create a count vector of k-mers in the two strings
-   `createDistanceMatrix`: create a pairwise distance matrix for the given set of strings
-   `euclideanDistnace`: calculate the Euclidean distance between two strings using their count vectors
-   `standardizedEuclidean`: calculate the standardized Euclidean distance between two strings using their count vectors, adjusted based on the overlapping capabilities of the k-mers
-   `neighborJoiningTree`: construct and plot phylogenetic tree using neighbor-joining
-   `upgmaTree`: construct and plot phylogenetic tree using UPGMA

As well as the following functions for plotting and visualizing the results

-   `neighborJoiningTree`: construct and plot phylogenetic tree using neighbor-joining
-   `upgmaTree`: construct and plot phylogenetic tree using UPGMA
-   `plotPairwiseDist`: plot pairwise distance matrix as a heatmap
-   `plotDistanceCorrelation`: plot the pairwise distance in the distance matrix against the pairwise distance on the given phylogenetic tree

## Contributions

This package is developed by Kevin Gao. The idea of using count vectors for distance calculation is common for alignment-free analysis but the initial inspiration was from Zielezinski et al. The algorithm used for computing the standardized Euclidean distance was based on the paper `A Measure of DNA Sequence Dissimilarity Based on Mahalanobis Distance between Frequencies of Words` by Wu, Burke, and Davison. Moreover, various functions in this package including the tree constructing functions and plotting functions relied on `ape` by Paradis and Schliep. The package also uses functions from the `stats` package for computation.

## References

Wu TJ, Burke JP, Davison DB. A measure of DNA sequence dissimilarity based on Mahalanobis distance between frequencies of words. Biometrics. 1997 Dec;53(4):1431-9.

Zielezinski, A., Vinga, S., Almeida, J. et al. Alignment-free sequence comparison: benefits, applications, and tools. Genome Biology. 2017; 18, 186.

Paradis E, Schliep K. "ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R." Bioinformatics. 2019; 35:526-528.

N Saitou, M Nei, The neighbor-joining method: a new method for reconstructing phylogenetic trees., Molecular Biology and Evolution. 1987 Jul; 4(4), 406--425.

## Acknowledgement

This package was developed as part of an assessment for 2022 BCB410H: Applied Bioinformat-\
ics course at the University of Toronto, Toronto, CANADA. `alfa` welcomes issues,\
enhancement requests, and other contributions. To submit an issue, use the GitHub issues.
