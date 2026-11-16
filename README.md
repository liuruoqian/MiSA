
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Identifying disease-associated dysbiosis in image-based spatial microbiome via MiSA

<!-- badges: start -->

<!-- badges: end -->

# Overview

This package (MiSA) conducts analysis for identifying disease-associated
spatial dysbiosis by testing whether the spatial organization of
microbial taxa differs between clinical groups, using image-based
microbiome data.

## Installation

Before setting up the package, users should have R version 3.5 or
higher.

To install MiSA in R, run:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("liuruoqian/MiSA", build_vignettes = TRUE)
```

Load MiSA and see the vignette for a step-by-step guide.

``` r
library(MiSA)
browseVignettes("MiSA")
```

## Tutorial

A typical analysis follows these steps:

1.  **Preprocessing:** identify the region of study (i.e., the area
    colonized by microbes) in each image and build observation windows
    (`kde_est()`, `alpha_windows()`).
2.  **Taxon-level testing:** characterize the spatial distribution of
    each taxon with the summary function *d(r)* and compare it between
    groups using the four MiSA variants.
3.  **Community-level testing:** combine the taxon-level results to
    assess overall spatial dysbiosis.

Find the MiSA vignette at:
<https://liuruoqian.github.io/MiSA/articles/MiSA_vignette.html> for a
step-by-step introduction.

See MiSA website: <https://liuruoqian.github.io/MiSA/> for details.
