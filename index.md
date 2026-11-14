# MiSA

MiSA identifies disease-associated spatial dysbiosis by testing whether
the spatial organization of microbial taxa differs between clinical
groups, using image-based microbiome data.

## Installation

You can install MiSA in R by running the code:

``` r

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("liuruoqian/MiSA", build_vignettes = TRUE)
```

Load MiSA in R and view the vignette for details for using MiSA.

``` r

library(MiSA)
browseVignettes("MiSA")
```

## Tutorial

See the [MiSA
vignette](https://liuruoqian.github.io/MiSA/articles/MiSA.html) for a
step-by-step introduction.
