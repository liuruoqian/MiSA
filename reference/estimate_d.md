# Estimate d-curves from L-function estimates for two groups

Computes the spatial summary function \\d(r)\\ for each taxon and image,
based on previously estimated \\L\\-function curves for two groups
(e.g., health vs. disease). The function \\d(r)\\ quantifies the
deviation of each point pattern from complete spatial randomness (CSR).

## Usage

``` r
estimate_d(groupA_L, groupB_L, taxa, r_vec)
```

## Arguments

- groupA_L:

  A list of estimated \\L\\-function curves for group A, as returned by
  [`estimate_L`](https://liuruoqian.github.io/MiSA/reference/estimate_L.md).

- groupB_L:

  A list of estimated \\L\\-function curves for group B, as returned by
  [`estimate_L`](https://liuruoqian.github.io/MiSA/reference/estimate_L.md).

- taxa:

  A character vector giving the taxa of interest for which \\d(r)\\ will
  be computed.

- r_vec:

  A numeric vector giving the sequence of radii at which \\d(r)\\ will
  be evaluated. This should match the sequence of radii used when
  estimating `groupA_L` and `groupB_L`.

## Value

A list containing the estimated \\d\\-curve and the corresponding group
label (`groupA` or `groupB`) for each image, for each taxon.

## Examples

``` r
if (FALSE) { # \dontrun{
r_all <- seq(0, 1500, by = 6)
taxa_all <- unique(healthImg.list[[1]]$taxon_name)
L_result <- estimate_L(healthImg.list, diseaseImg.list,
                        taxa = taxa_all, r_vec = r_all)
d_result <- estimate_d(groupA_L = est_L$groupA_L, groupB_L = est_L$groupB_L,
                        taxa = taxa_all, r_vec = r_all)
} # }
```
