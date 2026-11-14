# Example image for preprocessing

An example image used to demonstrate preprocessing step for identifying
the region of study in each image.

## Usage

``` r
example_image
```

## Format

A data frame with columns:

- x:

  x spatial coordinate of the microbial cell

- y:

  y spatial coordinate of the microbial cell

- taxon_name:

  taxon identity of the microbial cell

## Source

Derived from data described in Grodner et al. (2024),
<https://doi.org/10.1101/2024.07.24.604873>.

## Examples

``` r
head(example_image)
#>           x         y   taxon_name
#> 1 12.194920  18.12571    Treponema
#> 2  1.084668 155.86324 Leptotrichia
#> 3  3.083120 218.79823 Leptotrichia
#> 4  1.029882 248.80400 Leptotrichia
#> 5  5.074240 292.04723 Leptotrichia
#> 6  2.562239 303.77835 Leptotrichia
```
