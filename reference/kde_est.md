# Kernel density estimate of cell locations

Estimates the spatial density of segmented cells in an image, used to
identify the microbial colonization region before building the
observation window.

## Usage

``` r
kde_est(image)
```

## Arguments

- image:

  A data frame of segmented cells with columns `x`, `y` and
  `taxon_name`.

## Value

An object with class `"density"` whose underlying structure is a list.

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_image)
img_kde <- kde_est(example_image)} # }
```
