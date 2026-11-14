# Prepare d-curves for hypothesis testing

Prepares \\d\\-curves for downstream hypothesis testing by filtering out
images with limited spatial range and truncating the remaining curves to
a common maximum radius.

## Usage

``` r
MiSA_data(est_d, r_thresh = NULL)
```

## Arguments

- est_d:

  A list of estimated \\d\\-curves and group labels for each taxon and
  image, as returned by
  [`estimate_d`](https://liuruoqian.github.io/MiSA/reference/estimate_d.md).

- r_thresh:

  (optional) Numeric value giving the radius threshold used to filter
  out images with limited spatial range. Images whose \\d\\-curve does
  not have valid values up to `r_thresh` are excluded from further
  analysis. Defaults to `NULL`, in which case images whose \\d\\-curve
  does not reach 0.5 times the maximum radius considered while
  estimating \\d\\-curve. Users may set `r_thresh` based on the
  scientific question of interest; for example, `r_thresh = 500` will
  exclude \\d\\-curves that do not have valid estimates up to \\r =
  500\\.

## Value

A list of filtered and truncated \\d\\-curves, ready for input into
downstream hypothesis testing functions.

## Examples

``` r
if (FALSE) { # \dontrun{
data_pre <- MiSA_data(est_d)
data_pre <- MiSA_data(est_d, r_thresh = 500)
} # }
```
