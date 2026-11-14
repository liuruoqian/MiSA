# Plot Observation Window Regions

Draws the borders of all regions returned by
[`alpha_windows`](https://liuruoqian.github.io/MiSA/reference/alpha_windows.md)
on a single plot, with each region labeled by its name.

## Usage

``` r
plot_regions(wins, main = "Alpha-shape regions")
```

## Arguments

- wins:

  A named list of `owin` objects, as returned by
  [`alpha_windows`](https://liuruoqian.github.io/MiSA/reference/alpha_windows.md).

- main:

  Plot title. Default is `"Alpha-shape regions"`.

## Value

Invisibly returns the union of all regions as a single `owin`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_image)
wins <- alpha_windows(kde_est(example_image))
plot_regions(wins)} # }
```
