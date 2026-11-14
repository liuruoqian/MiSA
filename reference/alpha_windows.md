# Build observation windows from alpha shapes of high-density regions

Builds closed observation windows from alpha shapes of high-density cell
regions. Separated regions are kept apart, and each is returned as its
own `owin` object, sorted by area from largest to smallest.

## Usage

``` r
alpha_windows(image_kde, alpha = 200, cutoff = 0.002, fill_holes = TRUE)
```

## Arguments

- image_kde:

  Output of
  [`kde_est`](https://liuruoqian.github.io/MiSA/reference/kde_est.md).

- alpha:

  Positive number controlling how tightly the window follows the shape
  of the high-density cell region, in the same units as the coordinates.
  Smaller values give a tighter, more detailed boundary and may split
  the region into several pieces; larger values give a smoother boundary
  that bridges small gaps. As a rough guide, empty gaps wider than about
  `2 * alpha` separate regions or form holes. Default is 200.

- cutoff:

  Density threshold used to remove low-density areas. Default is 0.002.

- fill_holes:

  Logical. If `TRUE` (default), holes within each region are filled, so
  each window is bounded by its outer boundary only. If `FALSE`, holes
  are kept and excluded from the window.

## Value

A named list of `owin` objects (`region1`, `region2`, ...).

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_image)
wins <- alpha_windows(kde_est(example_image))
wins <- alpha_windows(kde_est(example_image), alpha = 100, cutoff = 2e-3, fill_holes = FALSE)
plot(wins$region1)} # }
```
