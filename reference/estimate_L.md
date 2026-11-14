# Estimate L curves for two groups of images

Estimates Ripley's \\L\\-function for each taxon of interest, across all
images in two groups (e.g., health vs. disease), over a pre-specified
sequence of radii.

## Usage

``` r
estimate_L(
  groupA_image_list,
  groupB_image_list,
  taxa,
  r_vec,
  groupA_window_list = NULL,
  groupB_window_list = NULL,
  border_correct = "border"
)
```

## Arguments

- groupA_image_list:

  A named list of data frames representing the images in group A. Each
  data frame must contain the columns `x`, `y`, and `taxon_name`,
  representing the x and y spatial coordinates and the corresponding
  taxon identity, respectively.

- groupB_image_list:

  A named list of data frames representing the images in group B, with
  the same structure as `groupA_image_list`.

- taxa:

  A character vector giving the taxa of interest for which the
  \\L\\-function will be estimated.

- r_vec:

  A numeric vector giving the sequence of radii at which the
  \\L\\-function will be evaluated. This set of radii is used
  consistently across all images.

- groupA_window_list:

  (optional) A named list specifying pre-defined study regions
  (observation windows) for the images in `groupA_image_list`. Each
  element must be an object of class `owin`, as required by the `window`
  argument of
  [`spatstat.geom::ppp()`](https://rdrr.io/pkg/spatstat.geom/man/ppp.html).
  Names must match the corresponding image names in `groupA_image_list`.
  A window may be specified for some images and left as `NULL` for
  others. If `NULL` is provided (either for a single image or the entire
  list, the default), the function derives a convex hull from the data
  as the default window.

- groupB_window_list:

  (optional) A named list specifying pre-defined study regions for the
  images in `groupB_image_list`, with the same structure and default
  behavior as `groupA_window_list`.

- border_correct:

  (optional) Character string specifying the edge correction method used
  when estimating the \\L\\-function. One of `"border"` (the default),
  `"iso"`, or `"trans"`. For spatial microbiome data, `"border"` is
  recommended, as the other two options substantially increase
  computation time.

## Value

A list containing the estimated \\L\\ curves for each image and taxon.

## Examples

``` r
if (FALSE) { # \dontrun{
r_all <- seq(0, 1500, by = 6)
taxa_all <- unique(healthImg.list[[1]]$taxon_name)

L_result <- estimate_L(groupA_image_list=healthImg.list, groupB_image_list=diseaseImg.list,
taxa=taxa_all, r_vec=r_all,groupA_window_list = NULL, groupB_window_list = NULL,
border_correct = "border")
} # }
```
