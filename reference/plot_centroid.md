# Plot cell centroids for a single image, colored by taxon

Creates a scatter plot of microbial cell centroids for a single image,
with points colored by taxon identity.

## Usage

``` r
plot_centroid(image.df, image.title = NULL)
```

## Arguments

- image.df:

  A data frame representing a single image, containing the columns `x`,
  `y`, and `taxon_name`, corresponding to the x and y spatial
  coordinates and the taxon identity of each cell, respectively.

- image.title:

  (optional) Character string giving the plot title. Defaults to `NULL`,
  in which case the title `"Scatter Plot by Taxon Name"` is used.

## Value

A `ggplot` object showing the spatial distribution of cell centroids,
colored by taxon.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_centroid(healthImg.list[['image1']])
plot_centroid(diseaseImg.list[['image1']], image.title = "Disease Image 1")
} # }
```
