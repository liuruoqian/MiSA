# Sample disease group imaging data

A named list of data frames containing example spatial microbiome
imaging data for the disease (e.g., mucositis) group, used to
demonstrate the MiSA workflow.

## Usage

``` r
diseaseImg.list
```

## Format

A named list, with one element per image. Each element is a data frame
with columns:

- x:

  x spatial coordinate of the microbial cell

- y:

  y spatial coordinate of the microbial cell

- taxon_name:

  taxon identity of the microbial cell

## Source

A simulated dataset generated using oral biofilm data as a template.
