# Check taxa for validity prior to hypothesis testing

Checks each taxon against a set of minimum requirements for hypothesis
testing, and drops any taxon that does not meet these requirements.

## Usage

``` r
MiSA_check(data_list)
```

## Arguments

- data_list:

  A named list, with one element per taxon, as returned by
  [`MiSA_data`](https://liuruoqian.github.io/MiSA/reference/MiSA_data.md).

## Value

A named list with the same structure as `data`, containing only the taxa
that passed all validity checks.

## Details

A taxon is dropped if any of the following conditions are not met:

1.  At least two groups remain after filtering.

2.  The total sample size is greater than 2.

3.  The corresponding `d` has more than one valid estimate.

## Examples

``` r
if (FALSE) { # \dontrun{
data_checked <- MiSA_check(data_list = data_pre)
} # }
```
