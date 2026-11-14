# Test for differential spatial patterns using MiSA

Performs hypothesis testing to detect differential spatial patterns
between two groups, for each taxon in `MiSA_data`, using one of four
available MiSA tests.

## Usage

``` r
MiSA(misa_data_list, method = "area", omnibus_r = NULL, omnibus_num = NULL)
```

## Arguments

- misa_data_list:

  A named list, with one element per taxon, as returned by
  [`MiSA_data`](https://liuruoqian.github.io/MiSA/reference/MiSA_data.md)
  and
  [`MiSA_check`](https://liuruoqian.github.io/MiSA/reference/MiSA_check.md).
  Each element must contain `group`, a vector of group labels, and `d`,
  a matrix of estimated \\d\\-curves.

- method:

  Character string specifying which test to perform. One of `"orig"`,
  `"abs"`, `"sup"`, or `"area"`. Defaults to `"area"`.

- omnibus_r:

  Optional numeric vector specifying the specific radii at which to
  evaluate the test statistic, used only when `method = "orig"` or
  `method = "abs"`. Defaults to `NULL`.

- omnibus_num:

  Optional integer specifying the number of radii, evenly spaced across
  the available range, at which to evaluate the test statistic, used
  only when `method = "orig"` or `method = "abs"` and `omnibus_r` is not
  specified. Defaults to `NULL`.

## Value

A list containing taxon-specific p-values and a community-level p-value
assessing overall spatial dysbiosis. For `method = "area"` or
`method = "sup"`, taxon-specific effect sizes are additionally returned.
For `method = "orig"` or `method = "abs"`, a `per_r_pval` element is
additionally returned, describing the testing results at each of the
specified radii.

## Details

Four tests are available via the `method` argument:

- `"orig"`:

  Evaluates the raw statistic \\d(r)\\ (capturing both the magnitude and
  direction of deviation from CSR) at a set of pre-selected radii, and
  combines the resulting p-values across radii using the harmonic mean
  p-value (HMP) method.

- `"abs"`:

  Similar to `"orig"`, but compares the absolute value of \\d(r)\\
  (capturing only the magnitude of deviation) at the selected radii
  instead.

- `"sup"`:

  Evaluates the maximum absolute deviation, \\\max_r \|d(r)\|\\,
  capturing signals that may arise at any spatial scale.

- `"area"`:

  Evaluates the integrated magnitude of deviation of \\d(r)\\ across the
  full range of radii, capturing signals that may arise at any spatial
  scale.

For `"orig"` and `"abs"`, the radii used can be specified via
`omnibus_r` or `omnibus_num`. If `omnibus_r` is specified, the function
always uses `omnibus_r`, regardless of the input for `omnibus_num`. If
`omnibus_r` is left at its default (`NULL`) and `omnibus_num` is
specified, the function uses `omnibus_num` evenly spaced radii instead.
If both `omnibus_r` and `omnibus_num` are left at their defaults
(`NULL`), the function sets `omnibus_num = min(5, c - 1)`, where `c` is
the number of radii available for that taxon.

This function uses whichever group is entered as `groupA` (during data
preparation) as the reference group. A negative effect size indicates
that the summary statistic is smaller in the reference group than in the
comparison group, indicating that the comparison group is more spatially
structured (e.g., exhibiting stronger clustering or regularity) relative
to the reference group.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- MiSA(misa_data_list = data_checked, method = "area")
result <- MiSA(misa_data_list = data_checked, method = "orig", omnibus_r = c(300, 600, 900))
result <- MiSA(misa_data_list = data_checked, method = "abs", omnibus_num = 6)
} # }
```
