# Volcano plot for differential methylation results

Produces a volcano plot from differential methylation results computed
by
[`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
and extracted with
[`results()`](https://carl-stone.github.io/commaKit/reference/results.md).
Sites are colored by significance category based on user-supplied
adjusted p-value and delta-beta thresholds.

## Usage

``` r
plot_volcano(
  results,
  delta_beta_threshold = NULL,
  padj_threshold = 0.05,
  facet = TRUE
)
```

## Arguments

- results:

  A `data.frame` returned by
  [`results()`](https://carl-stone.github.io/commaKit/reference/results.md),
  containing at minimum the columns `dm_delta_beta` (numeric, effect
  size as beta difference treatment minus control) and `dm_padj`
  (numeric, BH-adjusted p-value in \\\[0, 1\]\\).

- delta_beta_threshold:

  `NULL` (default) or a numeric scalar in (0, 1). When `NULL`,
  significance is determined by `padj_threshold` alone and no vertical
  lines are drawn. When numeric, sites must also satisfy
  `|dm_delta_beta| >= delta_beta_threshold` to be called significant,
  and dashed vertical lines are drawn at \\\pm\\`delta_beta_threshold`.

- padj_threshold:

  Numeric scalar in (0, 1). Adjusted p-value cutoff for significance.
  Default `0.05`.

- facet:

  Logical. Default `TRUE`. When `TRUE` and `results` contains a
  `mod_context` column with more than one level, the plot is faceted by
  `mod_context`. Set to `FALSE` to suppress faceting.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The x-axis shows `dm_delta_beta` (effect size), the y-axis shows
`-log10(dm_padj)` (significance). Sites are colored as
`"Hypermethylated"` (positive delta-beta, significant),
`"Hypomethylated"` (negative delta-beta, significant), or
`"Not significant"`. A dashed horizontal line marks `padj_threshold`;
dashed vertical lines at \\\pm\\`delta_beta_threshold` are added only
when that argument is non-`NULL`. When `facet = TRUE` and multiple
`mod_context` levels are present, panels are split by `mod_context` via
`facet_wrap`. Sites with `NA` adjusted p-value are excluded from the
plot.

## Details

For the y-axis calculation, adjusted p-values are clamped from below at
`.Machine$double.xmin` before applying `-log10()`. This keeps exact-zero
adjusted p-values finite in the plotted coordinates instead of producing
infinite values.

## See also

[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md),
[`results`](https://carl-stone.github.io/commaKit/reference/results.md),
[`filterResults`](https://carl-stone.github.io/commaKit/reference/filterResults.md)

## Examples

``` r
data(comma_example_data)
cd_dm <- diffMethyl(comma_example_data, ~condition)
#> diffMethyl: testing 'condition' -- 'treatment' vs 'control' (reference)
#> methylKit: comparing 'treatment' (treatment) vs 'control' (reference/control)
#> uniting...
#> methylKit: comparing 'treatment' (treatment) vs 'control' (reference/control)
#> uniting...
res <- results(cd_dm)
plot_volcano(res)


# Custom thresholds
plot_volcano(res, delta_beta_threshold = 0.3, padj_threshold = 0.01)

```
