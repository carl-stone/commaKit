# Heatmap of top differentially methylated sites

Produces a heatmap showing methylation beta values for the top
differentially methylated sites across all samples. Sites are ranked by
adjusted p-value and ordered vertically by effect size (delta beta) to
reveal condition-specific patterns.

## Usage

``` r
plot_heatmap(results, object, n_sites = 50L, annotation_cols = NULL)
```

## Arguments

- results:

  A `data.frame` returned by
  [`results()`](https://carl-stone.github.io/commaKit/reference/results.md),
  containing at minimum the columns `chrom`, `position`, `strand`,
  `mod_type`, `dm_padj`, and `dm_delta_beta`. Row names of the selected
  rows must be non-missing integer row indices into `object`, as
  produced by
  [`results`](https://carl-stone.github.io/commaKit/reference/results.md)
  and preserved by ordinary data-frame subsetting and sorting.

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object that was used to produce `results`. Used to extract the
  methylation matrix for selected sites.

- n_sites:

  Positive integer. The number of top sites (ranked by ascending
  `dm_padj`) to include in the heatmap. If fewer significant (non-`NA`
  padj) sites exist, all are shown. Default `50`.

- annotation_cols:

  Character vector naming columns from `sampleInfo(object)` to display
  as a colored annotation bar above the heatmap. If `NULL` (default),
  the `condition` column is used.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The x-axis shows sample names; the y-axis shows the selected
sites (y-axis labels suppressed for readability). The fill color encodes
methylation beta (blue = 0, white = 0.5, red = 1). `NA` values are shown
in light grey. An annotation strip above shows sample-level metadata
encoded by color.

## Details

`plot_heatmap()` maps differential methylation result rows back to
`object` through the row names on `results`: selected row names must be
integer row indices within `seq_len(nrow(object))`. Do not reset row
names or replace them with display labels such as site keys before
calling this function.

## See also

[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md),
[`results`](https://carl-stone.github.io/commaKit/reference/results.md),
[`plot_volcano`](https://carl-stone.github.io/commaKit/reference/plot_volcano.md)

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
plot_heatmap(res, cd_dm)


# Show only top 20 sites
plot_heatmap(res, cd_dm, n_sites = 20)

```
