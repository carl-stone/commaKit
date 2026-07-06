# Deprecated subset method for commaData objects

`subset(commaData)` is deprecated to avoid package-level masking of
[`base::subset()`](https://rdrr.io/r/base/subset.html). Use
[`filterSites`](https://carl-stone.github.io/commaKit/reference/filterSites.md)
for common commaData filters or `[` for index-based subsetting.

## Usage

``` r
subset.commaData(
  x,
  mod_type = NULL,
  condition = NULL,
  chrom = NULL,
  motif = NULL,
  mod_context = NULL,
  ...
)
```

## Arguments

- x:

  A `commaData` object.

- mod_type:

  Character vector or `NULL`. If provided, only sites with a matching
  modification type are kept (e.g., `"6mA"`).

- condition:

  Character vector or `NULL`. If provided, only samples matching the
  specified condition(s) are kept.

- chrom:

  Character vector or `NULL`. If provided, only sites on the specified
  chromosome(s) are kept.

- motif:

  Character vector or `NULL`. If provided, only sites with a matching
  sequence context motif are kept (e.g., `"GATC"`). Sites with `NA`
  motif values are excluded when this filter is active. Use
  [`motifs`](https://carl-stone.github.io/commaKit/reference/motifs.md)
  to see which motifs are present.

- mod_context:

  Character vector or `NULL`. If provided, only sites with a matching
  modification context are kept (e.g., `"6mA_GATC"`, `"5mC_CCWGG"`). A
  `mod_context` value is `paste(mod_type, motif, sep = "_")` when motif
  is available, or just `mod_type` for Dorado/Megalodon data. Use
  [`modContexts`](https://carl-stone.github.io/commaKit/reference/modContexts.md)
  to see which contexts are present. When provided, this filter is
  applied in addition to (ANDed with) any `mod_type` or `motif` filters.

- ...:

  Ignored.

## Value

A `commaData` object containing only the selected sites and samples.
