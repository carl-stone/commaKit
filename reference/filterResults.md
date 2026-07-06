# Filter differential methylation results by significance thresholds

A convenience wrapper around
[`results`](https://carl-stone.github.io/commaKit/reference/results.md)
that filters sites by adjusted p-value and absolute effect size.

## Usage

``` r
filterResults(object, ...)

# S4 method for class 'commaData'
filterResults(
  object,
  padj = 0.05,
  delta_beta = 0.1,
  mod_type = NULL,
  motif = NULL,
  mod_context = NULL,
  result = NULL,
  name = NULL,
  result_name = NULL,
  ...
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object on which
  [`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
  has been run.

- ...:

  Ignored.

- padj:

  Numeric. Maximum adjusted p-value threshold (inclusive). Default
  `0.05`.

- delta_beta:

  Numeric. Minimum absolute effect size threshold (\\\|\Delta\beta\|\\)
  (inclusive). Default `0.1`. Set to `0` to disable filtering on effect
  size.

- mod_type:

  Character vector or `NULL`. Passed to
  [`results`](https://carl-stone.github.io/commaKit/reference/results.md)
  for optional modification type filtering.

- motif:

  Character vector or `NULL`. Passed to
  [`results`](https://carl-stone.github.io/commaKit/reference/results.md)
  for optional sequence context motif filtering.

- mod_context:

  Character vector or `NULL`. Passed to
  [`results`](https://carl-stone.github.io/commaKit/reference/results.md)
  for optional modification context filtering (e.g., `"6mA_GATC"`).

- result:

  Character string or `NULL`. Name of a differential methylation result
  layer to filter. If `NULL` (default), the active default layer is
  used.

- name:

  Character string or `NULL`. Alias for `result`; provided for concise
  layer selection.

- result_name:

  Character string or `NULL`. Alias for `result`; provided for
  consistency with
  [`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md).

## Value

A `data.frame` (same format as
[`results`](https://carl-stone.github.io/commaKit/reference/results.md))
containing only sites where `dm_padj <= padj` **and**
`abs(dm_delta_beta) >= delta_beta`. Sites with `NA` values in either
column are excluded.

## See also

[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md),
[`results`](https://carl-stone.github.io/commaKit/reference/results.md)

## Examples

``` r
data(comma_example_data)
dm <- diffMethyl(comma_example_data, formula = ~condition, mod_type = "6mA")
#> diffMethyl: testing 'condition' -- 'treatment' vs 'control' (reference)
#> methylKit: comparing 'treatment' (treatment) vs 'control' (reference/control)
#> uniting...
sig <- filterResults(dm, padj = 0.05, delta_beta = 0.2)
nrow(sig)
#> [1] 31
```
