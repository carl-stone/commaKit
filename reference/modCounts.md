# Accessor for observed modified-read counts

Retrieves the sites × samples matrix of observed reads called as the
target modification. This assay is available for callers that report
count-like methylation evidence directly, such as modkit pileup. Older
objects or probability-only callers may contain `NA` values.

## Usage

``` r
modCounts(object)

# S4 method for class 'commaData'
modCounts(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

An integer matrix with rows corresponding to methylation sites and
columns corresponding to samples.

## See also

[`canonicalCounts`](https://carl-stone.github.io/commaKit/reference/canonicalCounts.md),
[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md)

## Examples

``` r
data(comma_example_data)
mod <- modCounts(comma_example_data)
dim(mod)
#> [1] 588   6
```
