# Accessor for observed canonical-read counts

Retrieves the sites × samples matrix of observed reads called as
canonical or unmodified for the site. Consult
[`assayProvenance`](https://carl-stone.github.io/commaKit/reference/assayProvenance.md)
for source details.

## Usage

``` r
canonicalCounts(object)

# S4 method for class 'commaData'
canonicalCounts(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

An integer matrix with rows corresponding to methylation sites and
columns corresponding to samples.

## See also

[`modCounts`](https://carl-stone.github.io/commaKit/reference/modCounts.md),
[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`otherModCounts`](https://carl-stone.github.io/commaKit/reference/otherModCounts.md),
[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md)

## Examples

``` r
data(comma_example_data)
canonical <- canonicalCounts(comma_example_data)
dim(canonical)
#> [1] 588   6
```
