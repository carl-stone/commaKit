# Accessor for observed non-target modified-read counts

Retrieves the sites × samples matrix of observed reads called as a
non-target modification at the same site. For modkit pileup this is the
`Nother_mod` column, and `coverage` is the denominator
`mod_counts + canonical_counts + other_mod_counts`. Consult
[`assayProvenance`](https://carl-stone.github.io/commaKit/reference/assayProvenance.md)
for source details.

## Usage

``` r
otherModCounts(object)

# S4 method for class 'commaData'
otherModCounts(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

An integer matrix with rows corresponding to methylation sites and
columns corresponding to samples.

## See also

[`modCounts`](https://carl-stone.github.io/commaKit/reference/modCounts.md),
[`canonicalCounts`](https://carl-stone.github.io/commaKit/reference/canonicalCounts.md),
[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md)
