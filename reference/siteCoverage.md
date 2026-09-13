# Accessor for the sequencing coverage (read depth) matrix

Retrieves the sites × samples matrix of read depth from a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object. This package-specific accessor avoids overloading
[`IRanges::coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html),
whose conventional Bioconductor meaning is genomic/Rle coverage
computation.

## Usage

``` r
siteCoverage(object)

# S4 method for class 'commaData'
siteCoverage(object)

# S4 method for class 'ANY'
siteCoverage(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

An integer matrix with rows corresponding to methylation sites and
columns corresponding to samples.

## See also

[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md),
[`siteInfo`](https://carl-stone.github.io/commaKit/reference/siteInfo.md)

## Examples

``` r
data(comma_example_data)
cov <- siteCoverage(comma_example_data)
summary(as.vector(cov))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   10.00   44.00   79.00   79.31  114.00  150.00 
```
