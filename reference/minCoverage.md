# Accessor for the minimum coverage threshold

Returns the minimum read depth threshold that was applied at
construction time. Sites with coverage below this threshold have their
beta value set to `NA`.

## Usage

``` r
minCoverage(object)

# S4 method for class 'commaData'
minCoverage(object)

# S4 method for class 'ANY'
minCoverage(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

An integer (the minimum coverage threshold), or `NA_integer_` if not
stored (e.g., objects created before min_coverage storage was
implemented).

## Examples

``` r
data(comma_example_data)
minCoverage(comma_example_data)
#> [1] 5
```
