# Return the modification types present in a commaData object

Returns the unique methylation modification types stored in a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object.

## Usage

``` r
modTypes(object)

# S4 method for class 'commaData'
modTypes(object)

# S4 method for class 'ANY'
modTypes(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A character vector of modification types present in
`rowData(object)$mod_type` (e.g., `c("6mA", "5mC")`).

## Examples

``` r
data(comma_example_data)
modTypes(comma_example_data)
#> [1] "5mC" "6mA"
```
