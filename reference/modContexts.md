# Return the modification contexts present in a commaData object

Returns the unique modification contexts stored in a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object. A `mod_context` is a composite string combining modification
type and sequence motif: `paste(mod_type, motif, sep = "_")` when motif
information is available (e.g., `"6mA_GATC"`, `"5mC_CCWGG"`), or just
`mod_type` when motif context is unavailable.

## Usage

``` r
modContexts(object)

# S4 method for class 'commaData'
modContexts(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A sorted character vector of unique `mod_context` strings present in
`rowData(object)$mod_context` (e.g., `c("5mC_CCWGG", "6mA_GATC")`).

## Details

All differential methylation analyses run independently per
`mod_context` group by default, preventing spurious pooling of
biologically distinct methylation events (e.g., 6mA at GATC motifs from
Dam methyltransferase versus any cytosine methylation detected at GATC
positions, which is likely artefactual).

## See also

[`modTypes`](https://carl-stone.github.io/commaKit/reference/modTypes.md),
[`motifs`](https://carl-stone.github.io/commaKit/reference/motifs.md),
[`subset`](https://rdrr.io/r/base/subset.html)

## Examples

``` r
data(comma_example_data)
modContexts(comma_example_data)
#> [1] "5mC_CCWGG" "6mA_GATC" 
```
