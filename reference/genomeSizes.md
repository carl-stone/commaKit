# Accessor for chromosome size information

Returns the chromosome sizes stored in a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object. Size information is stored in the `Seqinfo` attached to
`rowRanges(object)` and corresponds to `seqlengths(object)`. This
package-specific accessor avoids overloading `GenomeInfoDb::genome()`,
whose conventional Bioconductor meaning is genome build/name metadata
rather than chromosome lengths.

## Usage

``` r
genomeSizes(object)

# S4 method for class 'commaData'
genomeSizes(object)

# S4 method for class 'ANY'
genomeSizes(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A named integer vector of chromosome sizes (chromosome name -\> length
in bp), or `NULL` if no size information was provided at construction.

## Examples

``` r
data(comma_example_data)
genomeSizes(comma_example_data)
#> chr_sim 
#>  100000 
```
