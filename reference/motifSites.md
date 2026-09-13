# Accessor for motif site positions

Returns the
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html) of
all instances of the user-specified sequence motif in the genome, as
computed by
[`findMotifSites`](https://carl-stone.github.io/commaKit/reference/findMotifSites.md)
during object construction.

## Usage

``` r
motifSites(object)

# S4 method for class 'commaData'
motifSites(object)

# S4 method for class 'ANY'
motifSites(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object. May be empty (length 0) if no motif was specified at
construction.

## Examples

``` r
data(comma_example_data)
motifSites(comma_example_data)
#> GRanges object with 0 ranges and 0 metadata columns:
#>    seqnames    ranges strand
#>       <Rle> <IRanges>  <Rle>
#>   -------
#>   seqinfo: no sequences
```
