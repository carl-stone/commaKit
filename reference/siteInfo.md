# Accessor for per-site metadata

Returns the per-site metadata table from a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object. Reconstructs a flat `DataFrame` from the `rowRanges()` GRanges,
combining genomic coordinates (chrom, position, strand) with the mcols
columns (mod_type, motif, mod_context, plus any annotation/result
columns). This provides a backward-compatible interface to the
pre-Schema-v2 `rowData()` layout.

## Usage

``` r
siteInfo(object)

# S4 method for class 'commaData'
siteInfo(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with one row per methylation site. Always contains columns `chrom`,
`position`, `strand`, `mod_type`, `motif` (the sequence context),
`mod_context` (the composite modification context, e.g., `"6mA_GATC"`),
and `site_key` (a human-readable label with fixed
`"chrom:position:strand:mod_type:motif"` fields, e.g.,
`"chr1:512:+:6mA:GATC"`; computed on demand, not used for internal
matching). May contain additional annotation columns added by
[`annotateSites()`](https://carl-stone.github.io/commaKit/reference/annotateSites.md)
or result columns from
[`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md).

## See also

[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md),
[`modTypes`](https://carl-stone.github.io/commaKit/reference/modTypes.md)

## Examples

``` r
data(comma_example_data)
head(siteInfo(comma_example_data))
#> DataFrame with 6 rows and 8 columns
#>         chrom  position      strand mod_type       motif   is_diff mod_context
#>   <character> <integer> <character> <factor> <character> <logical> <character>
#> 1     chr_sim       443           +      6mA        GATC     FALSE    6mA_GATC
#> 2     chr_sim       512           +      6mA        GATC     FALSE    6mA_GATC
#> 3     chr_sim      1024           +      6mA        GATC     FALSE    6mA_GATC
#> 4     chr_sim      1073           +      6mA        GATC     FALSE    6mA_GATC
#> 5     chr_sim      1536           -      6mA        GATC     FALSE    6mA_GATC
#> 6     chr_sim      1602           +      6mA        GATC      TRUE    6mA_GATC
#>                 site_key
#>              <character>
#> 1 chr_sim:443:+:6mA:GATC
#> 2 chr_sim:512:+:6mA:GATC
#> 3 chr_sim:1024:+:6mA:G..
#> 4 chr_sim:1073:+:6mA:G..
#> 5 chr_sim:1536:-:6mA:G..
#> 6 chr_sim:1602:+:6mA:G..
```
