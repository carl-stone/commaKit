# List differential methylation result layers

Returns a registry of named differential methylation analyses stored in
a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object. Each call to
[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
can write a named result layer, allowing multiple analysis runs to
coexist while the default layer remains available through
[`results`](https://carl-stone.github.io/commaKit/reference/results.md).

## Usage

``` r
resultLayers(object)

# S4 method for class 'commaData'
resultLayers(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with one row per result layer and columns describing the layer name,
default status, statistical method, result columns, filters, and
provenance.

## See also

[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md),
[`results`](https://carl-stone.github.io/commaKit/reference/results.md),
[`filterResults`](https://carl-stone.github.io/commaKit/reference/filterResults.md)

## Examples

``` r
data(comma_example_data)
dm <- diffMethyl(comma_example_data, formula = ~ condition,
                 method = "quasi_f", result_name = "quasi_f.v1")
#> diffMethyl: testing 'condition' -- 'treatment' vs 'control' (reference)
resultLayers(dm)
#> DataFrame with 1 row and 18 columns
#>          name        role                   type      source is_default
#>   <character> <character>            <character> <character>  <logical>
#> 1  quasi_f.v1  diffMethyl differential_methyla..  diffMethyl       TRUE
#>        method     formula   reference   treatment     mod_context
#>   <character> <character> <character> <character> <CharacterList>
#> 1     quasi_f  ~condition     control   treatment                
#>          mod_type           motif p_adjust_method min_coverage     alpha
#>   <CharacterList> <CharacterList>     <character>    <integer> <numeric>
#> 1                                              BH            5       0.5
#>                           result_cols              timestamp package_version
#>                       <CharacterList>            <character>     <character>
#> 1 dm_pvalue,dm_padj,dm_delta_beta,... 2026-06-23 20:26:12 ..           0.2.0
```
