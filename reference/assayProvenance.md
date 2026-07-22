# Accessor for assay provenance metadata

Returns metadata describing how assay matrices were produced. The
constructor records whether count assays were observed directly from
modkit output or reconstructed for synthetic/test objects.

## Usage

``` r
assayProvenance(object)

# S4 method for class 'commaData'
assayProvenance(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A named `list`. Returns an empty list for legacy objects without assay
provenance metadata.

## See also

[`assayLayers`](https://carl-stone.github.io/commaKit/reference/assayLayers.md),
[`modCounts`](https://carl-stone.github.io/commaKit/reference/modCounts.md),
[`canonicalCounts`](https://carl-stone.github.io/commaKit/reference/canonicalCounts.md),
[`otherModCounts`](https://carl-stone.github.io/commaKit/reference/otherModCounts.md),
[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md)

## Examples

``` r
data(comma_example_data)
assayProvenance(comma_example_data)
#> $methylation
#> $methylation$type
#> [1] "filtered_beta"
#> 
#> $methylation$source
#> [1] "synthetic_example"
#> 
#> $methylation$role
#> [1] "methylation"
#> 
#> $methylation$parent_assays
#> [1] "coverage"
#> 
#> $methylation$method
#> [1] "simulation"
#> 
#> $methylation$params
#> $methylation$params$min_coverage
#> [1] 5
#> 
#> 
#> $methylation$default_for
#> [1] "methylation"
#> 
#> $methylation$timestamp
#> [1] NA
#> 
#> $methylation$package_version
#> [1] "0.2.0"
#> 
#> 
#> $coverage
#> $coverage$type
#> [1] "observed_total_coverage"
#> 
#> $coverage$source
#> [1] "synthetic_example"
#> 
#> $coverage$role
#> [1] "coverage"
#> 
#> $coverage$parent_assays
#> character(0)
#> 
#> $coverage$method
#> [1] "simulation"
#> 
#> $coverage$params
#> list()
#> 
#> $coverage$default_for
#> [1] "coverage"
#> 
#> $coverage$timestamp
#> [1] NA
#> 
#> $coverage$package_version
#> [1] "0.2.0"
#> 
#> 
#> $mod_counts
#> $mod_counts$type
#> [1] "reconstructed_counts"
#> 
#> $mod_counts$source
#> [1] "synthetic_example"
#> 
#> $mod_counts$role
#> [1] "mod_counts"
#> 
#> $mod_counts$parent_assays
#> [1] "methylation" "coverage"   
#> 
#> $mod_counts$method
#> [1] "round_beta_times_coverage"
#> 
#> $mod_counts$params
#> list()
#> 
#> $mod_counts$default_for
#> [1] "mod_counts"
#> 
#> $mod_counts$timestamp
#> [1] NA
#> 
#> $mod_counts$package_version
#> [1] "0.2.0"
#> 
#> 
#> $canonical_counts
#> $canonical_counts$type
#> [1] "reconstructed_counts"
#> 
#> $canonical_counts$source
#> [1] "synthetic_example"
#> 
#> $canonical_counts$role
#> [1] "canonical_counts"
#> 
#> $canonical_counts$parent_assays
#> [1] "coverage"   "mod_counts"
#> 
#> $canonical_counts$method
#> [1] "coverage_minus_mod_counts"
#> 
#> $canonical_counts$params
#> list()
#> 
#> $canonical_counts$default_for
#> [1] "canonical_counts"
#> 
#> $canonical_counts$timestamp
#> [1] NA
#> 
#> $canonical_counts$package_version
#> [1] "0.2.0"
#> 
#> 
```
