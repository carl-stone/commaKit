# List assay layers and defaults

Returns a compact registry of assay matrices stored in a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object. The registry is derived from `assayNames(object)`,
`metadata(object)$assay_provenance`, and
`metadata(object)$assay_defaults`. This mirrors layered-assay behavior
used by Bioconductor and single-cell workflows: raw layers remain
present, while derived layers are named explicitly and can be marked as
defaults for a role without overwriting the raw data.

## Usage

``` r
assayLayers(object)

# S4 method for class 'commaData'
assayLayers(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with one row per assay and columns `assay`, `role`, `type`, `source`,
`is_default`, `default_for`, `parent_assays`, `method`, `timestamp`, and
`package_version`.

## Details

commaKit's v1 layer policy is intentionally minimal:

- Raw evidence assays such as `methylation`, `coverage`, `mod_counts`,
  `canonical_counts`, and `other_mod_counts` are canonical input
  evidence. Package APIs do not mutate them in place, though users may
  still edit assays manually through normal `SummarizedExperiment`
  mechanisms.

- Transformations should be stored as explicitly named derived assay
  layers with provenance and parent assays rather than replacing raw
  assays.

- Filtering returns a subset `commaData` object with all assay layers
  subset to the same rows/samples. There are no hidden lazy views.

## See also

[`assayProvenance`](https://carl-stone.github.io/commaKit/reference/assayProvenance.md),
[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md),
[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`modCounts`](https://carl-stone.github.io/commaKit/reference/modCounts.md),
[`canonicalCounts`](https://carl-stone.github.io/commaKit/reference/canonicalCounts.md)

## Examples

``` r
data(comma_example_data)
assayLayers(comma_example_data)
#> DataFrame with 4 rows and 10 columns
#>              assay             role                   type            source
#>        <character>      <character>            <character>       <character>
#> 1      methylation      methylation          filtered_beta synthetic_example
#> 2         coverage         coverage observed_total_cover.. synthetic_example
#> 3       mod_counts       mod_counts   reconstructed_counts synthetic_example
#> 4 canonical_counts canonical_counts   reconstructed_counts synthetic_example
#>   is_default      default_for        parent_assays                 method
#>    <logical>  <CharacterList>      <CharacterList>            <character>
#> 1       TRUE      methylation             coverage             simulation
#> 2       TRUE         coverage                                  simulation
#> 3       TRUE       mod_counts methylation,coverage round_beta_times_cov..
#> 4       TRUE canonical_counts  coverage,mod_counts coverage_minus_mod_c..
#>     timestamp package_version
#>   <character>     <character>
#> 1          NA           0.2.0
#> 2          NA           0.2.0
#> 3          NA           0.2.0
#> 4          NA           0.2.0
```
