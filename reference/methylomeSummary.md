# Summarize per-sample methylation and coverage distributions

Computes per-sample summary statistics for methylation beta values and
sequencing coverage in a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object. Returns a tidy `data.frame` suitable for direct use with ggplot2
or for tabular reporting.

## Usage

``` r
methylomeSummary(object, mod_type = NULL, motif = NULL, mod_context = NULL)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object.

- mod_type:

  Character vector or `NULL`. If provided, only sites of the specified
  modification type(s) (e.g., `"6mA"`, `c("6mA", "5mC")`) are included
  in the summary. If `NULL` (default), all modification types are
  summarized together.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are included (e.g., `"GATC"`). If `NULL`
  (default), all motifs are included.

- mod_context:

  Character vector or `NULL`. If provided, only sites with a matching
  modification context are included (e.g., `"6mA:GATC"`). Applied after
  any `mod_type` and `motif` filters. Use
  [`modContexts`](https://carl-stone.github.io/commaKit/reference/modContexts.md)
  to see available values.

## Value

A `data.frame` with one row per sample, containing:

- `sample_name`:

  Sample identifier.

- `condition`:

  Experimental condition, from the optional `condition` column in
  `sampleInfo(object)`, or `NA` when that metadata is absent.

- `mod_type`:

  The modification type summarized (`"all"` if `mod_type = NULL`).

- `n_sites`:

  Total number of sites considered after filters.

- `n_covered`:

  Number of sites with non-`NA` methylation in this sample (i.e., sites
  above the coverage threshold); this is the denominator for beta
  summaries.

- `mean_beta`:

  Mean beta value across covered sites.

- `median_beta`:

  Median beta value across covered sites.

- `sd_beta`:

  Standard deviation of beta values across covered sites.

- `frac_methylated`:

  Fraction of covered sites with \\\beta \> 0.5\\ (broadly methylated).

- `mean_coverage`:

  Mean sequencing depth across non-missing coverage values for retained
  sites, including sites below the `min_coverage` threshold when
  coverage is available.

- `median_coverage`:

  Median sequencing depth across non-missing coverage values for
  retained sites, including sites below the `min_coverage` threshold
  when coverage is available.

- `min_coverage`:

  Minimum coverage threshold applied at construction, or `NA` if not
  stored.

## Details

Methylation beta summaries (`mean_beta`, `median_beta`, `sd_beta`, and
`frac_methylated`) are computed over covered sites: sites with non-`NA`
beta values after the object's coverage threshold has been applied.
Coverage summaries (`mean_coverage` and `median_coverage`) are computed
over non-missing coverage values for sites retained after filtering,
including sites whose beta value is `NA` in a sample because the site
did not meet the coverage threshold.

## See also

[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md),
[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`sampleInfo`](https://carl-stone.github.io/commaKit/reference/sampleInfo.md)

## Examples

``` r
data(comma_example_data)
ms <- methylomeSummary(comma_example_data)
ms
#>   sample_name condition mod_type n_sites n_covered mean_beta median_beta
#> 1      ctrl_1   control      all     588       588 0.8654839   0.8929141
#> 2      ctrl_2   control      all     588       588 0.8705692   0.8959600
#> 3      ctrl_3   control      all     588       588 0.8638033   0.8918851
#> 4     treat_1 treatment      all     588       588 0.8357998   0.8864176
#> 5     treat_2 treatment      all     588       588 0.8369054   0.8893089
#> 6     treat_3 treatment      all     588       588 0.8388398   0.8866568
#>      sd_beta frac_methylated mean_coverage median_coverage min_coverage
#> 1 0.10704958       0.9897959      79.25340            79.0            5
#> 2 0.09647906       0.9948980      79.83333            81.0            5
#> 3 0.10609612       0.9897959      82.67347            83.0            5
#> 4 0.17009987       0.9421769      76.49490            76.5            5
#> 5 0.16870207       0.9404762      78.44218            78.0            5
#> 6 0.16809920       0.9455782      79.17517            76.5            5
# Beta summaries use n_covered; coverage summaries use non-missing coverage.
ms[, c("sample_name", "n_sites", "n_covered", "mean_beta", "mean_coverage")]
#>   sample_name n_sites n_covered mean_beta mean_coverage
#> 1      ctrl_1     588       588 0.8654839      79.25340
#> 2      ctrl_2     588       588 0.8705692      79.83333
#> 3      ctrl_3     588       588 0.8638033      82.67347
#> 4     treat_1     588       588 0.8357998      76.49490
#> 5     treat_2     588       588 0.8369054      78.44218
#> 6     treat_3     588       588 0.8388398      79.17517

# Summarize only 6mA sites
ms_6mA <- methylomeSummary(comma_example_data, mod_type = "6mA")
ms_6mA[, c("sample_name", "condition", "mean_beta", "n_covered")]
#>   sample_name condition mean_beta n_covered
#> 1      ctrl_1   control 0.9036813       393
#> 2      ctrl_2   control 0.9029957       393
#> 3      ctrl_3   control 0.9014299       393
#> 4     treat_1 treatment 0.8557059       393
#> 5     treat_2 treatment 0.8491166       393
#> 6     treat_3 treatment 0.8527682       393
```
