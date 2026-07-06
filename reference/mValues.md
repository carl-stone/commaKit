# Compute M-values from a commaData object

Converts per-site beta values (methylation fractions) and read depths
into M-values using a pseudocount-offset logit transformation. M-values
are variance-stabilized relative to beta values and are better suited
for distance-based analyses such as PCA or hierarchical clustering.

## Usage

``` r
mValues(object, alpha = 0.5, mod_type = NULL, motif = NULL, mod_context = NULL)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object containing methylation beta values and coverage (read depth)
  assays.

- alpha:

  Positive numeric scalar. Pseudocount added to both the methylated and
  unmethylated read counts before log-transformation. Prevents infinite
  values at beta = 0 or beta = 1 and corresponds to a symmetric
  Beta(alpha, alpha) prior on the methylation fraction. Default `0.5`.

- mod_type:

  Character vector of modification types to include (e.g., `"6mA"`,
  `c("6mA", "5mC")`). If `NULL` (default), M-values are computed for all
  sites in `object`.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are included. If `NULL` (default), all
  motifs are included.

- mod_context:

  Character vector or `NULL`. If provided, only sites with a matching
  modification context are included (e.g., `"6mA_GATC"`). Applied after
  any `mod_type` and `motif` filters.

## Value

A numeric matrix of M-values with the same dimensions and `dimnames` as
`methylation(object)` (or the subset of rows matching `mod_type` if
specified). Positive values indicate hypermethylation; negative values
indicate hypomethylation; zero corresponds to a beta value of
approximately 0.5.

## Details

The M-value for a site in one sample is computed as:

\$\$M = \log_2\\\left(\frac{M\_{\mathrm{reads}} +
\alpha}{U\_{\mathrm{reads}} + \alpha}\right)\$\$

where \\M\_{\mathrm{reads}}\\ and \\U\_{\mathrm{reads}}\\ are observed
modified and canonical/unmodified read counts when count assays are
available. For legacy or probability-only objects, counts are
reconstructed from `round(beta * coverage)` as a documented fallback.
The pseudocount offset is \\\alpha\\.

Sites with zero coverage or `NA` beta values are returned as `NA`. The
pseudocount `alpha` must be strictly positive to avoid `-Inf` or `NaN`
values in the output.

## See also

[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md),
[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`modCounts`](https://carl-stone.github.io/commaKit/reference/modCounts.md),
[`canonicalCounts`](https://carl-stone.github.io/commaKit/reference/canonicalCounts.md),
[`plot_pca`](https://carl-stone.github.io/commaKit/reference/plot_pca.md)

## Examples

``` r
data(comma_example_data)

# Compute M-values for all modification types
m <- mValues(comma_example_data)
dim(m) # same as dim(methylation(comma_example_data))
#> [1] 588   6
range(m, na.rm = TRUE)
#> [1] -4.523562  7.851749

# Only 6mA sites
m6 <- mValues(comma_example_data, mod_type = "6mA")

# Use a smaller pseudocount
m_tight <- mValues(comma_example_data, alpha = 0.1)
```
