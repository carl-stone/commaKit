# Extract differential methylation results as a tidy data frame

Retrieves the per-site differential methylation statistics added by
[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
and returns them as a tidy `data.frame` suitable for downstream analysis
and plotting.

## Usage

``` r
results(object, ...)

# S4 method for class 'commaData'
results(
  object,
  mod_type = NULL,
  motif = NULL,
  mod_context = NULL,
  result = NULL,
  name = NULL,
  result_name = NULL,
  as = c("data.frame", "GRanges"),
  ...
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object on which
  [`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
  has been run.

- ...:

  Ignored (for S4 generic compatibility).

- mod_type:

  Character vector or `NULL`. If provided, only sites of the specified
  modification type(s) are returned. If `NULL` (default), results for
  all modification types are returned.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are returned. If `NULL` (default), all
  motifs are returned.

- mod_context:

  Character vector or `NULL`. If provided, only sites with a matching
  modification context (e.g., `"6mA_GATC"`) are returned. Use
  [`modContexts`](https://carl-stone.github.io/commaKit/reference/modContexts.md)`(object)`
  to see available values. Applied in addition to any `mod_type` or
  `motif` filters.

- result:

  Character string or `NULL`. Name of a differential methylation result
  layer to extract. If `NULL` (default), the active default layer is
  used.

- name:

  Character string or `NULL`. Alias for `result`; provided for concise
  layer selection.

- result_name:

  Character string or `NULL`. Alias for `result`; provided for
  consistency with
  [`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md).

- as:

  Character string. Output format: `"data.frame"` (default) or
  `"GRanges"`. `"GRanges"` returns filtered `rowRanges(object)` with the
  selected result columns in `mcols`.

## Value

A `data.frame` with one row per methylation site, containing:

- `chrom`:

  Chromosome name.

- `position`:

  1-based genomic position.

- `strand`:

  Strand (`"+"` or `"-"`).

- `mod_type`:

  Modification type (e.g., `"6mA"`).

- `dm_pvalue`:

  Raw p-value from the statistical test.

- `dm_padj`:

  Adjusted p-value (Benjamini-Hochberg by default).

- `dm_methylkit_qvalue`:

  methylKit q-value, present only when extracting a result layer
  produced by `method = "methylkit"`.

- `dm_delta_beta`:

  Effect size: mean methylation in the treatment group minus mean
  methylation in the reference group.

- `dm_mean_beta_<condition>`:

  One column per condition level with per-group mean beta values.

Any other annotation columns present in `rowData(object)` (e.g., from
[`annotateSites`](https://carl-stone.github.io/commaKit/reference/annotateSites.md))
are also included.

## See also

[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md),
[`filterResults`](https://carl-stone.github.io/commaKit/reference/filterResults.md)

## Examples

``` r
data(comma_example_data)
dm <- diffMethyl(comma_example_data, formula = ~condition, mod_type = "6mA")
#> diffMethyl: testing 'condition' -- 'treatment' vs 'control' (reference)
#> methylKit: comparing 'treatment' (treatment) vs 'control' (reference/control)
#> uniting...
res <- results(dm)
head(res[order(res$dm_padj), ])
#>       chrom position strand mod_type motif is_diff mod_context
#> 196 chr_sim    50176      -      6mA  GATC    TRUE    6mA_GATC
#> 287 chr_sim    70003      -      6mA  GATC    TRUE    6mA_GATC
#> 260 chr_sim    63550      +      6mA  GATC    TRUE    6mA_GATC
#> 249 chr_sim    61440      +      6mA  GATC    TRUE    6mA_GATC
#> 347 chr_sim    86016      +      6mA  GATC    TRUE    6mA_GATC
#> 9   chr_sim     2180      -      6mA  GATC    TRUE    6mA_GATC
#>                     site_key    dm_pvalue      dm_padj dm_methylkit_qvalue
#> 196 chr_sim:50176:-:6mA:GATC 4.705226e-78 1.849154e-75        1.259194e-75
#> 287 chr_sim:70003:-:6mA:GATC 1.982943e-70 3.896483e-68        2.653338e-68
#> 260 chr_sim:63550:+:6mA:GATC 3.822058e-68 5.006897e-66        3.409482e-66
#> 249 chr_sim:61440:+:6mA:GATC 1.199352e-66 1.178364e-64        8.024152e-65
#> 347 chr_sim:86016:+:6mA:GATC 4.658449e-64 3.661541e-62        2.493352e-62
#> 9    chr_sim:2180:-:6mA:GATC 6.144202e-62 4.024452e-60        2.740479e-60
#>     dm_delta_beta dm_mean_beta_control dm_mean_beta_treatment
#> 196    -0.7336497            0.8939269             0.16027720
#> 287    -0.7050844            0.8734775             0.16839309
#> 260    -0.7799241            0.8661683             0.08624415
#> 249    -0.7090099            0.9070218             0.19801189
#> 347    -0.6743832            0.8716203             0.19723706
#> 9      -0.7543758            0.9464297             0.19205391
```
