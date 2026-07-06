# Run differential methylation via methylKit

An internal wrapper that uses methylKit to test for differential
methylation. Called by
[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
when `method = "methylkit"`.

## Usage

``` r
.runMethylKit(
  methyl_mat,
  coverage_mat,
  site_df,
  coldata,
  formula,
  ref_level = NULL,
  design_info = NULL,
  mod_counts_mat = NULL,
  canonical_counts_mat = NULL,
  other_mod_counts_mat = NULL
)
```

## Arguments

- methyl_mat:

  Numeric matrix (sites × samples) of beta values.

- coverage_mat:

  Integer matrix (sites × samples) of read depths.

- site_df:

  Data frame with columns `chrom`, `position`, `strand`, `mod_type`,
  `motif` — one row per site.

- coldata:

  `data.frame` with at least one column matching the RHS variable in
  `formula`.

- formula:

  One-sided formula (e.g., `~ condition`).

- ref_level:

  Optional reference level for the two-level contrast.

- design_info:

  Optional precomputed design information from
  [`.resolveDiffMethylDesign()`](https://carl-stone.github.io/commaKit/reference/dot-resolveDiffMethylDesign.md).

- mod_counts_mat:

  Optional integer matrix of observed modified-read counts. If supplied,
  these counts are preferred over reconstructing from beta values.

- canonical_counts_mat:

  Optional integer matrix of observed canonical-read counts.

- other_mod_counts_mat:

  Optional integer matrix of observed non-target modified-read counts.
  When present with `canonical_counts_mat`, these counts are included in
  the non-target denominator so it matches coverage.

## Value

A `data.frame` with the same columns as other differential methylation
backends: `pvalue`, `delta_beta`, and one `mean_beta_<level>` column per
condition level. The methylKit-specific `qvalue` column is preserved for
callers that want backend-specific evidence in addition to commaKit's
adjusted p-values. Row names are site keys.

## Details

methylKit must be installed (it is listed in `Suggests`). If it is not
available, this function stops with an informative message.

The wrapper converts the methylation and coverage matrices from a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object into the format expected by
[`methylKit::methylRawList`](https://rdrr.io/pkg/methylKit/man/methylRawList-class.html),
runs
[`methylKit::unite()`](https://rdrr.io/pkg/methylKit/man/unite-methods.html)
and
[`methylKit::calculateDiffMeth()`](https://rdrr.io/pkg/methylKit/man/calculateDiffMeth-methods.html),
and returns results in the same standardised format as
`.betaBinomialTest()`.

Only the first RHS variable of `formula` is used as the grouping
variable. Complex formulas with interactions or batch terms are not
currently supported by this wrapper.
