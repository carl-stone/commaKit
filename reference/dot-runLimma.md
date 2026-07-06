# Per-site moderated t-test via limma eBayes for differential methylation

An internal wrapper that uses limma's empirical Bayes moderated t-test
to identify differentially methylated sites. Called by
[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
when `method = "limma"`.

## Usage

``` r
.runLimma(
  methyl_mat,
  coverage_mat,
  site_df,
  coldata,
  formula,
  alpha = 0.5,
  ref_level = NULL,
  design_info = NULL,
  mod_counts_mat = NULL,
  canonical_counts_mat = NULL,
  other_mod_counts_mat = NULL
)
```

## Arguments

- methyl_mat:

  Numeric matrix (sites × samples) of beta values. `NA` indicates
  below-coverage sites.

- coverage_mat:

  Integer matrix (sites × samples) of read depths.

- site_df:

  Data frame with columns `chrom`, `position`, `strand`, `mod_type`,
  `motif` — one row per site.

- coldata:

  `data.frame` with at least one column matching the RHS variable in
  `formula` (typically `condition`).

- formula:

  One-sided formula specifying the design (e.g., `~ condition`).

- alpha:

  Positive numeric pseudocount added to modified and unmodified read
  counts before log-transformation. Default `0.5`.

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

A `data.frame` with one row per site (same row order as `methyl_mat`),
containing:

- `pvalue`:

  Moderated t-test p-value from `eBayes`. `NA` for sites with any
  missing data.

- `delta_beta`:

  Effect size (treatment mean beta minus reference mean beta) on the 0–1
  scale. `NA` where group means cannot be computed.

- `mean_beta_<level>`:

  One column per condition level containing the per-group observed mean
  beta value.

## Details

limma must be installed (it is listed in `Suggests`). If it is not
available, this function stops with an informative message.

Beta values are first transformed to M-values: \$\$M = \log_2\\\left(
\frac{n\_{\mathrm{mod}} + \alpha}{n\_{\mathrm{unmod}} + \alpha}
\right)\$\$ where \\\alpha\\ is a pseudocount (default 0.5). M-values
are approximately normally distributed and homoscedastic, making OLS
appropriate.

A linear model is fitted across all sites simultaneously with
[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html).
[`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html) then estimates an
empirical Bayes prior on the residual variance across all sites and
computes a moderated posterior variance per site — shrinking the noisy
per-site estimate toward the genome-wide mean. P-values are derived from
a moderated t-statistic with posterior degrees of freedom \\d_0 +
df\_{\mathrm{residual}}\\.

Only sites where all samples have non-NA M-values (i.e., non-zero
coverage after `min_coverage` thresholding) are passed to limma. Sites
with any `NA` retain `NA` in all result columns.

Effect sizes (`delta_beta`) and per-group means are reported on the
original beta (0–1) scale for interpretability, not back-transformed
from M-value coefficients.
