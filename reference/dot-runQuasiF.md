# Per-site quasi-likelihood F-test for differential methylation

An internal wrapper that combines the per-site quasibinomial GLM of
`.betaBinomialTest` with empirical Bayes shrinkage of the per-site
dispersion estimates. Called by
[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
when `method = "quasi_f"`.

## Usage

``` r
.runQuasiF(
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

  Moderated quasi-F p-value. `NA` for untestable sites.

- `delta_beta`:

  Effect size (treatment mean beta minus reference mean beta) on the 0–1
  scale.

- `mean_beta_<level>`:

  One column per condition level containing the per-group observed mean
  beta value.

## Details

limma must be installed (it is listed in `Suggests`).

The method runs in three passes:

**Pass 1 — per-site GLM.** The same quasibinomial model as
`.betaBinomialTest` is fitted at each site:
\$\$\mathrm{glm}(\mathrm{cbind}(n\_{\mathrm{mod}},\\
n\_{\mathrm{unmod}}) \sim \mathrm{condition},\\ \mathrm{family} =
\mathrm{quasibinomial}())\$\$ For each site \\j\\, three quantities are
collected:

- \\\hat\phi_j = \\`fit\$dispersion` — Pearson chi-squared dispersion
  estimate

- \\df_j = \\`fit\$df.residual` — residual degrees of freedom

- \\\tilde{t}\_j^{(0)} = t_j \times \sqrt{\hat\phi_j}\\ — the "unscaled"
  t-statistic (independent of \\\hat\phi_j\\), where \\t_j\\ is the Wald
  t-statistic from `coef(summary(fit))`

**Pass 2 — empirical Bayes dispersion shrinkage.**
[`squeezeVar`](https://rdrr.io/pkg/limma/man/squeezeVar.html) pools the
\\\\\hat\phi_j\\\\ estimates across all testable sites, fits a
log-normal prior, and returns posterior dispersion estimates
\\\\\tilde\phi_j\\\\ and a prior degrees-of-freedom scalar \\d_0\\.

**Pass 3 — moderated test statistic.** The posterior t-statistic and
p-value for each site are: \$\$\tilde{t}\_j =
\frac{\tilde{t}\_j^{(0)}}{\sqrt{\tilde\phi_j}}, \quad p_j = 2\\P(T \leq
-\|\tilde{t}\_j\|),\\ T \sim t(d_0 + df_j)\$\$ The additional \\d_0\\
degrees of freedom are the power gain over the unadjusted quasibinomial
test.

This procedure is methodologically equivalent to the quasi-likelihood
F-test of edgeR (`glmQLFTest`), adapted for methylation proportions
(quasibinomial) rather than RNA-seq counts (quasi-negative-binomial).
