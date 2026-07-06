# Identify differentially methylated sites between conditions

The main function for differential methylation analysis in commaKit.
Analogous to `DESeq2::DESeq()`, `diffMethyl()` accepts a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object, fits a statistical model to each methylation site, and returns
the same object enriched with per-site test results stored as new
columns in `rowData`.

## Usage

``` r
diffMethyl(
  object,
  formula = ~condition,
  reference = NULL,
  method = c("methylkit", "limma", "quasi_f"),
  mod_context = NULL,
  mod_type = NULL,
  motif = NULL,
  min_coverage = 5L,
  alpha = 0.5,
  p_adjust_method = "BH",
  result_name = NULL,
  overwrite = NULL,
  make_default = TRUE,
  ...
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object with at least two samples in distinct conditions.

- formula:

  A one-sided formula specifying one two-level comparison variable, e.g.
  `~ condition`. The RHS variable must match a column in
  `sampleInfo(object)` and must have exactly 2 distinct levels.
  Multi-factor formulas, interactions, offsets, transformed terms, and
  continuous covariates are intentionally not interpreted in v1; subset
  or define one two-group comparison per `diffMethyl()` call. Default is
  `~ condition`.

- reference:

  Character string or `NULL`. The reference (control) level for the
  primary formula variable. When provided, it must match one of the two
  values present in the corresponding `colData` column. When `NULL`
  (default), the reference level is determined automatically: if the
  column is a factor, its first factor level is used; otherwise the
  alphabetically first value is used (matching R's default contrast
  behaviour). The reported `dm_delta_beta` is the non-reference level
  minus the reference level.

- method:

  Character string selecting the statistical backend. `"methylkit"`
  (default) wraps
  [`methylKit::calculateDiffMeth()`](https://rdrr.io/pkg/methylKit/man/calculateDiffMeth-methods.html)
  with logistic regression. The raw methylKit p-value feeds `dm_pvalue`;
  `dm_padj` is computed by commaKit using `p_adjust_method`; methylKit
  q-values are preserved separately as `dm_methylkit_qvalue`. Requires
  methylKit (Bioconductor). `"quasi_f"` applies empirical Bayes
  shrinkage of quasibinomial dispersions via
  [`squeezeVar`](https://rdrr.io/pkg/limma/man/squeezeVar.html)
  (quasi-likelihood F-test; count-data EB), making it a good
  general-purpose alternative for bacterial methylomes. Requires limma.
  `"limma"` applies empirical Bayes variance shrinkage via
  [`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html) on
  M-value-transformed data; recommended when replicates are few (n \< 3
  per group). Requires limma.

- mod_context:

  Character vector or `NULL`. Modification context(s) to test (e.g.,
  `"6mA_GATC"`, `c("6mA_GATC", "5mC_CCWGG")`). A `mod_context` value is
  `paste(mod_type, motif, sep = "_")` when motif information is
  available, or just `mod_type` for Dorado/Megalodon data. Use
  [`modContexts`](https://carl-stone.github.io/commaKit/reference/modContexts.md)`(object)`
  to see which contexts are present. If `NULL` (default), all contexts
  present in `object` are tested independently. When provided, takes
  precedence over the `mod_type` and `motif` arguments.

- mod_type:

  Character vector or `NULL`. Modification type(s) to test (e.g.,
  `"6mA"`, `c("6mA", "5mC")`). If `NULL` (default), all modification
  types present in `object` are tested. Ignored when `mod_context` is
  provided.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are tested (e.g., `"GATC"`). Uses
  [`motifs`](https://carl-stone.github.io/commaKit/reference/motifs.md)
  to validate the requested values. If `NULL` (default), all motifs
  (including `NA`) are included. Ignored when `mod_context` is provided.

- min_coverage:

  Integer. Minimum per-sample read depth required to include a site in
  testing. Sites where any sample has coverage below this threshold are
  treated as `NA` in that sample. Default `5L`.

- alpha:

  Positive numeric pseudocount used to compute M-values when
  `method = "limma"`: \\M = \log_2((n\_{\mathrm{mod}} + \alpha) /
  (n\_{\mathrm{unmod}} + \alpha))\\. Default `0.5` (a Beta(0.5, 0.5)
  prior). Ignored for other methods.

- p_adjust_method:

  Character string. Multiple testing correction method, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`
  (Benjamini-Hochberg). Other options: `"bonferroni"`, `"holm"`, `"BY"`,
  `"none"`.

- result_name:

  Character string or `NULL`. Name for the differential methylation
  result layer. If `NULL` (default), results are written to the
  compatibility layer `"diffMethyl"`, replacing that layer on repeated
  unnamed calls. Provide explicit names such as `"quasi_f.min5"` or
  `"limma.alpha05"` to keep multiple runs.

- overwrite:

  Logical or `NULL`. Whether an existing named result layer may be
  replaced. If `NULL` (default), unnamed compatibility runs may replace
  `"diffMethyl"`, while explicit `result_name` values must be unique
  unless `overwrite = TRUE`.

- make_default:

  Logical. If `TRUE` (default), the new result layer becomes the default
  used by
  [`results()`](https://carl-stone.github.io/commaKit/reference/results.md)
  and is mirrored into `rowData` as bare `dm_*` columns.

- ...:

  Additional arguments (reserved for future use).

## Value

The input `commaData` object with additional columns in `rowData`:
`dm_pvalue`, `dm_padj`, `dm_delta_beta`, and one
`dm_mean_beta_<condition>` column per condition level. The methylKit
backend also adds `dm_methylkit_qvalue`. The `metadata` slot is updated
with analysis parameters, result column names, and a named result layer
registry.

## Details

**Default method (`method = "methylkit"`):** Wraps
[`methylKit::calculateDiffMeth()`](https://rdrr.io/pkg/methylKit/man/calculateDiffMeth-methods.html),
which uses logistic regression and returns raw p-values plus
methylKit-specific q-values. commaKit stores the raw methylKit p-value
as `dm_pvalue`, computes `dm_padj` with the selected `p_adjust_method`
across the full `diffMethyl()` testing family, and preserves methylKit
q-values separately as `dm_methylkit_qvalue`. Requires methylKit
(`BiocManager::install("methylKit")`).

**Alternative model (`method = "quasi_f"`):** A per-site quasibinomial
GLM with empirical Bayes shrinkage of dispersion estimates, analogous to
the quasi-likelihood F-test of edgeR. Fits the quasibinomial GLM per
site, collects per-site dispersion \\\hat\phi_j\\ and residual df
\\df_j\\, then calls
[`squeezeVar`](https://rdrr.io/pkg/limma/man/squeezeVar.html) to compute
posterior dispersion estimates \\\\\tilde\phi_j\\\\. T-statistics are
evaluated against a t-distribution with \\d_0 + df_j\\ degrees of
freedom. Requires limma.

**Alternative model (`method = "limma"`):** Beta values are transformed
to M-values via \\M = \log_2((n\_{\mathrm{mod}} + \alpha) /
(n\_{\mathrm{unmod}} + \alpha))\\, then
[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) fits an OLS model
per site and [`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html)
applies empirical Bayes variance shrinkage, borrowing information across
all sites to stabilize the per-site variance estimate. Recommended when
replicates are few (n \< 3 per group). Requires limma
(`BiocManager::install("limma")`). Effect sizes are reported on the
original beta scale.

**Method selection guidance:** Use the default `"methylkit"` backend
when you want compatibility with established methylKit workflows or need
results that closely follow methylKit's logistic-regression conventions.
Use `"quasi_f"` as a good general-purpose alternative for bacterial
methylomes when you want a count-aware model with empirical Bayes
dispersion shrinkage and genome-wide multiple-testing correction handled
entirely inside commaKit. It is often a good first alternative if
methylKit convergence warnings, zero-variance sites, or runtime become
distracting. Use `"limma"` when you want the familiar limma
empirical-Bayes linear-model workflow on M-values, especially for
complete datasets with few replicates per group. All three backends
report `dm_delta_beta` on the original beta scale, so effect sizes
remain comparable even when p-values differ.

**Multiple mod contexts:** When `mod_context = NULL` (default), all
modification contexts (mod_type x motif combinations) present in the
object are tested independently and results are combined. Sites not
belonging to a tested context receive `NA` in all `dm_*` columns.
Testing by `mod_context` rather than just `mod_type` prevents spurious
pooling of biologically distinct methylation events (e.g., 6mA at GATC
from Dam methyltransferase versus cytosine modifications at GATC, which
are likely artefactual).

**Small-sample note:** Differential methylation testing with very few
replicates (e.g., n = 1 per group) is mathematically possible but has
extremely low statistical power. Treat such results as exploratory only.

**Result columns added to `rowData`:**

- `dm_pvalue`:

  Raw p-value from the GLM Wald test.

- `dm_padj`:

  Adjusted p-value (Benjamini-Hochberg by default).

- `dm_methylkit_qvalue`:

  methylKit q-value, present only for result layers produced with
  `method = "methylkit"`.

- `dm_delta_beta`:

  Effect size: mean methylation in the treatment group minus mean
  methylation in the reference (control) group. Positive values indicate
  hypermethylation in treatment.

- `dm_mean_beta_<condition>`:

  One column per condition level (named after the actual condition
  values), containing the per-group mean beta value for each site.

Analysis parameters and result column names are stored in
`metadata(object)$diffMethyl_params` and
`metadata(object)$diffMethyl_result_cols`. Named result layers are
stored in `metadata(object)$diffMethyl_results` and listed by
[`resultLayers()`](https://carl-stone.github.io/commaKit/reference/resultLayers.md).
By default, the active result layer is also mirrored into the legacy
`dm_*` columns in `rowData` so existing
[`results()`](https://carl-stone.github.io/commaKit/reference/results.md)
and plotting workflows keep working.

## See also

[`results`](https://carl-stone.github.io/commaKit/reference/results.md)
to extract the test results as a tidy `data.frame`;
[`filterResults`](https://carl-stone.github.io/commaKit/reference/filterResults.md)
to filter by significance thresholds.

## Examples

``` r
data(comma_example_data)
# Test for differential 6mA methylation between conditions
dm <- diffMethyl(comma_example_data, formula = ~condition, mod_type = "6mA")
#> diffMethyl: testing 'condition' -- 'treatment' vs 'control' (reference)
#> methylKit: comparing 'treatment' (treatment) vs 'control' (reference/control)
#> uniting...

# How many sites are significant?
rd <- as.data.frame(SummarizedExperiment::rowData(dm))
sum(rd$dm_padj < 0.05 & abs(rd$dm_delta_beta) >= 0.2, na.rm = TRUE)
#> [1] 31
```
