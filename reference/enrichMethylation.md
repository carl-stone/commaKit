# Gene set enrichment analysis of differential methylation results

Maps per-site differential methylation statistics up to the gene level
and runs GO and/or KEGG enrichment analysis using clusterProfiler.
Supports both over-representation analysis (ORA) and gene set enrichment
analysis (GSEA). Works with model organisms via standard Bioconductor
databases, or with any organism via a custom `TERM2GENE` mapping.

## Usage

``` r
enrichMethylation(
  object,
  method = "ora",
  OrgDb = NULL,
  keyType = "SYMBOL",
  ont = "BP",
  organism = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  kegg_term2gene = NULL,
  kegg_term2name = NULL,
  gene_col = "feature_names",
  feature_type = "gene",
  gene_role = c("target", "regulator", "both"),
  overlap_only = NULL,
  padj_threshold = 0.05,
  delta_beta_threshold = 0.1,
  score_metric = "combined",
  gene_score_agg = "max",
  mod_type = NULL,
  mod_context = NULL,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10L,
  maxGSSize = 500L
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object on which
  [`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
  and
  [`annotateSites`](https://carl-stone.github.io/commaKit/reference/annotateSites.md)
  have been run, **or** a `data.frame` produced by
  [`results()`](https://carl-stone.github.io/commaKit/reference/results.md).

- method:

  Character vector; one or both of `"ora"` (over- representation
  analysis) and `"gsea"` (gene set enrichment analysis). Default
  `"ora"`.

- OrgDb:

  A Bioconductor `OrgDb` annotation object for GO enrichment (e.g.,
  `org.EcK12.eg.db`). `NULL` skips GO analysis unless `TERM2GENE` is
  provided.

- keyType:

  Character string; key type for gene IDs in `OrgDb` (e.g., `"SYMBOL"`,
  `"ENTREZID"`). Default `"SYMBOL"`.

- ont:

  Character string; GO ontology. One of `"BP"`, `"MF"`, `"CC"`, or
  `"ALL"`. Default `"BP"`.

- organism:

  Character string; KEGG organism code (e.g., `"eco"`). `NULL` skips
  KEGG. Requires internet access and fires one HTTP request per KEGG
  pathway, which can exceed the API rate limit for organisms with many
  pathways. Prefer `kegg_term2gene` (built once via
  [`buildKEGGTermGene()`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md))
  for reliable KEGG analysis.

- TERM2GENE:

  A two-column `data.frame` with columns `term` and `gene`. When
  provided, takes precedence over `OrgDb` for GO analysis; results
  appear in the `$go` slot. To supply a pre-built KEGG mapping, use
  `kegg_term2gene` instead so that results land in the `$kegg` slot.

- TERM2NAME:

  Optional two-column `data.frame` with columns `term` and `name`.
  Paired with `TERM2GENE` for GO enrichment.

- kegg_term2gene:

  A two-column `data.frame` with columns `term` and `gene` containing
  pre-fetched KEGG pathway-gene associations. Build it once per organism
  with
  [`buildKEGGTermGene()`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md)
  and reuse across analyses — no live KEGG API calls are made. Results
  appear in the `$kegg` slot, consistent with the `organism` path. Takes
  precedence over `organism` when both are supplied.

- kegg_term2name:

  Optional two-column `data.frame` with columns `term` and `name`
  containing KEGG pathway descriptions. Returned by
  [`buildKEGGTermGene()`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md)
  alongside `kegg_term2gene`.

- gene_col:

  Character string; the `rowData` column containing gene identifiers per
  site. Default `"feature_names"`.

- feature_type:

  Character vector or `NULL`. Feature type(s) to analyse. When more than
  one type is given, each is run separately and results are returned as
  a named list. Default `"gene"`. Set `NULL` to include all annotated
  sites.

- gene_role:

  Character string; which role to test. One of `"target"` (default),
  `"regulator"`, or `"both"`. See Details.

- overlap_only:

  Logical, `NULL`, or a named logical vector. When `TRUE`, only sites
  where `rel_position == 0` (inside the feature) contribute. `NULL`
  (default) auto-detects: `TRUE` for large region features (gene, CDS,
  etc.), `FALSE` for small regulatory features (promoter, binding sites,
  etc.).

- padj_threshold:

  Numeric; ORA significance threshold. Default `0.05`.

- delta_beta_threshold:

  Numeric; minimum \\\|\Delta\beta\|\\ for ORA. Default `0.1`.

- score_metric:

  Character string; GSEA ranking metric. One of `"combined"` (default),
  `"padj"`, or `"delta_beta"`.

- gene_score_agg:

  Character string; aggregation across sites per gene. Either `"max"`
  (default) or `"mean"`.

- mod_type:

  Character vector or `NULL`; modification-type filter passed to
  [`results`](https://carl-stone.github.io/commaKit/reference/results.md).
  Ignored when `object` is a `data.frame`.

- mod_context:

  Character string or `NULL`; mod-context filter passed to
  [`results`](https://carl-stone.github.io/commaKit/reference/results.md).
  Ignored for `data.frame` input.

- pvalueCutoff:

  Numeric; p-value cutoff for clusterProfiler. Default `0.05`.

- qvalueCutoff:

  Numeric; q-value cutoff (ORA only). Default `0.2`.

- minGSSize:

  Integer; minimum gene set size. Default `10`.

- maxGSSize:

  Integer; maximum gene set size. Default `500`.

## Value

When a single `feature_type` and `gene_role != "both"` are requested: a
named list `list(go = ..., kegg = ...)` (backward-compatible). Each
element is a clusterProfiler `enrichResult` (ORA) or `gseaResult`
(GSEA), or `NULL`. When both ORA and GSEA are requested, each element is
a nested list `list(ora = ..., gsea = ...)`. When multiple
`feature_type` values are given, a named list per feature type is
returned. When `gene_role = "both"`, each feature type's result is a
list `list(target = ..., regulator = ...)`.

## Details

Different feature types ask different biological questions, and the
appropriate comparison universe changes accordingly. Use `gene_role` to
specify which perspective to test:

- `"target"`:

  (default) Which target genes have methylated sites in/near this
  feature type? Universe = all genes with site coverage.

- `"regulator"`:

  Which regulators (sigma factors, TFs, RNA regulators) have methylated
  binding sites? Universe = all regulators of the same type found in the
  annotation for this feature type.

- `"both"`:

  Run target and regulator enrichments separately, returning a named
  sub-list `list(target=..., regulator=...)`.

## Prerequisites

Before calling `enrichMethylation()`, you must:

1.  Run
    [`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
    to compute per-site `dm_padj` and `dm_delta_beta` values.

2.  Run
    [`annotateSites`](https://carl-stone.github.io/commaKit/reference/annotateSites.md)
    to assign feature identifiers to each site. For regulatory feature
    types, pass
    `metadata_cols = c("feature_subtype", "transcription_unit")` to
    capture sigma factor identity and target gene information.

## Gene-to-pathway mapping

Supply at least one of the following:

- `kegg_term2gene`:

  (Recommended for KEGG) A two-column `data.frame` pre-built by
  [`buildKEGGTermGene()`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md).
  No network access is required at analysis time; results appear in the
  `$kegg` slot.

- `TERM2GENE`:

  A two-column `data.frame` for custom GO enrichment. Results appear in
  the `$go` slot.

- `OrgDb`:

  A Bioconductor `OrgDb` annotation object (e.g., `org.EcK12.eg.db`) for
  GO enrichment. Ignored when `TERM2GENE` is supplied.

- `organism`:

  A KEGG organism code (e.g., `"eco"`). Makes one live HTTP request per
  KEGG pathway; may exceed the API rate limit. Ignored when
  `kegg_term2gene` is provided.

## Feature-type-specific parsing

Gene names are extracted from the annotation using feature-type-aware
rules:

- `"gene"`, CDS, etc.:

  Feature name is used directly.

- `"promoter"`, `"minus_10_signal"`, `"minus_35_signal"`,
  `"transcription_factor_binding_site"`:

  Strip trailing `p` or `p`\\N\\ suffix, then split on `"-"` for
  operonic promoters.

- `"protein_binding_site"`:

  Target gene from `transcription_unit` attribute (if available via
  `metadata_cols`) or from *"... of {promoter}"* pattern. Regulator gene
  from binding-protein name.

- `"RNA_binding_site"`:

  Target from `transcription_unit` or *"regulating {gene}"* pattern.
  Regulator from first word if it matches a gene-name pattern.

- `"terminator"`:

  Strip *" terminator"* suffix.

## GSEA ranking

When `method` includes `"gsea"`, a per-gene score is computed from the
site-level statistics and genes are ranked in decreasing order. Three
metrics are available via `score_metric`:

- `"combined"`:

  (default) \\-\log\_{10}(\text{padj}) \times
  \text{sign}(\Delta\beta)\\; captures both significance and direction.

- `"padj"`:

  \\-\log\_{10}(\text{padj})\\; significance only, direction-agnostic.

- `"delta_beta"`:

  \\\Delta\beta\\ directly; effect-size ranking.

When multiple sites map to the same gene, the score is aggregated across
sites using `gene_score_agg`: `"max"` (default; selects the site with
the largest absolute score, preserving its sign) or `"mean"`.

## See also

[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md),
[`annotateSites`](https://carl-stone.github.io/commaKit/reference/annotateSites.md),
[`results`](https://carl-stone.github.io/commaKit/reference/results.md),
[`filterResults`](https://carl-stone.github.io/commaKit/reference/filterResults.md)

## Examples

``` r
# \donttest{
# Requires clusterProfiler and a custom TERM2GENE mapping
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  data(comma_example_data)
  dm <- diffMethyl(comma_example_data, formula = ~condition, mod_type = "6mA")
  ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")

  # Custom TERM2GENE (works without network access or OrgDb)
  fake_t2g <- data.frame(
    term = c("PATH:01", "PATH:01", "PATH:02"),
    gene = c("geneA", "geneB", "geneC")
  )
  res <- enrichMethylation(ann, TERM2GENE = fake_t2g, method = c("ora", "gsea"))
  str(res, max.level = 2)
}
#> diffMethyl: testing 'condition' -- 'treatment' vs 'control' (reference)
#> methylKit: comparing 'treatment' (treatment) vs 'control' (reference/control)
#> uniting...
#> No gene sets have size between 10 and 500 ...
#> --> return NULL...
#> using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
#> preparing geneSet collections...
#> GSEA analysis...
#> no term enriched under specific pvalueCutoff...
#> List of 2
#>  $ go  :List of 2
#>   ..$ ora : NULL
#>   ..$ gsea:Formal class 'gseaResult' [package "DOSE"] with 13 slots
#>  $ kegg:List of 2
#>   ..$ ora : NULL
#>   ..$ gsea: NULL
# }
```
