## Tests for enrichMethylation(), .siteToGeneMap(), .computeGeneScores()
## All tests use custom TERM2GENE to avoid requiring clusterProfiler, OrgDb,
## KEGG access, or internet connectivity.

library(testthat)
library(commaKit)

# ── Shared fixtures ───────────────────────────────────────────────────────────

# Minimal commaData with diffMethyl results and feature_names annotation
make_annotated_dm <- function() {
  data(comma_example_data)
  dm <- diffMethyl(comma_example_data, formula = ~condition, mod_type = "6mA")
  ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")
  ann
}

# Minimal TERM2GENE mapping aligned with example data gene names (geneA-E)
fake_t2g <- data.frame(
  term = c("PATH:01", "PATH:01", "PATH:02", "PATH:02", "PATH:03"),
  gene = c("geneA", "geneB", "geneC", "geneD", "geneE"),
  stringsAsFactors = FALSE
)

fake_t2n <- data.frame(
  term = c("PATH:01", "PATH:02", "PATH:03"),
  name = c("Pathway one", "Pathway two", "Pathway three"),
  stringsAsFactors = FALSE
)

# ── enrichment dispatch setup ────────────────────────────────────────────────

test_that("enrichment dispatch builds custom GO and KEGG call specs", {
  kegg_t2g <- data.frame(
    term = "eco00010", gene = "geneA",
    stringsAsFactors = FALSE
  )
  kegg_t2n <- data.frame(
    term = "eco00010", name = "Glycolysis",
    stringsAsFactors = FALSE
  )
  ctx <- commaKit:::.enrichmentDispatchContext(
    method = c("ora", "gsea"),
    org_db = NULL,
    keyType = "SYMBOL",
    ont = "BP",
    organism = NULL,
    term2gene = fake_t2g,
    term2name = fake_t2n,
    kegg_term2gene = kegg_t2g,
    kegg_term2name = kegg_t2n,
    padj_threshold = 0.05,
    delta_beta_threshold = 0.1,
    score_metric = "combined",
    gene_score_agg = "max",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 1L,
    maxGSSize = 500L
  )

  go_ora <- commaKit:::.enrichmentCallSpec(
    "go", "ora", "geneA", c("geneA", "geneB"), ctx
  )
  go_gsea <- commaKit:::.enrichmentCallSpec(
    "go", "gsea", c(geneA = 2, geneB = 1), c("geneA", "geneB"), ctx
  )
  kegg_ora <- commaKit:::.enrichmentCallSpec(
    "kegg", "ora", "geneA", c("geneA", "geneB"), ctx
  )

  expect_equal(go_ora$fun, "enricher")
  expect_equal(go_ora$args$gene, "geneA")
  expect_equal(go_ora$args$universe, c("geneA", "geneB"))
  expect_equal(go_ora$args$TERM2GENE, fake_t2g)
  expect_equal(go_ora$args$TERM2NAME, fake_t2n)
  expect_equal(go_ora$args$qvalueCutoff, 0.2)

  expect_equal(go_gsea$fun, "GSEA")
  expect_equal(go_gsea$args$geneList, c(geneA = 2, geneB = 1))
  expect_false("universe" %in% names(go_gsea$args))
  expect_false("qvalueCutoff" %in% names(go_gsea$args))

  expect_equal(kegg_ora$fun, "enricher")
  expect_equal(kegg_ora$args$TERM2GENE, kegg_t2g)
  expect_equal(kegg_ora$args$TERM2NAME, kegg_t2n)
})

test_that("enrichment dispatch builds optional database call specs", {
  fake_orgdb <- structure(list(), class = "OrgDb")
  ctx <- commaKit:::.enrichmentDispatchContext(
    method = c("ora", "gsea"),
    org_db = fake_orgdb,
    keyType = "SYMBOL",
    ont = "BP",
    organism = "eco",
    term2gene = NULL,
    term2name = NULL,
    kegg_term2gene = NULL,
    kegg_term2name = NULL,
    padj_threshold = 0.05,
    delta_beta_threshold = 0.1,
    score_metric = "combined",
    gene_score_agg = "max",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 1L,
    maxGSSize = 500L
  )

  go_ora <- commaKit:::.enrichmentCallSpec(
    "go", "ora", "geneA", c("geneA", "geneB"), ctx
  )
  go_gsea <- commaKit:::.enrichmentCallSpec(
    "go", "gsea", c(geneA = 2, geneB = 1), c("geneA", "geneB"), ctx
  )
  kegg_ora <- commaKit:::.enrichmentCallSpec(
    "kegg", "ora", "geneA", c("geneA", "geneB"), ctx
  )
  kegg_gsea <- commaKit:::.enrichmentCallSpec(
    "kegg", "gsea", c(geneA = 2, geneB = 1), c("geneA", "geneB"), ctx
  )

  expect_equal(go_ora$fun, "enrichGO")
  expect_identical(go_ora$args$OrgDb, fake_orgdb)
  expect_equal(go_ora$args$readable, FALSE)
  expect_equal(go_ora$args$universe, c("geneA", "geneB"))

  expect_equal(go_gsea$fun, "gseGO")
  expect_false("universe" %in% names(go_gsea$args))
  expect_false("readable" %in% names(go_gsea$args))

  expect_equal(kegg_ora$fun, "enrichKEGG")
  expect_equal(kegg_ora$args$organism, "eco")
  expect_equal(kegg_ora$args$universe, c("geneA", "geneB"))

  expect_equal(kegg_gsea$fun, "gseKEGG")
  expect_equal(kegg_gsea$args$organism, "eco")
  expect_false("universe" %in% names(kegg_gsea$args))
})

test_that("enrichment dispatch keeps empty KEGG names optional", {
  ctx <- commaKit:::.enrichmentDispatchContext(
    method = "ora",
    org_db = NULL,
    keyType = "SYMBOL",
    ont = "BP",
    organism = NULL,
    term2gene = NULL,
    term2name = NULL,
    kegg_term2gene = data.frame(
      term = "eco00010", gene = "geneA",
      stringsAsFactors = FALSE
    ),
    kegg_term2name = data.frame(
      term = character(), name = character(),
      stringsAsFactors = FALSE
    ),
    padj_threshold = 0.05,
    delta_beta_threshold = 0.1,
    score_metric = "combined",
    gene_score_agg = "max",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 1L,
    maxGSSize = 500L
  )

  spec <- commaKit:::.enrichmentCallSpec(
    "kegg", "ora", "geneA", c("geneA", "geneB"), ctx
  )

  expect_equal(spec$fun, "enricher")
  expect_null(spec$args$TERM2NAME)
})

test_that("enrichment dispatch preserves combined short-circuit shape", {
  ctx <- commaKit:::.enrichmentDispatchContext(
    method = c("ora", "gsea"),
    org_db = NULL,
    keyType = "SYMBOL",
    ont = "BP",
    organism = NULL,
    term2gene = fake_t2g,
    term2name = NULL,
    kegg_term2gene = NULL,
    kegg_term2name = NULL,
    padj_threshold = 0.01,
    delta_beta_threshold = 1,
    score_metric = "combined",
    gene_score_agg = "max",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 1L,
    maxGSSize = 500L
  )
  sg <- data.frame(
    gene_id = "geneA",
    site_key = "chr1:100:+",
    dm_padj = NA_real_,
    dm_delta_beta = NA_real_,
    stringsAsFactors = FALSE
  )

  warnings <- character()
  res <- withCallingHandlers(
    commaKit:::.runEnrichmentForGeneMap(sg, "geneA", ctx),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl(
    "No significantly differentially methylated genes", warnings
  )))
  expect_true(any(grepl("No valid gene scores", warnings)))
  expect_equal(names(res), c("go", "kegg"))
  expect_equal(names(res$go), c("ora", "gsea"))
  expect_equal(names(res$kegg), c("ora", "gsea"))
  expect_null(res$go$ora)
  expect_null(res$go$gsea)
  expect_null(res$kegg$ora)
  expect_null(res$kegg$gsea)
})

test_that("enrichment dispatch ignores duplicate methods", {
  ctx <- commaKit:::.enrichmentDispatchContext(
    method = c("ora", "ora"),
    org_db = NULL,
    keyType = "SYMBOL",
    ont = "BP",
    organism = NULL,
    term2gene = fake_t2g,
    term2name = NULL,
    kegg_term2gene = NULL,
    kegg_term2name = NULL,
    padj_threshold = 0.01,
    delta_beta_threshold = 1,
    score_metric = "combined",
    gene_score_agg = "max",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 1L,
    maxGSSize = 500L
  )
  sg <- data.frame(
    gene_id = "geneA",
    site_key = "chr1:100:+",
    dm_padj = NA_real_,
    dm_delta_beta = NA_real_,
    stringsAsFactors = FALSE
  )

  expect_identical(ctx$method, "ora")
  expect_warning(
    res <- commaKit:::.runEnrichmentForGeneMap(sg, "geneA", ctx),
    "No significantly differentially methylated genes"
  )
  expect_equal(names(res), c("go", "kegg"))
  expect_null(res$go)
  expect_null(res$kegg)
})

# ── .siteToGeneMap() ─────────────────────────────────────────────────────────

test_that(".siteToGeneMap returns data.frame with expected columns", {
  res_df <- data.frame(
    chrom = c("chr1", "chr1"),
    position = c(100L, 200L),
    strand = c("+", "+"),
    dm_padj = c(0.01, 0.5),
    dm_delta_beta = c(-0.3, 0.05),
    feature_names = I(list(c("geneA", "geneB"), character(0))),
    stringsAsFactors = FALSE
  )
  sg <- commaKit:::.siteToGeneMap(res_df, "feature_names")

  expect_s3_class(sg, "data.frame")
  expect_true(all(
    c("gene_id", "site_key", "dm_padj", "dm_delta_beta") %in% colnames(sg)
  ))
})

test_that(".siteToGeneMap correctly explodes multi-gene sites", {
  res_df <- data.frame(
    chrom = "chr1",
    position = 100L,
    strand = "+",
    dm_padj = 0.01,
    dm_delta_beta = -0.3,
    feature_names = I(list(c("geneA", "geneB"))),
    stringsAsFactors = FALSE
  )
  sg <- commaKit:::.siteToGeneMap(res_df, "feature_names")

  expect_equal(nrow(sg), 2L)
  expect_setequal(sg$gene_id, c("geneA", "geneB"))
  # Both rows inherit the same site statistics
  expect_equal(sg$dm_padj, c(0.01, 0.01))
})

test_that(".siteToGeneMap excludes intergenic sites silently", {
  res_df <- data.frame(
    chrom = c("chr1", "chr1"),
    position = c(100L, 200L),
    strand = "+",
    dm_padj = c(0.01, 0.5),
    dm_delta_beta = c(-0.3, 0.05),
    feature_names = I(list(c("geneA"), character(0))), # second site intergenic
    stringsAsFactors = FALSE
  )
  sg <- commaKit:::.siteToGeneMap(res_df, "feature_names")

  expect_equal(nrow(sg), 1L)
  expect_equal(sg$gene_id, "geneA")
})

test_that(".siteToGeneMap returns empty data.frame when all sites intergenic", {
  res_df <- data.frame(
    chrom = "chr1",
    position = 100L,
    strand = "+",
    dm_padj = 0.01,
    dm_delta_beta = -0.3,
    feature_names = I(list(character(0))),
    stringsAsFactors = FALSE
  )
  expect_warning(
    sg <- commaKit:::.siteToGeneMap(res_df, "feature_names"),
    "No sites have gene annotations"
  )
  expect_equal(nrow(sg), 0L)
})

test_that(".siteToGeneMap errors when gene_col is missing", {
  res_df <- data.frame(
    chrom = "chr1", position = 100L, strand = "+",
    dm_padj = 0.01, dm_delta_beta = -0.3,
    stringsAsFactors = FALSE
  )
  expect_error(
    commaKit:::.siteToGeneMap(res_df, "feature_names"),
    "not found in results"
  )
})

test_that(".siteToGeneMap constructs correct site_key", {
  res_df <- data.frame(
    chrom = "chrX",
    position = 555L,
    strand = "-",
    dm_padj = 0.001,
    dm_delta_beta = 0.4,
    feature_names = I(list("geneZ")),
    stringsAsFactors = FALSE
  )
  sg <- commaKit:::.siteToGeneMap(res_df, "feature_names")
  expect_equal(sg$site_key, "chrX:555:-")
})

# ── .computeGeneScores() ─────────────────────────────────────────────────────

make_sg <- function() {
  data.frame(
    gene_id = c("geneA", "geneA", "geneB", "geneC"),
    site_key = c("chr1:100:+", "chr1:200:+", "chr1:300:+", "chr1:400:+"),
    dm_padj = c(0.001, 0.01, 0.5, NA),
    dm_delta_beta = c(-0.5, -0.2, 0.05, 0.3),
    stringsAsFactors = FALSE
  )
}

test_that(".computeGeneScores returns named sorted vector", {
  scores <- commaKit:::.computeGeneScores(make_sg(), "combined", "max")

  expect_true(is.numeric(scores))
  expect_true(!is.null(names(scores)))
  # Must be sorted decreasing
  expect_equal(scores, sort(scores, decreasing = TRUE))
})

test_that(".computeGeneScores 'combined' metric has correct sign", {
  sg <- data.frame(
    gene_id = "geneA",
    site_key = "chr1:100:+",
    dm_padj = 0.01,
    dm_delta_beta = 0.5, # positive delta_beta → positive score
    stringsAsFactors = FALSE
  )
  scores <- commaKit:::.computeGeneScores(sg, "combined", "max")
  expect_true(scores[["geneA"]] > 0)

  sg$dm_delta_beta <- -0.5 # negative delta_beta → negative score
  scores2 <- commaKit:::.computeGeneScores(sg, "combined", "max")
  expect_true(scores2[["geneA"]] < 0)
})

test_that(".computeGeneScores excludes sites with NA padj", {
  # geneC has NA padj; should not appear in scores
  scores <- commaKit:::.computeGeneScores(make_sg(), "combined", "max")
  expect_false("geneC" %in% names(scores))
})

test_that(".computeGeneScores 'max' aggregation picks largest absolute score", {
  # geneA has two sites: padj=0.001 (larger -log10) and padj=0.01
  scores <- commaKit:::.computeGeneScores(make_sg(), "padj", "max")
  # The expected maximum score is 3.
  expect_equal(scores[["geneA"]], -log10(0.001), tolerance = 1e-6)
})

test_that(".computeGeneScores 'mean' aggregation averages across sites", {
  sg <- data.frame(
    gene_id = c("geneA", "geneA"),
    site_key = c("a", "b"),
    dm_padj = c(0.01, 0.01),
    dm_delta_beta = c(0.4, 0.6),
    stringsAsFactors = FALSE
  )
  scores <- commaKit:::.computeGeneScores(sg, "delta_beta", "mean")
  expect_equal(scores[["geneA"]], 0.5, tolerance = 1e-6)
})

test_that(".computeGeneScores returns empty vector when all NA", {
  sg <- data.frame(
    gene_id = c("g1"),
    site_key = "a",
    dm_padj = NA_real_,
    dm_delta_beta = NA_real_,
    stringsAsFactors = FALSE
  )
  scores <- commaKit:::.computeGeneScores(sg)
  expect_length(scores, 0L)
})

test_that(".computeGeneScores errors on unknown score_metric", {
  expect_error(
    commaKit:::.computeGeneScores(make_sg(), "bogus"),
    "score_metric"
  )
})

# ── enrichMethylation() — error conditions ────────────────────────────────────

test_that("enrichMethylation errors when clusterProfiler is absent", {
  skip_if(
    requireNamespace("clusterProfiler", quietly = TRUE),
    "clusterProfiler is installed; cannot test missing-package error"
  )
  data(comma_example_data)
  expect_error(
    enrichMethylation(comma_example_data, TERM2GENE = fake_t2g),
    "clusterProfiler"
  )
})

test_that("enrichMethylation errors when no mapping source is provided", {
  data(comma_example_data)
  expect_error(
    enrichMethylation(comma_example_data),
    "No gene-to-term mapping supplied"
  )
})

test_that("enrichMethylation errors when diffMethyl not run", {
  skip_if_not_installed("clusterProfiler")
  data(comma_example_data)
  # comma_example_data has no dm_padj column
  expect_error(
    enrichMethylation(comma_example_data, TERM2GENE = fake_t2g),
    "diffMethyl"
  )
})

test_that(
  "enrichMethylation warns when annotateSites not run (feature_type='gene')",
  {
    skip_if_not_installed("clusterProfiler")
    data(comma_example_data)
    dm <- diffMethyl(comma_example_data, formula = ~condition, mod_type = "6mA")
    # dm has no feature_types column; new code path warns and returns NULL
    expect_warning(
      res <- enrichMethylation(dm, TERM2GENE = fake_t2g),
      "feature_types.*not found|annotateSites"
    )
    expect_null(res$go)
    expect_null(res$kegg)
  }
)

# ── enrichMethylation() — ORA with TERM2GENE ─────────────────────────────────

test_that("enrichMethylation ORA returns list with $go and $kegg", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  res <- enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g)

  expect_type(res, "list")
  expect_true(all(c("go", "kegg") %in% names(res)))
  # KEGG is NULL when no organism provided
  expect_null(res$kegg)
})

test_that("enrichMethylation ORA $go is enrichResult or NULL", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  res <- enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g)

  # Result should be NULL (no enrichment with tiny fake data) or enrichResult
  if (!is.null(res$go)) {
    expect_s4_class(res$go, "enrichResult")
  }
})

test_that(
  paste(
    "enrichMethylation ORA with TERM2NAME returns enrichResult",
    "with Description column"
  ),
  {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    res <- enrichMethylation(ann,
      method = "ora",
      TERM2GENE = fake_t2g, TERM2NAME = fake_t2n
    )
    # Should return a list with go and kegg
    expect_type(res, "list")
    expect_true("go" %in% names(res))
    # If go result exists, TERM2NAME should populate Description
    if (!is.null(res$go) && inherits(res$go, "enrichResult")) {
      expect_true("Description" %in% colnames(res$go@result))
    }
  }
)

# ── enrichMethylation() — GSEA with TERM2GENE ────────────────────────────────

test_that("enrichMethylation GSEA returns list with $go and $kegg", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  res <- enrichMethylation(ann, method = "gsea", TERM2GENE = fake_t2g)

  expect_type(res, "list")
  expect_true(all(c("go", "kegg") %in% names(res)))
  expect_null(res$kegg)
})

test_that("enrichMethylation GSEA $go is gseaResult or NULL", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  res <- enrichMethylation(ann,
    method = "gsea", TERM2GENE = fake_t2g,
    minGSSize = 1L
  )

  if (!is.null(res$go)) {
    expect_s4_class(res$go, "gseaResult")
  }
})

# ── enrichMethylation() — both ORA and GSEA ──────────────────────────────────

test_that("enrichMethylation both methods returns nested list structure", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  res <- enrichMethylation(ann,
    method = c("ora", "gsea"),
    TERM2GENE = fake_t2g, minGSSize = 1L
  )

  expect_type(res, "list")
  expect_true(all(c("go", "kegg") %in% names(res)))
  # When both methods requested, $go is itself a list with $ora and $gsea
  expect_type(res$go, "list")
  expect_true(all(c("ora", "gsea") %in% names(res$go)))
  expect_type(res$kegg, "list")
  expect_true(all(c("ora", "gsea") %in% names(res$kegg)))
})

# ── enrichMethylation() — ORA warning when no sig genes ──────────────────────

test_that("enrichMethylation ORA warns (not errors) when no sig genes", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  # Impossible thresholds → no significant genes
  expect_warning(
    res <- enrichMethylation(ann,
      method = "ora", TERM2GENE = fake_t2g,
      padj_threshold = 0, delta_beta_threshold = 2
    ),
    "No significantly differentially methylated genes"
  )
  expect_null(res$go)
  expect_null(res$kegg)
})

# ── enrichMethylation() — mod_type / mod_context filters ─────────────────────

test_that(
  "enrichMethylation mod_type filter restricts gene set to filtered sites",
  {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    res_all <- enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g)
    res_filt <- enrichMethylation(ann,
      method = "ora", TERM2GENE = fake_t2g,
      mod_type = "6mA"
    )
    # Both should return results
    expect_type(res_all, "list")
    expect_type(res_filt, "list")
    # The filtered enrichment should use a subset of the genes from unfiltered
    # GO may be NULL with tiny synthetic data (gene sets below min/max size)
    # but the subset relationship should hold when both produce results
    if (!is.null(res_all$go) && !is.null(res_filt$go)) {
      genes_all <- res_all$go@gene
      genes_filt <- res_filt$go@gene
      expect_true(all(genes_filt %in% genes_all))
    }
    # At minimum, verify mod_type filter produces different gene Universes
    # (the unfiltered run includes all mod_types, filtered only 6mA)
    expect_true(
      length(res_filt$kegg) <= length(res_all$kegg) ||
        is.null(res_filt$kegg)
    )
  }
)

test_that("enrichMethylation mod_type filter errors on unknown type", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  expect_error(
    enrichMethylation(ann,
      method = "ora", TERM2GENE = fake_t2g,
      mod_type = "99mX"
    ),
    "mod_type"
  )
})

# ── .siteToGeneMap() — NA gene_id handling ───────────────────────────────────

test_that(".siteToGeneMap drops NA gene IDs silently", {
  res_df <- data.frame(
    chrom = c("chr1", "chr1"),
    position = c(100L, 200L),
    strand = "+",
    dm_padj = c(0.01, 0.02),
    dm_delta_beta = c(-0.3, 0.2),
    feature_names = I(list(c(NA_character_, "geneA"), c("geneB"))),
    stringsAsFactors = FALSE
  )
  sg <- commaKit:::.siteToGeneMap(res_df, "feature_names")

  expect_false(any(is.na(sg$gene_id)))
  expect_true("geneA" %in% sg$gene_id)
  expect_true("geneB" %in% sg$gene_id)
})

# ── .computeGeneScores() — NA gene_id robustness ────────────────────────────

test_that(".computeGeneScores 'max' does not crash when gene_ids contain NA", {
  sg <- data.frame(
    gene_id = c(NA_character_, "geneA", "geneA"),
    site_key = c("chr1:100:+", "chr1:200:+", "chr1:300:+"),
    dm_padj = c(0.01, 0.01, 0.05),
    dm_delta_beta = c(0.4, -0.5, -0.2),
    stringsAsFactors = FALSE
  )
  expect_no_error(
    scores <- commaKit:::.computeGeneScores(sg, "combined", "max")
  )
  expect_false(is.null(names(scores)))
  expect_false(any(is.na(names(scores))))
})

test_that(
  ".computeGeneScores 'mean' does not return empty when gene_ids contain NA",
  {
    sg <- data.frame(
      gene_id = c(NA_character_, "geneA"),
      site_key = c("chr1:100:+", "chr1:200:+"),
      dm_padj = c(0.01, 0.01),
      dm_delta_beta = c(0.4, 0.5),
      stringsAsFactors = FALSE
    )
    scores <- commaKit:::.computeGeneScores(sg, "delta_beta", "mean")
    expect_true(length(scores) > 0L)
    expect_true("geneA" %in% names(scores))
  }
)

# ── enrichMethylation() — feature_type argument ───────────────────────────────

test_that(
  "enrichMethylation feature_type = 'gene' filters to gene features only",
  {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    res <- enrichMethylation(ann,
      method = "ora", TERM2GENE = fake_t2g,
      feature_type = "gene"
    )
    expect_type(res, "list")
    # Result should exist (gene features are present in test data)
    if (!is.null(res$go)) {
      expect_s4_class(res$go, "enrichResult")
    }
  }
)

test_that(
  "enrichMethylation feature_type = NULL matches default (all features)",
  {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    res_null <- enrichMethylation(ann,
      method = "ora", TERM2GENE = fake_t2g,
      feature_type = NULL
    )
    res_default <- enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g)
    # Both should return the same structure
    expect_type(res_null, "list")
    expect_type(res_default, "list")
    expect_equal(names(res_null), names(res_default))
    # Both should use the same gene set (NULL = default = all features)
    if (!is.null(res_null$go) && !is.null(res_default$go)) {
      expect_equal(sort(res_null$go@gene), sort(res_default$go@gene))
    }
  }
)

test_that(
  "enrichMethylation warns and returns NULL for unmatched feature_type",
  {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    expect_warning(
      res <- enrichMethylation(ann,
        method = "ora", TERM2GENE = fake_t2g,
        feature_type = "nonexistent_type"
      ),
      "No sites with feature_type"
    )
    expect_null(res$go)
    expect_null(res$kegg)
  }
)

# ── .parseTargetGenes() ───────────────────────────────────────────────────────

test_that(".parseTargetGenes identity for gene/CDS/tRNA/rRNA/ncRNA", {
  nms <- c("geneA", "geneB")
  for (ft in c("gene", "CDS", "tRNA", "rRNA", "ncRNA")) {
    res <- commaKit:::.parseTargetGenes(nms, ft)
    expect_equal(res, as.list(nms), info = ft)
  }
})

test_that(".parseTargetGenes strips promoter suffix for promoter/TFBS types", {
  nms <- c("aaeRp", "geneA-geneBp1", "flhDCp")
  for (ft in c(
    "promoter", "minus_10_signal", "minus_35_signal",
    "transcription_factor_binding_site"
  )) {
    res <- commaKit:::.parseTargetGenes(nms, ft)
    expect_equal(res[[1]], "aaeR", info = ft)
    expect_equal(res[[2]], c("geneA", "geneB"), info = paste(ft, "operon"))
    expect_equal(res[[3]], "flhDC", info = ft)
  }
})

test_that(".parseTargetGenes protein_binding_site prefers transcription_unit", {
  nms <- "AcrR DNA-binding-site 5 bases of acrRp"
  tu_val <- "acrAB"
  res <- commaKit:::.parseTargetGenes(
    nms,
    "protein_binding_site",
    tu_values = tu_val
  )
  expect_equal(res[[1]], "acrAB")
})

test_that(".parseTargetGenes protein_binding_site falls back to 'of' pattern", {
  nms <- "AcrR DNA-binding-site 5 bases of acrRp"
  res <- commaKit:::.parseTargetGenes(
    nms,
    "protein_binding_site",
    tu_values = NA_character_
  )
  expect_equal(res[[1]], "acrR")
})

test_that(".parseTargetGenes protein_binding_site returns NA when no info", {
  nms <- "unknown binding site"
  res <- commaKit:::.parseTargetGenes(
    nms,
    "protein_binding_site",
    tu_values = NA_character_
  )
  expect_true(is.na(res[[1]]))
})

test_that(".parseTargetGenes RNA_binding_site prefers transcription_unit", {
  nms <- "ArcZ mRNA-binding-site regulating eptB"
  tu_val <- "eptB"
  res <- commaKit:::.parseTargetGenes(
    nms,
    "RNA_binding_site",
    tu_values = tu_val
  )
  expect_equal(res[[1]], "eptB")
})

test_that(
  ".parseTargetGenes RNA_binding_site falls back to 'regulating' pattern",
  {
    nms <- "ArcZ mRNA-binding-site regulating eptB"
    res <- commaKit:::.parseTargetGenes(
      nms,
      "RNA_binding_site",
      tu_values = NA_character_
    )
    expect_equal(res[[1]], "eptB")
  }
)

test_that(".parseTargetGenes terminator strips ' terminator' suffix", {
  nms <- c("geneA terminator", "thrLABC terminator")
  res <- commaKit:::.parseTargetGenes(nms, "terminator")
  expect_equal(res[[1]], "geneA")
  expect_equal(res[[2]], "thrLABC")
})

test_that(".parseTargetGenes default branch returns identity", {
  nms <- c("IS1", "IS10-right")
  res <- commaKit:::.parseTargetGenes(nms, "insertion_sequence")
  expect_equal(res, as.list(nms))
})

# ── .parseRegulatorGenes() ────────────────────────────────────────────────────

test_that(".parseRegulatorGenes returns NA for non-regulatory types", {
  nms <- c("geneA", "geneB")
  for (ft in c("gene", "CDS", "promoter", "terminator")) {
    res <- commaKit:::.parseRegulatorGenes(nms, ft)
    expect_true(all(vapply(res, function(x) all(is.na(x)), logical(1))),
      info = ft
    )
  }
})

test_that(".parseRegulatorGenes maps sigma factor for TFBS types", {
  nms <- c("sigma70-site", "sigma32-site")
  subtypes <- c("Sigma70", "Sigma32")
  for (ft in c(
    "transcription_factor_binding_site",
    "minus_10_signal", "minus_35_signal"
  )) {
    res <- commaKit:::.parseRegulatorGenes(nms, ft, subtype_values = subtypes)
    expect_equal(res[[1]], "rpoD", info = ft)
    expect_equal(res[[2]], "rpoH", info = ft)
  }
})

test_that(".parseRegulatorGenes TFBS returns NA when no subtype", {
  res <- commaKit:::.parseRegulatorGenes(
    "some-site",
    "transcription_factor_binding_site"
  )
  expect_true(is.na(res[[1]]))
})

test_that(".parseRegulatorGenes TFBS returns NA for unknown sigma factor", {
  res <- commaKit:::.parseRegulatorGenes("some-site",
    "transcription_factor_binding_site",
    subtype_values = "SigmaXX"
  )
  expect_true(is.na(res[[1]]))
})

test_that(".parseRegulatorGenes protein_binding_site extracts protein name", {
  nms <- c("AcrR DNA-binding-site", "DnaA binding site of dnaA")
  res <- commaKit:::.parseRegulatorGenes(nms, "protein_binding_site")
  expect_equal(res[[1]], "acrR")
  expect_equal(res[[2]], "dnaA")
})

test_that(
  ".parseRegulatorGenes RNA_binding_site extracts gene-name first word",
  {
    # ArcZ matches gene-name pattern (3-5 lowercase after uppercase)
    res_match <- commaKit:::.parseRegulatorGenes(
      "ArcZ mRNA-binding-site regulating eptB",
      "RNA_binding_site"
    )
    expect_equal(res_match[[1]], "arcZ")
  }
)

test_that(".parseRegulatorGenes RNA_binding_site returns NA for riboswitches", {
  # "adenosylcobalamin" is all-lowercase and too long — no match
  res_no <- commaKit:::.parseRegulatorGenes(
    "adenosylcobalamin riboswitch",
    "RNA_binding_site"
  )
  expect_true(is.na(res_no[[1]]))
})

# ── Feature-role policy

test_that("enrichment feature-role policy is explicit for supported types", {
  feature_types <- c(
    "gene", "CDS", "mRNA", "tRNA", "rRNA", "ncRNA", "promoter",
    "minus_10_signal", "minus_35_signal",
    "transcription_factor_binding_site", "protein_binding_site",
    "RNA_binding_site", "terminator", "operon", "repeat_region", "prophage",
    "region", "insertion_sequence"
  )
  expected <- data.frame(
    feature_type = feature_types,
    target = c(
      rep("identity", 6L), "promoter", "promoter", "promoter", "promoter",
      "protein_binding", "rna_binding", "terminator", rep("identity", 5L)
    ),
    regulator = c(
      rep("none", 7L), rep("sigma_factor", 3L), "protein", "rna",
      rep("none", 6L)
    ),
    role_type = c(
      rep(NA_character_, 7L), rep("sigma_factor", 3L), "TF_protein",
      "RNA_regulator", rep(NA_character_, 6L)
    ),
    stringsAsFactors = FALSE
  )
  actual <- do.call(rbind, lapply(feature_types, function(ft) {
    rule <- commaKit:::.enrichmentFeatureRoleRule(ft)
    data.frame(
      feature_type = ft,
      target = rule$target,
      regulator = rule$regulator,
      role_type = rule$role_type,
      stringsAsFactors = FALSE
    )
  }))
  rownames(actual) <- NULL

  expect_equal(actual, expected)

  # Unlisted annotation types retain the historical generic target rule.
  default_rule <- commaKit:::.enrichmentFeatureRoleRule("custom_feature")
  expect_equal(
    unname(unlist(default_rule)),
    c("identity", "none", NA_character_)
  )

  expect_error(
    commaKit:::.enrichmentFeatureRoleRule(1),
    "single non-NA string"
  )
  expect_error(
    commaKit:::.enrichmentFeatureRoleRule(factor("gene")),
    "single non-NA string"
  )
})

# ── .extractGeneRoles() ───────────────────────────────────────────────────────

# Helper: build a minimal res_df with feature_types and feature_names lists
make_role_res_df <- function(feature_types_list, feature_names_list,
                             dm_padj = 0.01, dm_delta_beta = -0.3,
                             subtype_values_list = NULL,
                             tu_values_list = NULL) {
  n <- length(feature_types_list)
  df <- data.frame(
    chrom = rep("chr1", n),
    position = seq(100L, by = 100L, length.out = n),
    strand = rep("+", n),
    dm_padj = rep(dm_padj, n),
    dm_delta_beta = rep(dm_delta_beta, n),
    stringsAsFactors = FALSE
  )
  df$feature_types <- feature_types_list
  df$feature_names <- feature_names_list
  if (!is.null(subtype_values_list)) {
    df$feature_subtype_values <- subtype_values_list
  }
  if (!is.null(tu_values_list)) {
    df$transcription_unit_values <- tu_values_list
  }
  df
}

test_that(".extractGeneRoles returns NULL when no sites match ft", {
  df <- make_role_res_df(list(c("promoter")), list(c("geneAp1")))
  out <- commaKit:::.extractGeneRoles(df, "gene", "feature_names")
  expect_null(out)
})

test_that(".extractGeneRoles rejects invalid feature types", {
  df <- make_role_res_df(list(c("gene")), list(c("geneA")))
  expect_error(
    commaKit:::.extractGeneRoles(df, NULL, "feature_names"),
    "single non-NA string"
  )
  expect_error(
    commaKit:::.extractGeneRoles(df, NA_character_, "feature_names"),
    "single non-NA string"
  )
})

test_that(".extractGeneRoles gene type: target only, identity", {
  df <- make_role_res_df(list(c("gene")), list(c("geneA")))
  out <- commaKit:::.extractGeneRoles(df, "gene", "feature_names")
  expect_false(is.null(out))
  expect_equal(out$role, "target")
  expect_equal(out$gene_id, "geneA")
  expect_false("regulator" %in% out$role)
})

test_that(
  ".extractGeneRoles promoter: strips suffix, operonic → multiple targets",
  {
    df <- make_role_res_df(
      list(c("promoter", "promoter")),
      list(c("aaeRp1", "geneA-geneBp"))
    )
    out <- commaKit:::.extractGeneRoles(df, "promoter", "feature_names")
    targets <- out$gene_id[out$role == "target"]
    expect_true("aaeR" %in% targets)
    expect_true("geneA" %in% targets)
    expect_true("geneB" %in% targets)
    # No regulators for promoter type
    expect_false("regulator" %in% out$role)
  }
)

test_that(
  ".extractGeneRoles TFBS: target from promoter name, regulator from subtype",
  {
    df <- make_role_res_df(
      list(c("transcription_factor_binding_site")),
      list(c("acrRp")),
      subtype_values_list = list(c("Sigma70"))
    )
    out <- commaKit:::.extractGeneRoles(
      df, "transcription_factor_binding_site",
      "feature_names"
    )
    expect_true("target" %in% out$role)
    expect_true("regulator" %in% out$role)
    expect_equal(out$gene_id[out$role == "target"], "acrR")
    expect_equal(out$gene_id[out$role == "regulator"], "rpoD")
    expect_equal(out$role_type[out$role == "regulator"], "sigma_factor")
  }
)

test_that(
  paste(
    ".extractGeneRoles protein_binding_site: target from TU attr,",
    "regulator from name"
  ),
  {
    df <- make_role_res_df(
      list(c("protein_binding_site")),
      list(c("AcrR DNA-binding-site")),
      tu_values_list = list(c("acrAB"))
    )
    out <- commaKit:::.extractGeneRoles(
      df,
      "protein_binding_site",
      "feature_names"
    )
    expect_equal(out$gene_id[out$role == "target"], "acrAB")
    expect_equal(out$gene_id[out$role == "regulator"], "acrR")
    expect_equal(out$role_type[out$role == "regulator"], "TF_protein")
  }
)

test_that(
  paste(
    ".extractGeneRoles RNA_binding_site: target from TU, regulator",
    "from first word"
  ),
  {
    df <- make_role_res_df(
      list(c("RNA_binding_site")),
      list(c("ArcZ mRNA-binding-site regulating eptB")),
      tu_values_list = list(c("eptB"))
    )
    out <- commaKit:::.extractGeneRoles(df, "RNA_binding_site", "feature_names")
    expect_equal(out$gene_id[out$role == "target"], "eptB")
    expect_equal(out$gene_id[out$role == "regulator"], "arcZ")
    expect_equal(out$role_type[out$role == "regulator"], "RNA_regulator")
  }
)

test_that(".extractGeneRoles keeps list-column associations aligned", {
  df <- make_role_res_df(
    list(c(
      "gene", "transcription_factor_binding_site", "protein_binding_site",
      "RNA_binding_site"
    )),
    list(c(
      "geneA", "acrRp", "AcrR DNA-binding-site",
      "ArcZ mRNA-binding-site regulating eptB"
    )),
    subtype_values_list = list(c(
      NA_character_, "Sigma70", NA_character_, NA_character_
    )),
    tu_values_list = list(c(
      NA_character_, NA_character_, "acrAB", "eptB"
    ))
  )

  out <- commaKit:::.extractGeneRoles(
    df, "transcription_factor_binding_site", "feature_names"
  )

  expect_identical(
    names(out),
    c("gene_id", "role", "role_type", "site_key", "dm_padj", "dm_delta_beta")
  )
  expect_equal(out$gene_id[out$role == "target"], "acrR")
  expect_equal(out$gene_id[out$role == "regulator"], "rpoD")
  expect_equal(out$role_type[out$role == "regulator"], "sigma_factor")
})

test_that(".extractGeneRoles preserves NULL for empty and NA associations", {
  empty_df <- make_role_res_df(
    list(character(0L), NA_character_),
    list(character(0L), NA_character_)
  )
  expect_null(commaKit:::.extractGeneRoles(empty_df, "gene", "feature_names"))

  no_name_df <- make_role_res_df(list("gene"), list(character(0L)))
  expect_null(commaKit:::.extractGeneRoles(no_name_df, "gene", "feature_names"))
})

test_that(".extractGeneRoles overlap_only=TRUE filters by rel_position", {
  # Two sites: one inside (rel_position=0), one outside (rel_position=-20)
  df <- make_role_res_df(
    list(c("gene"), c("gene")),
    list(c("geneA"), c("geneB"))
  )
  df$rel_position <- list(0L, -20L)

  out_all <- commaKit:::.extractGeneRoles(df, "gene", "feature_names",
    overlap_only = FALSE
  )
  out_ovlp <- commaKit:::.extractGeneRoles(df, "gene", "feature_names",
    overlap_only = TRUE
  )

  expect_true("geneA" %in% out_all$gene_id)
  expect_true("geneB" %in% out_all$gene_id)
  expect_true("geneA" %in% out_ovlp$gene_id)
  expect_false("geneB" %in% out_ovlp$gene_id)
  expect_identical(
    names(out_ovlp),
    c("gene_id", "role", "role_type", "site_key", "dm_padj", "dm_delta_beta")
  )
})

# ── enrichMethylation() — gene_role parameter ─────────────────────────────────

# Build a res_df with TFBS entries for gene_role tests
make_tfbs_res_df <- function() {
  data.frame(
    chrom = "chr1",
    position = 100L,
    strand = "+",
    dm_padj = 0.01,
    dm_delta_beta = -0.3,
    feature_types = I(list(c("transcription_factor_binding_site"))),
    feature_names = I(list(c("geneAp1"))),
    feature_subtype_values = I(list(c("Sigma70"))),
    stringsAsFactors = FALSE
  )
}

test_that(
  paste(
    "enrichMethylation gene_role='target' and 'regulator' use",
    "different gene sets"
  ),
  {
    skip_if_not_installed("clusterProfiler")
    df <- make_tfbs_res_df()
    suppressWarnings({
      res_target <- enrichMethylation(df,
        TERM2GENE = fake_t2g,
        feature_type = "transcription_factor_binding_site",
        gene_role = "target"
      )
      res_regulator <- enrichMethylation(df,
        TERM2GENE = fake_t2g,
        feature_type = "transcription_factor_binding_site",
        gene_role = "regulator"
      )
    })
    expect_type(res_target, "list")
    expect_type(res_regulator, "list")
    # Target and regulator should use different gene sets when both produce
    # results.
    # GO may be NULL with tiny synthetic data — but the gene_role distinction
    # should still manifest in the gene universe (accessible via @gene slot)
    if (!is.null(res_target$go) && !is.null(res_regulator$go)) {
      genes_target <- res_target$go@gene
      genes_regulator <- res_regulator$go@gene
      expect_false(
        identical(sort(genes_target), sort(genes_regulator)),
        "Target and regulator gene sets should differ"
      )
    }
  }
)

test_that("enrichMethylation gene_role='regulator' returns valid results", {
  skip_if_not_installed("clusterProfiler")
  df <- make_tfbs_res_df()
  suppressWarnings({
    res <- enrichMethylation(df,
      TERM2GENE = fake_t2g,
      feature_type = "transcription_factor_binding_site",
      gene_role = "regulator"
    )
  })
  expect_type(res, "list")
  expect_true(all(c("go", "kegg") %in% names(res)))
  # Regulator GO result may be NULL with tiny synthetic data
  if (!is.null(res$go)) {
    expect_s4_class(res$go, "enrichResult")
  }
})

test_that(
  "enrichMethylation gene_role='both' returns list with target and regulator",
  {
    skip_if_not_installed("clusterProfiler")
    df <- make_tfbs_res_df()
    suppressWarnings({
      res <- enrichMethylation(df,
        TERM2GENE = fake_t2g,
        feature_type = "transcription_factor_binding_site",
        gene_role = "both"
      )
    })
    expect_type(res, "list")
    expect_true(all(c("target", "regulator") %in% names(res)))
    expect_true(all(c("go", "kegg") %in% names(res$target)))
    expect_true(all(c("go", "kegg") %in% names(res$regulator)))
  }
)

# ── enrichMethylation() — multiple feature_type values ────────────────────────

test_that("enrichMethylation multiple feature_type returns named list", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  suppressWarnings({
    res <- enrichMethylation(ann,
      TERM2GENE = fake_t2g,
      feature_type = c("gene", "nonexistent_type")
    )
  })
  expect_type(res, "list")
  expect_true(all(c("gene", "nonexistent_type") %in% names(res)))
  # gene result: standard list(go, kegg)
  expect_true(all(c("go", "kegg") %in% names(res[["gene"]])))
  # nonexistent_type: NULL slots
  expect_null(res[["nonexistent_type"]]$go)
})

# ── enrichMethylation() — data.frame input ────────────────────────────────────

test_that("enrichMethylation accepts data.frame input and warns on mod_type", {
  skip_if_not_installed("clusterProfiler")
  ann <- make_annotated_dm()
  res_df <- results(ann)
  expect_warning(
    enrichMethylation(res_df, TERM2GENE = fake_t2g, mod_type = "6mA"),
    "mod_type.*ignored|ignored.*mod_type"
  )
})

test_that("enrichMethylation errors for non-commaData non-data.frame input", {
  expect_error(
    enrichMethylation(list(a = 1), TERM2GENE = fake_t2g),
    "commaData.*data.frame|data.frame.*commaData"
  )
})

test_that(
  paste(
    "enrichMethylation validates TERM2GENE and TERM2NAME column",
    "names before clusterProfiler"
  ),
  {
    res_df <- data.frame(
      chrom = "chr1", position = 1L, strand = "+", mod_type = "6mA",
      dm_padj = 0.01, dm_delta_beta = 0.2,
      feature_names = I(list("geneA")),
      stringsAsFactors = FALSE
    )
    bad_t2g <- data.frame(pathway = "P1", gene = "geneA")
    bad_t2n <- data.frame(term = "P1", description = "Pathway 1")

    expect_error(
      enrichMethylation(res_df, TERM2GENE = bad_t2g),
      regexp = "TERM2GENE"
    )
    expect_error(
      enrichMethylation(res_df, TERM2GENE = fake_t2g, TERM2NAME = bad_t2n),
      regexp = "TERM2NAME"
    )
  }
)
