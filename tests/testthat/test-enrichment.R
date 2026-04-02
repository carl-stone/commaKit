## tests/testthat/test-enrichment.R
## Tests for enrichMethylation(), .siteToGeneMap(), .computeGeneScores()
## All tests use custom TERM2GENE to avoid requiring clusterProfiler, OrgDb,
## KEGG access, or internet connectivity.

library(testthat)
library(comma)

# ── Shared fixtures ───────────────────────────────────────────────────────────

# Minimal commaData with diffMethyl results and feature_names annotation
make_annotated_dm <- function() {
    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    ann <- annotateSites(dm, annotation(comma_example_data), type = "overlap")
    ann
}

# Minimal TERM2GENE mapping aligned with example data gene names (geneA-E)
fake_t2g <- data.frame(
    term = c("PATH:01", "PATH:01", "PATH:02", "PATH:02", "PATH:03"),
    gene = c("geneA",   "geneB",   "geneC",   "geneD",   "geneE"),
    stringsAsFactors = FALSE
)

fake_t2n <- data.frame(
    term = c("PATH:01", "PATH:02", "PATH:03"),
    name = c("Pathway one", "Pathway two", "Pathway three"),
    stringsAsFactors = FALSE
)

# ── .siteToGeneMap() ─────────────────────────────────────────────────────────

test_that(".siteToGeneMap returns data.frame with expected columns", {
    res_df <- data.frame(
        chrom         = c("chr1", "chr1"),
        position      = c(100L, 200L),
        strand        = c("+", "+"),
        dm_padj       = c(0.01, 0.5),
        dm_delta_beta = c(-0.3, 0.05),
        feature_names = I(list(c("geneA", "geneB"), character(0))),
        stringsAsFactors = FALSE
    )
    sg <- comma:::.siteToGeneMap(res_df, "feature_names")

    expect_s3_class(sg, "data.frame")
    expect_true(all(c("gene_id", "site_key", "dm_padj", "dm_delta_beta") %in% colnames(sg)))
})

test_that(".siteToGeneMap correctly explodes multi-gene sites", {
    res_df <- data.frame(
        chrom         = "chr1",
        position      = 100L,
        strand        = "+",
        dm_padj       = 0.01,
        dm_delta_beta = -0.3,
        feature_names = I(list(c("geneA", "geneB"))),
        stringsAsFactors = FALSE
    )
    sg <- comma:::.siteToGeneMap(res_df, "feature_names")

    expect_equal(nrow(sg), 2L)
    expect_setequal(sg$gene_id, c("geneA", "geneB"))
    # Both rows inherit the same site statistics
    expect_equal(sg$dm_padj, c(0.01, 0.01))
})

test_that(".siteToGeneMap excludes intergenic sites silently", {
    res_df <- data.frame(
        chrom         = c("chr1", "chr1"),
        position      = c(100L, 200L),
        strand        = "+",
        dm_padj       = c(0.01, 0.5),
        dm_delta_beta = c(-0.3, 0.05),
        feature_names = I(list(c("geneA"), character(0))),  # second site intergenic
        stringsAsFactors = FALSE
    )
    sg <- comma:::.siteToGeneMap(res_df, "feature_names")

    expect_equal(nrow(sg), 1L)
    expect_equal(sg$gene_id, "geneA")
})

test_that(".siteToGeneMap returns empty data.frame when all sites intergenic", {
    res_df <- data.frame(
        chrom         = "chr1",
        position      = 100L,
        strand        = "+",
        dm_padj       = 0.01,
        dm_delta_beta = -0.3,
        feature_names = I(list(character(0))),
        stringsAsFactors = FALSE
    )
    expect_warning(
        sg <- comma:::.siteToGeneMap(res_df, "feature_names"),
        "No sites have gene annotations"
    )
    expect_equal(nrow(sg), 0L)
})

test_that(".siteToGeneMap errors when gene_col is missing", {
    res_df <- data.frame(chrom="chr1", position=100L, strand="+",
                         dm_padj=0.01, dm_delta_beta=-0.3,
                         stringsAsFactors=FALSE)
    expect_error(
        comma:::.siteToGeneMap(res_df, "feature_names"),
        "not found in results"
    )
})

test_that(".siteToGeneMap constructs correct site_key", {
    res_df <- data.frame(
        chrom         = "chrX",
        position      = 555L,
        strand        = "-",
        dm_padj       = 0.001,
        dm_delta_beta = 0.4,
        feature_names = I(list("geneZ")),
        stringsAsFactors = FALSE
    )
    sg <- comma:::.siteToGeneMap(res_df, "feature_names")
    expect_equal(sg$site_key, "chrX:555:-")
})

# ── .computeGeneScores() ─────────────────────────────────────────────────────

make_sg <- function() {
    data.frame(
        gene_id       = c("geneA", "geneA", "geneB", "geneC"),
        site_key      = c("chr1:100:+", "chr1:200:+", "chr1:300:+", "chr1:400:+"),
        dm_padj       = c(0.001, 0.01, 0.5, NA),
        dm_delta_beta = c(-0.5, -0.2, 0.05, 0.3),
        stringsAsFactors = FALSE
    )
}

test_that(".computeGeneScores returns named sorted vector", {
    scores <- comma:::.computeGeneScores(make_sg(), "combined", "max")

    expect_true(is.numeric(scores))
    expect_true(!is.null(names(scores)))
    # Must be sorted decreasing
    expect_equal(scores, sort(scores, decreasing = TRUE))
})

test_that(".computeGeneScores 'combined' metric has correct sign", {
    sg <- data.frame(
        gene_id       = "geneA",
        site_key      = "chr1:100:+",
        dm_padj       = 0.01,
        dm_delta_beta = 0.5,   # positive delta_beta → positive score
        stringsAsFactors = FALSE
    )
    scores <- comma:::.computeGeneScores(sg, "combined", "max")
    expect_true(scores[["geneA"]] > 0)

    sg$dm_delta_beta <- -0.5   # negative delta_beta → negative score
    scores2 <- comma:::.computeGeneScores(sg, "combined", "max")
    expect_true(scores2[["geneA"]] < 0)
})

test_that(".computeGeneScores excludes sites with NA padj", {
    # geneC has NA padj; should not appear in scores
    scores <- comma:::.computeGeneScores(make_sg(), "combined", "max")
    expect_false("geneC" %in% names(scores))
})

test_that(".computeGeneScores 'max' aggregation picks largest absolute score", {
    # geneA has two sites: padj=0.001 (larger -log10) and padj=0.01
    scores <- comma:::.computeGeneScores(make_sg(), "padj", "max")
    # -log10(0.001) = 3 > -log10(0.01) = 2
    expect_equal(scores[["geneA"]], -log10(0.001), tolerance = 1e-6)
})

test_that(".computeGeneScores 'mean' aggregation averages across sites", {
    sg <- data.frame(
        gene_id       = c("geneA", "geneA"),
        site_key      = c("a", "b"),
        dm_padj       = c(0.01, 0.01),
        dm_delta_beta = c(0.4, 0.6),
        stringsAsFactors = FALSE
    )
    scores <- comma:::.computeGeneScores(sg, "delta_beta", "mean")
    expect_equal(scores[["geneA"]], 0.5, tolerance = 1e-6)
})

test_that(".computeGeneScores returns empty vector when all NA", {
    sg <- data.frame(
        gene_id=c("g1"), site_key="a", dm_padj=NA_real_, dm_delta_beta=NA_real_,
        stringsAsFactors=FALSE
    )
    scores <- comma:::.computeGeneScores(sg)
    expect_length(scores, 0L)
})

test_that(".computeGeneScores errors on unknown score_metric", {
    expect_error(
        comma:::.computeGeneScores(make_sg(), "bogus"),
        "score_metric"
    )
})

# ── enrichMethylation() — error conditions ────────────────────────────────────

test_that("enrichMethylation errors when clusterProfiler is absent", {
    skip_if(requireNamespace("clusterProfiler", quietly=TRUE),
            "clusterProfiler is installed; cannot test missing-package error")
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

test_that("enrichMethylation errors when annotateSites not run", {
    skip_if_not_installed("clusterProfiler")
    data(comma_example_data)
    dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    # dm has no feature_names column
    expect_error(
        enrichMethylation(dm, TERM2GENE = fake_t2g),
        "not found in results|annotateSites"
    )
})

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

test_that("enrichMethylation ORA with TERM2NAME works without error", {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    expect_no_error(
        enrichMethylation(ann, method = "ora",
                          TERM2GENE = fake_t2g, TERM2NAME = fake_t2n)
    )
})

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
    res <- enrichMethylation(ann, method = "gsea", TERM2GENE = fake_t2g,
                             minGSSize = 1L)

    if (!is.null(res$go)) {
        expect_s4_class(res$go, "gseaResult")
    }
})

# ── enrichMethylation() — both ORA and GSEA ──────────────────────────────────

test_that("enrichMethylation both methods returns nested list structure", {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    res <- enrichMethylation(ann, method = c("ora", "gsea"),
                             TERM2GENE = fake_t2g, minGSSize = 1L)

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
        res <- enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g,
                                 padj_threshold = 0, delta_beta_threshold = 2),
        "No significantly differentially methylated genes"
    )
    expect_null(res$go)
    expect_null(res$kegg)
})

# ── enrichMethylation() — mod_type / mod_context filters ─────────────────────

test_that("enrichMethylation mod_type filter passes through to results()", {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    # Should work without error (filters to 6mA sites only)
    expect_no_error(
        enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g,
                          mod_type = "6mA")
    )
})

test_that("enrichMethylation mod_type filter errors on unknown type", {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    expect_error(
        enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g,
                          mod_type = "99mX"),
        "mod_type"
    )
})

# ── .siteToGeneMap() — NA gene_id handling ───────────────────────────────────

test_that(".siteToGeneMap drops NA gene IDs silently", {
    res_df <- data.frame(
        chrom         = c("chr1", "chr1"),
        position      = c(100L, 200L),
        strand        = "+",
        dm_padj       = c(0.01, 0.02),
        dm_delta_beta = c(-0.3, 0.2),
        feature_names = I(list(c(NA_character_, "geneA"), c("geneB"))),
        stringsAsFactors = FALSE
    )
    sg <- comma:::.siteToGeneMap(res_df, "feature_names")

    expect_false(any(is.na(sg$gene_id)))
    expect_true("geneA" %in% sg$gene_id)
    expect_true("geneB" %in% sg$gene_id)
})

# ── .computeGeneScores() — NA gene_id robustness ────────────────────────────

test_that(".computeGeneScores 'max' does not crash when gene_ids contain NA", {
    sg <- data.frame(
        gene_id       = c(NA_character_, "geneA", "geneA"),
        site_key      = c("chr1:100:+", "chr1:200:+", "chr1:300:+"),
        dm_padj       = c(0.01, 0.01, 0.05),
        dm_delta_beta = c(0.4, -0.5, -0.2),
        stringsAsFactors = FALSE
    )
    expect_no_error(
        scores <- comma:::.computeGeneScores(sg, "combined", "max")
    )
    expect_false(is.null(names(scores)))
    expect_false(any(is.na(names(scores))))
})

test_that(".computeGeneScores 'mean' does not return empty when gene_ids contain NA", {
    sg <- data.frame(
        gene_id       = c(NA_character_, "geneA"),
        site_key      = c("chr1:100:+", "chr1:200:+"),
        dm_padj       = c(0.01, 0.01),
        dm_delta_beta = c(0.4, 0.5),
        stringsAsFactors = FALSE
    )
    scores <- comma:::.computeGeneScores(sg, "delta_beta", "mean")
    expect_true(length(scores) > 0L)
    expect_true("geneA" %in% names(scores))
})

# ── enrichMethylation() — feature_type argument ───────────────────────────────

test_that("enrichMethylation feature_type = 'gene' runs without error", {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    expect_no_error(
        enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g,
                          feature_type = "gene")
    )
})

test_that("enrichMethylation feature_type = NULL includes all features", {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    expect_no_error(
        enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g,
                          feature_type = NULL)
    )
})

test_that("enrichMethylation warns and returns NULL for unmatched feature_type", {
    skip_if_not_installed("clusterProfiler")
    ann <- make_annotated_dm()
    expect_warning(
        res <- enrichMethylation(ann, method = "ora", TERM2GENE = fake_t2g,
                                 feature_type = "nonexistent_type"),
        "No sites with feature_type"
    )
    expect_null(res$go)
    expect_null(res$kegg)
})
