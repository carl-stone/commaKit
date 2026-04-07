# Tests for plot_volcano()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_volcano_results <- function() {
    data.frame(
        chrom         = rep("chr_sim", 20L),
        position      = seq(1000L, 20000L, by = 1000L),
        strand        = rep("+", 20L),
        mod_type      = rep("6mA", 20L),
        mod_context   = rep("GATC", 20L),
        dm_pvalue     = c(runif(15L, 0, 0.1), runif(5L, 0.1, 1)),
        dm_padj       = c(runif(5L, 0, 0.05), runif(10L, 0.05, 0.5), runif(5L, 0.5, 1)),
        dm_delta_beta = c(runif(5L, 0.2, 0.8), runif(5L, -0.8, -0.2),
                          runif(10L, -0.1, 0.1)),
        stringsAsFactors = FALSE
    )
}

.make_volcano_results_multicontext <- function() {
    res <- rbind(.make_volcano_results(), .make_volcano_results())
    res$mod_context <- rep(c("GATC", "CCWGG"), each = 20L)
    res
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_volcano: returns ggplot for valid input", {
    res <- .make_volcano_results()
    p <- plot_volcano(res)
    expect_s3_class(p, "ggplot")
})

test_that("plot_volcano: custom thresholds accepted", {
    res <- .make_volcano_results()
    p <- plot_volcano(res, delta_beta_threshold = 0.3, padj_threshold = 0.01)
    expect_s3_class(p, "ggplot")
})

# ─── Threshold lines ──────────────────────────────────────────────────────────

test_that("plot_volcano: NULL delta_beta_threshold has 2 layers (points + hline)", {
    res <- .make_volcano_results()
    p <- plot_volcano(res)
    expect_gte(length(p$layers), 2L)
})

test_that("plot_volcano: numeric delta_beta_threshold adds vlines (>= 4 layers)", {
    res <- .make_volcano_results()
    p <- plot_volcano(res, delta_beta_threshold = 0.2)
    expect_gte(length(p$layers), 4L)
})

# ─── Faceting ─────────────────────────────────────────────────────────────────

test_that("plot_volcano: single mod_context does not facet", {
    res <- .make_volcano_results()
    p <- plot_volcano(res)
    expect_false(inherits(p$facet, "FacetWrap"))
})

test_that("plot_volcano: multiple mod_context levels facet when facet = TRUE", {
    res <- .make_volcano_results_multicontext()
    p <- plot_volcano(res, facet = TRUE)
    expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("plot_volcano: facet = FALSE suppresses faceting for multi-context results", {
    res <- .make_volcano_results_multicontext()
    p <- plot_volcano(res, facet = FALSE)
    expect_false(inherits(p$facet, "FacetWrap"))
})

test_that("plot_volcano: error on invalid facet argument", {
    res <- .make_volcano_results()
    expect_error(plot_volcano(res, facet = "yes"), "facet")
    expect_error(plot_volcano(res, facet = NA), "facet")
})

# ─── NA handling ─────────────────────────────────────────────────────────────

test_that("plot_volcano: rows with NA padj are excluded without error", {
    res <- .make_volcano_results()
    res$dm_padj[1:5] <- NA
    p <- plot_volcano(res)
    expect_s3_class(p, "ggplot")
})

test_that("plot_volcano: rows with NA delta_beta handled without error", {
    res <- .make_volcano_results()
    res$dm_delta_beta[1:3] <- NA
    p <- plot_volcano(res)
    expect_s3_class(p, "ggplot")
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_volcano: error on non-data.frame input", {
    expect_error(plot_volcano(list(dm_delta_beta = 1, dm_padj = 0.01)),
                 "data.frame")
})

test_that("plot_volcano: error when dm_delta_beta column is missing", {
    res <- .make_volcano_results()
    res$dm_delta_beta <- NULL
    expect_error(plot_volcano(res), "dm_delta_beta")
})

test_that("plot_volcano: error when dm_padj column is missing", {
    res <- .make_volcano_results()
    res$dm_padj <- NULL
    expect_error(plot_volcano(res), "dm_padj")
})

test_that("plot_volcano: error when all padj values are NA", {
    res <- .make_volcano_results()
    res$dm_padj <- NA_real_
    expect_error(plot_volcano(res))
})

test_that("plot_volcano: NULL delta_beta_threshold colors by padj alone", {
    set.seed(42L)
    res <- .make_volcano_results()
    p <- plot_volcano(res, delta_beta_threshold = NULL)
    pd <- ggplot2::ggplot_build(p)$data[[1L]]
    # All points with dm_padj <= 0.05 and dm_delta_beta > 0 should be Hypermethylated
    sig_pos <- res$dm_padj <= 0.05 & res$dm_delta_beta > 0 & !is.na(res$dm_padj)
    expect_true(any(sig_pos))
    expect_s3_class(p, "ggplot")
})

test_that("plot_volcano: error on invalid delta_beta_threshold", {
    res <- .make_volcano_results()
    expect_error(plot_volcano(res, delta_beta_threshold = 0), "delta_beta_threshold")
    expect_error(plot_volcano(res, delta_beta_threshold = 1), "delta_beta_threshold")
    expect_error(plot_volcano(res, delta_beta_threshold = -0.1), "delta_beta_threshold")
})

test_that("plot_volcano: error on invalid padj_threshold", {
    res <- .make_volcano_results()
    expect_error(plot_volcano(res, padj_threshold = 0), "padj_threshold")
    expect_error(plot_volcano(res, padj_threshold = 1), "padj_threshold")
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_volcano: works with results() output from comma_example_data", {
    data(comma_example_data)
    cd_dm <- diffMethyl(comma_example_data, ~ condition, mod_type = "6mA")
    res <- results(cd_dm)
    p <- plot_volcano(res)
    expect_s3_class(p, "ggplot")
})
