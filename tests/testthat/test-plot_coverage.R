# Tests for plot_coverage()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_cov_data <- function() {
    n_sites   <- 10L
    positions <- seq(1000L, 10000L, by = 1000L)
    site_keys <- paste0("chr_sim:", positions, ":+:6mA:GATC")
    set.seed(7L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(site_keys, c("samp1", "samp2"))
    )
    depths <- matrix(
        sample(5L:100L, n_sites * 2L, replace = TRUE),
        nrow = n_sites, ncol = 2L,
        dimnames = dimnames(betas)
    )
    rd <- S4Vectors::DataFrame(
        chrom       = rep("chr_sim", n_sites),
        position    = positions,
        strand      = rep("+", n_sites),
        mod_type    = rep("6mA", n_sites),
        motif       = rep("GATC", n_sites),
        mod_context = rep("6mA_GATC", n_sites),
        row.names   = site_keys
    )
    cd <- S4Vectors::DataFrame(
        sample_name = c("samp1", "samp2"),
        condition   = c("ctrl", "treat"),
        replicate   = 1:2,
        row.names   = c("samp1", "samp2")
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = betas, coverage = depths),
        rowData = rd,
        colData = cd
    )
    new("commaData", se,
        genomeInfo = c(chr_sim = 100000L),
        annotation = GenomicRanges::GRanges(),
        motifSites = GenomicRanges::GRanges())
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_coverage: returns ggplot for valid input", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj)
    expect_s3_class(p, "ggplot")
})

test_that("plot_coverage: per_sample = FALSE returns ggplot", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj, per_sample = FALSE)
    expect_s3_class(p, "ggplot")
})

test_that("plot_coverage: mod_type filter accepted", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj, mod_type = "6mA")
    expect_s3_class(p, "ggplot")
})

# ─── Faceting ─────────────────────────────────────────────────────────────────

test_that("plot_coverage: per_sample = TRUE produces faceted plot", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj, per_sample = TRUE)
    expect_false(inherits(p$facet, "FacetNull"))
})

test_that("plot_coverage: per_sample = FALSE produces unfaceted plot", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj, per_sample = FALSE)
    expect_true(inherits(p$facet, "FacetNull"))
})

# ─── Axis labels ─────────────────────────────────────────────────────────────

test_that("plot_coverage: x-axis label mentions coverage", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj)
    expect_true(grepl("[Cc]overage", p$labels$x))
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_coverage: error on non-commaData input", {
    expect_error(plot_coverage(data.frame(x = 1)), "commaData")
})

test_that("plot_coverage: error on invalid mod_type", {
    obj <- .make_cov_data()
    expect_error(plot_coverage(obj, mod_type = "4mC"), "not found")
})

test_that("plot_coverage: error on invalid per_sample", {
    obj <- .make_cov_data()
    expect_error(plot_coverage(obj, per_sample = "yes"), "per_sample")
})

# ─── Single sample ────────────────────────────────────────────────────────────

test_that("plot_coverage: works with a single-sample object", {
    obj <- .make_cov_data()
    obj_1samp <- obj[, 1L, drop = FALSE]
    p <- plot_coverage(obj_1samp)
    expect_s3_class(p, "ggplot")
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_coverage: works with comma_example_data", {
    data(comma_example_data)
    p <- plot_coverage(comma_example_data)
    expect_s3_class(p, "ggplot")
})
