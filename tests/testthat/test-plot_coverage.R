# Tests for plot_coverage()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_cov_data <- function() {
    n_sites   <- 10L
    positions <- seq(1000L, 10000L, by = 1000L)
    set.seed(7L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(NULL, c("samp1", "samp2"))
    )
    depths <- matrix(
        sample(5L:100L, n_sites * 2L, replace = TRUE),
        nrow = n_sites, ncol = 2L,
        dimnames = dimnames(betas)
    )
    sample_info <- data.frame(
        sample_name = c("samp1", "samp2"),
        condition   = c("ctrl", "treat"),
        replicate   = 1:2,
        stringsAsFactors = FALSE
    )
    .make_commaData_fixture(betas, depths, sample_info, positions)
}

# ─── Data mapping ─────────────────────────────────────────────────────────────

test_that("plot_coverage: maps exact coverage depths to x aesthetic with correct sample names", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj)
    expect_s3_class(p, "ggplot")
    # p$data should have depth, sample_name, condition columns
    expect_true("depth" %in% colnames(p$data))
    expect_true("sample_name" %in% colnames(p$data))
    # 10 sites * 2 samples = 20 rows (no NAs in this fixture)
    expect_equal(nrow(p$data), 20L)
    # Sample names should match the object
    expect_equal(sort(unique(p$data$sample_name)), c("samp1", "samp2"))
    # Depths should match the coverage assay values exactly
    cov_mat <- siteCoverage(obj)
    expect_equal(sort(p$data$depth), sort(as.vector(cov_mat)))
})

test_that("plot_coverage: mod_type filter to 6mA produces identical data (all sites are 6mA)", {
    obj <- .make_cov_data()
    p_all <- plot_coverage(obj)
    p_filt <- plot_coverage(obj, mod_type = "6mA")
    # All sites are 6mA in this fixture, so data should be identical
    expect_equal(nrow(p_filt$data), nrow(p_all$data))
    expect_equal(p_filt$data$depth, p_all$data$depth)
})

# ─── Faceting ─────────────────────────────────────────────────────────────────

test_that("plot_coverage: per_sample = TRUE facets by sample_name", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj, per_sample = TRUE)
    expect_false(inherits(p$facet, "FacetNull"))
    # FacetWrap: verify the facet variable is sample_name
    # Extract from the facet formula without importing rlang
    facet_formula <- p$facet$vars()
    expect_true("sample_name" %in% facet_formula)
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

# ─── Median vline ─────────────────────────────────────────────────────────────

test_that("plot_coverage: per_sample vline positions match computed medians per sample", {
    obj <- .make_cov_data()
    p <- plot_coverage(obj, per_sample = TRUE)
    # Find the vline layer
    layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    vline_idx <- which(layer_classes == "GeomVline")
    expect_gte(length(vline_idx), 1L)
    # The vline layer data should have median_depth per sample
    vline_data <- p$layers[[vline_idx[1]]]$data
    expect_s3_class(vline_data, "data.frame")
    expect_true("median_depth" %in% colnames(vline_data))
    # Verify exact median values match stats::median of coverage per sample
    cov_mat <- siteCoverage(obj)
    for (samp in vline_data$sample_name) {
        expected_med <- stats::median(cov_mat[, samp])
        actual_med <- vline_data$median_depth[vline_data$sample_name == samp]
        expect_equal(actual_med, expected_med)
    }
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

test_that("plot_coverage: single-sample object plots only that sample's depths", {
    obj <- .make_cov_data()
    obj_1samp <- obj[, 1L, drop = FALSE]
    p <- plot_coverage(obj_1samp)
    expect_s3_class(p, "ggplot")
    # Only one sample in the data
    expect_equal(unique(p$data$sample_name), "samp1")
    expect_equal(nrow(p$data), 10L)  # 10 sites, 1 sample
    # Depths should match the single column exactly
    expect_equal(sort(p$data$depth), sort(as.vector(siteCoverage(obj_1samp))))
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_coverage: comma_example_data has correct row count", {
    data(comma_example_data)
    p <- plot_coverage(comma_example_data)
    expect_s3_class(p, "ggplot")
    # 588 sites * 6 samples = 3528 rows
    expect_equal(nrow(p$data), 588L * 6L)
})

test_that("plot_coverage: mod_type vector filters to matching sites exactly", {
    data(comma_example_data)
    p_6ma <- plot_coverage(comma_example_data, mod_type = "6mA")
    p_both <- plot_coverage(comma_example_data, mod_type = c("6mA", "5mC"))
    # 6mA only: 393 sites * 6 samples
    expect_equal(nrow(p_6ma$data), 393L * 6L)
    # Both: 588 sites * 6 samples
    expect_equal(nrow(p_both$data), 588L * 6L)
})

test_that("plot_coverage: mod_type vector with invalid value gives error", {
    data(comma_example_data)
    expect_error(
        plot_coverage(comma_example_data, mod_type = c("6mA", "invalid")),
        "not found in object"
    )
})
