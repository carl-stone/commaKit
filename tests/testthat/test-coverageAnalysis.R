test_that("coverageDepth: returns a data.frame", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L)
    expect_s3_class(result, "data.frame")
})

test_that("coverageDepth: has required columns", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L)
    expect_true(all(c("chrom", "window_start", "window_end", "sample_name", "depth") %in%
                        colnames(result)))
})

test_that("coverageDepth: all sample names present in output", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L)
    expect_setequal(unique(result$sample_name), sampleInfo(comma_example_data)$sample_name)
})

test_that("coverageDepth: window_start and window_end are consistent", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L)
    expect_true(all(result$window_end >= result$window_start))
    expect_true(all(result$window_end - result$window_start < 10000L))
})

test_that("coverageDepth: log2_transform adds log2_depth column", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L, log2_transform = TRUE)
    expect_true("log2_depth" %in% colnames(result))
})

test_that("coverageDepth: log2_depth is log2(depth+1)", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L, log2_transform = TRUE)
    non_na <- !is.na(result$depth) & !is.na(result$log2_depth)
    expect_equal(result$log2_depth[non_na], log2(result$depth[non_na] + 1), tolerance = 1e-10)
})

test_that("coverageDepth: method='median' works without error", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L, method = "median")
    expect_s3_class(result, "data.frame")
    expect_true("depth" %in% colnames(result))
})

test_that("coverageDepth: depth values are non-negative where not NA", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L)
    depths <- result$depth[!is.na(result$depth)]
    expect_true(all(depths >= 0))
})

test_that("coverageDepth: window_start begins at 1", {
    data(comma_example_data)
    result <- coverageDepth(comma_example_data, window = 10000L)
    expect_true(1L %in% result$window_start)
})

test_that("coverageDepth: error on non-commaData input", {
    expect_error(coverageDepth(data.frame(x = 1), window = 1000), "'object' must be a commaData")
})

test_that("coverageDepth: error on missing window", {
    data(comma_example_data)
    expect_error(coverageDepth(comma_example_data))
})

test_that("varianceByDepth: returns a data.frame", {
    data(comma_example_data)
    result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
    expect_s3_class(result, "data.frame")
})

test_that("varianceByDepth: has required columns", {
    data(comma_example_data)
    result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
    expect_true(all(c("coverage", "sample_name", "variance", "n_sites") %in% colnames(result)))
})

test_that("varianceByDepth: all sample names present", {
    data(comma_example_data)
    result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
    expect_setequal(unique(result$sample_name), sampleInfo(comma_example_data)$sample_name)
})

test_that("varianceByDepth: coverage levels match coverage_bins", {
    data(comma_example_data)
    bins <- 5:30
    result <- varianceByDepth(comma_example_data, coverage_bins = bins)
    expect_setequal(unique(result$coverage), bins)
})

test_that("varianceByDepth: variance is non-negative or NA", {
    data(comma_example_data)
    result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
    variances <- result$variance[!is.na(result$variance)]
    expect_true(all(variances >= 0))
})

test_that("varianceByDepth: NA variance when fewer than 2 sites at a level", {
    data(comma_example_data)
    result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
    # Where n_sites < 2, variance should be NA
    low_n <- result[!is.na(result$n_sites) & result$n_sites < 2, ]
    if (nrow(low_n) > 0) {
        expect_true(all(is.na(low_n$variance)))
    }
})

test_that("varianceByDepth: mod_type filtering works", {
    data(comma_example_data)
    result_6mA <- varianceByDepth(comma_example_data, coverage_bins = 5:30, mod_type = "6mA")
    result_5mC <- varianceByDepth(comma_example_data, coverage_bins = 5:30, mod_type = "5mC")
    # n_sites sums differ between mod types
    expect_false(
        identical(
            sum(result_6mA$n_sites, na.rm = TRUE),
            sum(result_5mC$n_sites, na.rm = TRUE)
        )
    )
})

test_that("varianceByDepth: NULL coverage_bins uses all unique levels", {
    data(comma_example_data)
    result <- varianceByDepth(comma_example_data, coverage_bins = NULL)
    # Should have multiple coverage levels
    expect_true(length(unique(result$coverage)) > 1)
})

test_that("varianceByDepth: error on non-commaData input", {
    expect_error(varianceByDepth(data.frame(x = 1)), "'object' must be a commaData")
})

test_that("varianceByDepth: error on unknown mod_type", {
    data(comma_example_data)
    expect_error(varianceByDepth(comma_example_data, mod_type = "9mZ"), "not found in object")
})
