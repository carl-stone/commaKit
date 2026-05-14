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
    # Use a wide bin range; with 300 sites across 5:30 coverage levels some bins
    # will have 0 or 1 sites, guaranteeing low_n rows exist.
    result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
    low_n <- result[!is.na(result$n_sites) & result$n_sites < 2L, ]
    expect_true(nrow(low_n) > 0L)           # condition must actually be tested
    expect_true(all(is.na(low_n$variance))) # variance must be NA for those rows
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
