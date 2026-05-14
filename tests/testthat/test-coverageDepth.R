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
