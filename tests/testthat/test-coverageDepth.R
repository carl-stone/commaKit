.make_coverage_depth_object <- function() {
    positions <- c(1L, 5L, 10L, 21L, 30L)
    methyl <- matrix(
        c(0.1, 0.2, 0.3, 0.4, 0.5,
          0.2, 0.3, 0.4, 0.5, 0.6),
        nrow = length(positions),
        dimnames = list(NULL, c("s1", "s2"))
    )
    cov <- matrix(
        c(10L, 20L, 35L, 40L, 60L,
          5L, 15L, 25L, 35L, 45L),
        nrow = length(positions),
        dimnames = list(NULL, c("s1", "s2"))
    )
    site_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr_test", length(positions)),
        ranges = IRanges::IRanges(start = positions, width = 1L),
        strand = rep("+", length(positions)),
        mod_type = factor(rep("6mA", length(positions)),
                          levels = c("4mC", "5mC", "6mA")),
        motif = rep("GATC", length(positions))
    )
    GenomeInfoDb::seqinfo(site_gr) <- GenomeInfoDb::Seqinfo(
        seqnames = "chr_test",
        seqlengths = 40L,
        isCircular = FALSE
    )
    col_data <- S4Vectors::DataFrame(
        sample_name = c("s1", "s2"),
        condition = c("control", "treatment"),
        replicate = c(1L, 1L),
        row.names = c("s1", "s2")
    )
    rse <- SummarizedExperiment::SummarizedExperiment(
        assays = list(methylation = methyl, coverage = cov),
        rowRanges = site_gr,
        colData = col_data
    )
    new("commaData", rse)
}

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

test_that("coverageDepth: method='median' and 'mean' both return valid results with correct row counts", {
    data(comma_example_data)
    result_median <- coverageDepth(comma_example_data, window = 10000L, method = "median")
    result_mean   <- coverageDepth(comma_example_data, window = 10000L, method = "mean")
    # Both methods should return valid data.frames with same structure
    expect_s3_class(result_median, "data.frame")
    expect_s3_class(result_mean, "data.frame")
    expect_true("depth" %in% colnames(result_median))
    expect_true("depth" %in% colnames(result_mean))
    # Both should have the same number of rows (same windows)
    expect_equal(nrow(result_median), nrow(result_mean))
    # Same sample names
    expect_equal(sort(unique(result_median$sample_name)),
                 sort(unique(result_mean$sample_name)))
    # Both should produce finite numeric depths where coverage exists
    med_depths <- result_median$depth[!is.na(result_median$depth)]
    mean_depths <- result_mean$depth[!is.na(result_mean$depth)]
    expect_true(length(med_depths) > 0)
    expect_true(length(mean_depths) > 0)
    expect_true(all(is.finite(med_depths)))
    expect_true(all(is.finite(mean_depths)))
})

test_that("coverageDepth: exact window aggregation keeps empty windows as NA", {
    object <- .make_coverage_depth_object()

    result_mean <- coverageDepth(object, window = 10L, method = "mean")
    result_median <- coverageDepth(object, window = 10L, method = "median")

    expect_equal(nrow(result_mean), 8L)
    expect_equal(
        result_mean[result_mean$sample_name == "s1", "depth"],
        c(mean(c(10, 20, 35)), NA_real_, mean(c(40, 60)), NA_real_)
    )
    expect_equal(
        result_mean[result_mean$sample_name == "s2", "depth"],
        c(mean(c(5, 15, 25)), NA_real_, mean(c(35, 45)), NA_real_)
    )
    expect_equal(
        result_median[result_median$sample_name == "s1", "depth"],
        c(20, NA_real_, 50, NA_real_)
    )
    expect_equal(
        result_median[result_median$sample_name == "s2", "depth"],
        c(15, NA_real_, 40, NA_real_)
    )
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
