test_that("slidingWindow: returns a data.frame", {
    data(comma_example_data)
    result <- slidingWindow(comma_example_data, window = 5000L)
    expect_s3_class(result, "data.frame")
})

test_that("slidingWindow: output has required columns", {
    data(comma_example_data)
    result <- slidingWindow(comma_example_data, window = 5000L)
    expect_true(all(c("chrom", "position", "sample_name", "window_median") %in% colnames(result)))
})

test_that("slidingWindow: stat='mean' produces window_mean column", {
    data(comma_example_data)
    result <- slidingWindow(comma_example_data, window = 5000L, stat = "mean")
    expect_true("window_mean" %in% colnames(result))
    expect_false("window_median" %in% colnames(result))
})

test_that("slidingWindow: number of rows equals genome_size * n_samples", {
    data(comma_example_data)
    gi <- genome(comma_example_data)
    n_pos <- sum(gi)
    n_samp <- ncol(comma_example_data)
    result <- slidingWindow(comma_example_data, window = 5000L)
    expect_equal(nrow(result), n_pos * n_samp)
})

test_that("slidingWindow: positions span 1 to chromosome size", {
    data(comma_example_data)
    gi <- genome(comma_example_data)
    chr_size <- gi[1]
    result <- slidingWindow(comma_example_data, window = 5000L)
    chr_result <- result[result$chrom == names(gi)[1], ]
    expect_equal(min(chr_result$position[chr_result$sample_name == chr_result$sample_name[1]]), 1L)
    expect_equal(max(chr_result$position[chr_result$sample_name == chr_result$sample_name[1]]),
                 as.integer(chr_size))
})

test_that("slidingWindow: all sample names present in output", {
    data(comma_example_data)
    result <- slidingWindow(comma_example_data, window = 5000L)
    expect_setequal(unique(result$sample_name), sampleInfo(comma_example_data)$sample_name)
})

test_that("slidingWindow: median and mean give different results", {
    data(comma_example_data)
    r_med  <- slidingWindow(comma_example_data, window = 5000L, stat = "median")
    r_mean <- slidingWindow(comma_example_data, window = 5000L, stat = "mean")
    # They can differ; compare non-NA values
    v_med  <- r_med$window_median[!is.na(r_med$window_median)]
    v_mean <- r_mean$window_mean[!is.na(r_mean$window_mean)]
    expect_false(isTRUE(all.equal(v_med, v_mean, check.attributes = FALSE)))
})

test_that("slidingWindow: mod_type filtering works", {
    data(comma_example_data)
    result_6mA <- slidingWindow(comma_example_data, window = 5000L, mod_type = "6mA")
    result_5mC <- slidingWindow(comma_example_data, window = 5000L, mod_type = "5mC")
    # Both still produce full-genome output (genome size × n_samples)
    gi <- genome(comma_example_data)
    n_samp <- ncol(comma_example_data)
    expect_equal(nrow(result_6mA), sum(gi) * n_samp)
    expect_equal(nrow(result_5mC), sum(gi) * n_samp)
    # The smoothed values should differ because different sites are included
    v6 <- result_6mA$window_median[!is.na(result_6mA$window_median)]
    v5 <- result_5mC$window_median[!is.na(result_5mC$window_median)]
    expect_false(isTRUE(all.equal(v6, v5, check.attributes = FALSE)))
})

test_that("slidingWindow: values are in [0,1] range (ignoring NA)", {
    data(comma_example_data)
    result <- slidingWindow(comma_example_data, window = 5000L)
    vals <- result$window_median[!is.na(result$window_median)]
    expect_true(all(vals >= 0 & vals <= 1))
})

test_that("slidingWindow: circular=FALSE works without error", {
    data(comma_example_data)
    result <- slidingWindow(comma_example_data, window = 5000L, circular = FALSE)
    expect_s3_class(result, "data.frame")
    expect_true("window_median" %in% colnames(result))
})

test_that("slidingWindow: circular=TRUE and FALSE give different edge results", {
    data(comma_example_data)
    r_circ   <- slidingWindow(comma_example_data, window = 5000L, circular = TRUE)
    r_linear <- slidingWindow(comma_example_data, window = 5000L, circular = FALSE)
    gi   <- genome(comma_example_data)
    chr  <- names(gi)[1]
    samp <- sampleInfo(comma_example_data)$sample_name[1]
    # Compare all positions for one sample on the chromosome.
    # With 300 sites across a 100kb genome and a 5000bp window, circular wrapping
    # at the chromosome boundary must produce at least one position with a different
    # smoothed value.
    v_circ   <- r_circ$window_median[r_circ$chrom == chr & r_circ$sample_name == samp]
    v_linear <- r_linear$window_median[r_linear$chrom == chr & r_linear$sample_name == samp]
    expect_false(isTRUE(all.equal(v_circ, v_linear)))
})

test_that("slidingWindow: error on non-commaData input", {
    expect_error(slidingWindow(data.frame(x = 1), window = 100), "'object' must be a commaData")
})

test_that("slidingWindow: error on missing window argument", {
    data(comma_example_data)
    expect_error(slidingWindow(comma_example_data))
})

test_that("slidingWindow: error when genome is NULL", {
    data(comma_example_data)
    obj_no_genome <- new("commaData",
        as(comma_example_data, "SummarizedExperiment"),
        genomeInfo = NULL,
        annotation = comma_example_data@annotation,
        motifSites = comma_example_data@motifSites
    )
    expect_error(slidingWindow(obj_no_genome, window = 1000L), "genome\\(object\\)")
})

test_that("slidingWindow: error when window exceeds chromosome size", {
    data(comma_example_data)
    gi <- genome(comma_example_data)   # chr_sim = 100000 bp
    expect_error(
        slidingWindow(comma_example_data, window = 200000L),
        "exceeds the smallest chromosome size"
    )
})

test_that("slidingWindow: error on invalid mod_type", {
    data(comma_example_data)
    expect_error(
        slidingWindow(comma_example_data, window = 1000L, mod_type = "invalid_type"),
        "No sites remain"
    )
})

test_that("slidingWindow: NaN values converted to NA", {
    data(comma_example_data)
    result <- slidingWindow(comma_example_data, window = 5000L)
    expect_false(any(is.nan(result$window_median)))
})

test_that("slidingWindow: known smoothed value for simple input", {
    # Single chromosome, 3 sites at known positions, window=3
    # We verify the smoothed value at a position adjacent to a site
    gi <- c(chr_test = 20L)
    # Create minimal commaData
    sites_gr <- GenomicRanges::GRanges(
        seqnames = "chr_test",
        ranges   = IRanges::IRanges(start = c(5L, 10L, 15L), width = 1L),
        strand   = "+"
    )
    methyl_mat <- matrix(c(0.2, 0.8, 0.5), nrow = 3, ncol = 1,
                         dimnames = list(NULL, "samp1"))
    cov_mat    <- matrix(10L, nrow = 3, ncol = 1,
                         dimnames = list(NULL, "samp1"))
    rd <- S4Vectors::DataFrame(
        chrom    = "chr_test",
        position = c(5L, 10L, 15L),
        strand   = "+",
        mod_type = "6mA"
    )
    cd <- S4Vectors::DataFrame(
        sample_name = "samp1",
        condition   = "ctrl",
        replicate   = 1L
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = methyl_mat, coverage = cov_mat),
        rowData = rd,
        colData = cd
    )
    obj <- new("commaData", se, genomeInfo = gi,
               annotation = GenomicRanges::GRanges(),
               motifSites = GenomicRanges::GRanges())

    # window=5 centered on position 10 covers positions 8:12.
    # Only the site at position 10 (beta=0.8) falls in [8:12]; sites at 5 and 15
    # are outside. So median({0.8}) == 0.8 exactly.
    result  <- slidingWindow(obj, window = 5L, circular = FALSE)
    val_p10 <- result$window_median[result$position == 10L]
    expect_false(is.na(val_p10))
    expect_equal(val_p10, 0.8, tolerance = 1e-9)
})
