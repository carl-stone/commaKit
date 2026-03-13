test_that("methylomeSummary: returns a data.frame", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    expect_s3_class(result, "data.frame")
})

test_that("methylomeSummary: one row per sample", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    expect_equal(nrow(result), ncol(comma_example_data))
})

test_that("methylomeSummary: has required columns", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    expected_cols <- c("sample_name", "condition", "mod_type", "n_sites",
                       "n_covered", "mean_beta", "median_beta", "sd_beta",
                       "frac_methylated", "mean_coverage", "median_coverage")
    expect_true(all(expected_cols %in% colnames(result)))
})

test_that("methylomeSummary: sample names match object", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    expect_setequal(result$sample_name, sampleInfo(comma_example_data)$sample_name)
})

test_that("methylomeSummary: condition column populated from sampleInfo", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    si <- sampleInfo(comma_example_data)
    for (samp in result$sample_name) {
        expected_cond <- si$condition[si$sample_name == samp]
        actual_cond   <- result$condition[result$sample_name == samp]
        expect_equal(actual_cond, expected_cond)
    }
})

test_that("methylomeSummary: mod_type='all' when not filtered", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    expect_true(all(result$mod_type == "all"))
})

test_that("methylomeSummary: mod_type filtering works", {
    data(comma_example_data)
    result_6mA <- methylomeSummary(comma_example_data, mod_type = "6mA")
    result_5mC <- methylomeSummary(comma_example_data, mod_type = "5mC")
    expect_true(all(result_6mA$mod_type == "6mA"))
    expect_true(all(result_5mC$mod_type == "5mC"))
    # n_sites differs between mod types
    expect_false(result_6mA$n_sites[1] == result_5mC$n_sites[1])
})

test_that("methylomeSummary: n_sites correct after mod_type filter", {
    data(comma_example_data)
    n_6mA <- sum(siteInfo(comma_example_data)$mod_type == "6mA")
    result <- methylomeSummary(comma_example_data, mod_type = "6mA")
    expect_equal(result$n_sites[1], n_6mA)
})

test_that("methylomeSummary: mean_beta in [0,1]", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    mb <- result$mean_beta[!is.na(result$mean_beta)]
    expect_true(all(mb >= 0 & mb <= 1))
})

test_that("methylomeSummary: frac_methylated in [0,1]", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    fm <- result$frac_methylated[!is.na(result$frac_methylated)]
    expect_true(all(fm >= 0 & fm <= 1))
})

test_that("methylomeSummary: mean_beta matches manual computation", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    samp <- sampleInfo(comma_example_data)$sample_name[1]
    betas <- methylation(comma_example_data)[, samp]
    expected_mean <- mean(betas, na.rm = TRUE)
    actual_mean   <- result$mean_beta[result$sample_name == samp]
    expect_equal(actual_mean, expected_mean, tolerance = 1e-10)
})

test_that("methylomeSummary: n_covered matches non-NA count", {
    data(comma_example_data)
    result <- methylomeSummary(comma_example_data)
    samp <- sampleInfo(comma_example_data)$sample_name[1]
    betas <- methylation(comma_example_data)[, samp]
    expected_covered <- sum(!is.na(betas))
    actual_covered   <- result$n_covered[result$sample_name == samp]
    expect_equal(actual_covered, expected_covered)
})

test_that("methylomeSummary: error on non-commaData input", {
    expect_error(methylomeSummary(data.frame(x = 1)), "'object' must be a commaData")
})

test_that("methylomeSummary: error on unknown mod_type", {
    data(comma_example_data)
    expect_error(methylomeSummary(comma_example_data, mod_type = "9mX"), "not found in object")
})

test_that("methylomeSummary: all-NA sample gives NA mean_beta and median_beta, not error", {
    data(comma_example_data)
    methyl <- methylation(comma_example_data)
    methyl[, "ctrl_1"] <- NA_real_
    se_mod <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = methyl,
                       coverage    = coverage(comma_example_data)),
        rowData = SummarizedExperiment::rowData(comma_example_data),
        colData = SummarizedExperiment::colData(comma_example_data)
    )
    obj_mod <- new("commaData", se_mod,
                   genomeInfo = comma_example_data@genomeInfo,
                   annotation = comma_example_data@annotation,
                   motifSites = comma_example_data@motifSites)
    result    <- methylomeSummary(obj_mod)
    ctrl1_row <- result[result$sample_name == "ctrl_1", ]
    expect_true(is.na(ctrl1_row$mean_beta))
    expect_true(is.na(ctrl1_row$median_beta))
})
