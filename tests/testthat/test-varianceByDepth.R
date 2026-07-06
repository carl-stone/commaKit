test_that("varianceByDepth: computes known variances and counts by filter", {
  beta <- rbind(
    c(s1 = 0.10, s2 = 0.20),
    c(s1 = 0.30, s2 = 0.60),
    c(s1 = 0.50, s2 = 0.40),
    c(s1 = 0.90, s2 = 0.80),
    c(s1 = 0.20, s2 = 0.10),
    c(s1 = 0.60, s2 = 0.70)
  )
  depth <- rbind(
    c(s1 = 10L, s2 = 10L),
    c(s1 = 10L, s2 = 10L),
    c(s1 = 20L, s2 = 20L),
    c(s1 = 20L, s2 = 30L),
    c(s1 = 10L, s2 = 10L),
    c(s1 = 20L, s2 = 20L)
  )
  object <- .make_commaData_fixture(
    beta = beta,
    coverage = depth,
    mod_type = c("6mA", "6mA", "6mA", "6mA", "5mC", "5mC"),
    motif = c("GATC", "GATC", "GATC", "GATC", "CCWGG", "CCWGG")
  )

  unfiltered <- varianceByDepth(
    object,
    coverage_bins = c(10L, 20L, 30L)
  )
  row.names(unfiltered) <- NULL
  expect_equal(
    unfiltered,
    data.frame(
      coverage = rep(c(10L, 20L, 30L), 2L),
      sample_name = c(rep("s1", 3L), rep("s2", 3L)),
      variance = c(
        0.01,
        0.0433333333333333,
        NA_real_,
        0.07,
        0.045,
        NA_real_
      ),
      n_sites = c(3L, 3L, 0L, 3L, 2L, 1L),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-10
  )

  only_6ma <- varianceByDepth(
    object,
    coverage_bins = c(10L, 20L, 30L),
    mod_type = "6mA"
  )
  row.names(only_6ma) <- NULL
  expect_equal(
    only_6ma,
    data.frame(
      coverage = rep(c(10L, 20L, 30L), 2L),
      sample_name = c(rep("s1", 3L), rep("s2", 3L)),
      variance = c(
        0.02,
        0.08,
        NA_real_,
        0.08,
        NA_real_,
        NA_real_
      ),
      n_sites = c(2L, 2L, 0L, 2L, 1L, 1L),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-10
  )
})

test_that("varianceByDepth: returns a data.frame", {
  data(comma_example_data)
  result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
  expect_s3_class(result, "data.frame")
})

test_that("varianceByDepth: has required columns", {
  data(comma_example_data)
  result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
  required_cols <- c("coverage", "sample_name", "variance", "n_sites")
  expect_true(all(required_cols %in% colnames(result)))
})

test_that("varianceByDepth: all sample names present", {
  data(comma_example_data)
  result <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
  expect_setequal(
    unique(result$sample_name),
    sampleInfo(comma_example_data)$sample_name
  )
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
  expect_true(nrow(low_n) > 0L) # condition must actually be tested
  expect_true(all(is.na(low_n$variance))) # variance must be NA for those rows
})

test_that("varianceByDepth: mod_type filtering works", {
  data(comma_example_data)
  result_6ma <- varianceByDepth(
    comma_example_data,
    coverage_bins = 5:30,
    mod_type = "6mA"
  )
  result_5mc <- varianceByDepth(
    comma_example_data,
    coverage_bins = 5:30,
    mod_type = "5mC"
  )
  # n_sites sums differ between mod types
  expect_false(
    identical(
      sum(result_6ma$n_sites, na.rm = TRUE),
      sum(result_5mc$n_sites, na.rm = TRUE)
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
  expect_error(
    varianceByDepth(data.frame(x = 1)),
    "'object' must be a commaData"
  )
})

test_that("varianceByDepth: error on unknown mod_type", {
  data(comma_example_data)
  expect_error(
    varianceByDepth(comma_example_data, mod_type = "9mZ"),
    "not found in object"
  )
})
