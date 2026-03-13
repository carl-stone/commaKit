# calculateMethylSiteDepth is an internal (non-exported) function.
# Access it via CoMMA:::calculateMethylSiteDepth.

test_that("returns a dataframe with 'position' and 'coverage' columns", {
  result <- CoMMA:::calculateMethylSiteDepth(fixture_cov_df, position_col = "Position",
                                             cov_col = "coverage", w_size = 10L)
  expect_true(is.data.frame(result))
  expect_true("position" %in% names(result))
  expect_true("coverage" %in% names(result))
})

test_that("number of rows equals ceiling(max_position / w_size)", {
  # fixture_cov_df has positions 1:100; w_size=10 → ceiling(100/10) = 10 rows
  result <- CoMMA:::calculateMethylSiteDepth(fixture_cov_df, position_col = "Position",
                                             cov_col = "coverage", w_size = 10L)
  expect_equal(nrow(result), 10L)
})

test_that("uniform coverage is returned correctly for each window", {
  # All positions have coverage 10; every window mean should be 10.
  result <- CoMMA:::calculateMethylSiteDepth(fixture_cov_df, position_col = "Position",
                                             cov_col = "coverage", w_size = 10L)
  expect_true(all(result$coverage == 10))
})

test_that("calc_log2 = TRUE adds a 'log2_coverage' column", {
  result <- CoMMA:::calculateMethylSiteDepth(fixture_cov_df, position_col = "Position",
                                             cov_col = "coverage", w_size = 10L,
                                             calc_log2 = TRUE)
  expect_true("log2_coverage" %in% names(result))
})

test_that("log2_coverage values are log2 of coverage values", {
  result <- CoMMA:::calculateMethylSiteDepth(fixture_cov_df, position_col = "Position",
                                             cov_col = "coverage", w_size = 10L,
                                             calc_log2 = TRUE)
  expect_equal(result$log2_coverage, log2(result$coverage))
})

test_that("calc_log2 = FALSE does not add 'log2_coverage' column", {
  result <- CoMMA:::calculateMethylSiteDepth(fixture_cov_df, position_col = "Position",
                                             cov_col = "coverage", w_size = 10L,
                                             calc_log2 = FALSE)
  expect_false("log2_coverage" %in% names(result))
})

test_that("window size equal to total range produces one output row", {
  result <- CoMMA:::calculateMethylSiteDepth(fixture_cov_df, position_col = "Position",
                                             cov_col = "coverage", w_size = 100L)
  expect_equal(nrow(result), 1L)
})

test_that("window size of 1 produces one row per position", {
  small_df <- data.frame(Position = 1:5, coverage = c(2L, 4L, 6L, 8L, 10L),
                         stringsAsFactors = FALSE)
  result <- CoMMA:::calculateMethylSiteDepth(small_df, position_col = "Position",
                                             cov_col = "coverage", w_size = 1L)
  expect_equal(nrow(result), 5L)
})
