# methylRollingMean is an internal (non-exported) function.
# Access it via CoMMA:::methylRollingMean.

make_mean_df <- function(positions = c(10L, 30L, 50L, 70L, 90L),
                         values    = c(0.1,  0.3,  0.5,  0.7,  0.9)) {
  data.frame(Position = positions, beta = values, stringsAsFactors = FALSE)
}

test_that("returns a dataframe with 'position' and 'mean_methyl' columns", {
  df <- make_mean_df()
  result <- CoMMA:::methylRollingMean(df, position_col = "Position", methyl_col = "beta",
                                     w_size = 20L, genome_size = 100L)
  expect_true(is.data.frame(result))
  expect_true("position"    %in% names(result))
  expect_true("mean_methyl" %in% names(result))
})

test_that("returns one row per input site", {
  df <- make_mean_df()
  result <- CoMMA:::methylRollingMean(df, position_col = "Position", methyl_col = "beta",
                                     w_size = 20L, genome_size = 100L)
  expect_equal(nrow(result), nrow(df))
})

test_that("mean_methyl values are in [0, 1] for methylation inputs in [0, 1]", {
  df <- make_mean_df()
  result <- CoMMA:::methylRollingMean(df, position_col = "Position", methyl_col = "beta",
                                     w_size = 20L, genome_size = 100L)
  expect_true(all(result$mean_methyl >= 0 & result$mean_methyl <= 1))
})

test_that("single-site input returns that site's methylation value as the mean", {
  df <- data.frame(Position = 50L, beta = 0.75, stringsAsFactors = FALSE)
  result <- CoMMA:::methylRollingMean(df, position_col = "Position", methyl_col = "beta",
                                     w_size = 10L, genome_size = 100L)
  expect_equal(nrow(result), 1L)
  expect_equal(result$mean_methyl, 0.75)
})

test_that("output positions match the input positions", {
  df <- make_mean_df()
  result <- CoMMA:::methylRollingMean(df, position_col = "Position", methyl_col = "beta",
                                     w_size = 20L, genome_size = 100L)
  expect_equal(sort(result$position), sort(df$Position))
})

test_that("non-default genome_size parameter is accepted without error", {
  df <- data.frame(Position = c(100L, 300L, 500L), beta = c(0.2, 0.5, 0.8),
                   stringsAsFactors = FALSE)
  expect_no_error(
    CoMMA:::methylRollingMean(df, position_col = "Position", methyl_col = "beta",
                              w_size = 50L, genome_size = 600L)
  )
})

test_that("verbose = FALSE suppresses output", {
  df <- make_mean_df()
  expect_silent(
    CoMMA:::methylRollingMean(df, position_col = "Position", methyl_col = "beta",
                              w_size = 20L, genome_size = 100L, verbose = FALSE)
  )
})

test_that("KNOWN ISSUE: hardcoded default genome_size of 4641652 is E. coli-specific", {
  expect_equal(
    formals(CoMMA:::methylRollingMean)$genome_size,
    4641652
  )
})
