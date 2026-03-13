# Helper: build a small test dataframe with known positions and methylation
make_rolling_df <- function(genome_size = 100L, positions = c(10L, 30L, 50L, 70L, 90L),
                            values = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  data.frame(Position = positions, beta = values, stringsAsFactors = FALSE)
}

test_that("method = 'exact' returns a dataframe with correct columns", {
  df <- make_rolling_df()
  result <- methylRollingMedian(df, position_col = "Position", methyl_col = "beta",
                                w_size = 10L, genome_size = 100L, method = "exact")
  expect_true(is.data.frame(result))
  expect_true("position"   %in% names(result))
  expect_true("med_methyl" %in% names(result))
})

test_that("method = 'exact' returns exactly genome_size rows", {
  df <- make_rolling_df()
  result <- methylRollingMedian(df, position_col = "Position", methyl_col = "beta",
                                w_size = 10L, genome_size = 100L, method = "exact")
  expect_equal(nrow(result), 100L)
})

test_that("method = 'exact' positions run from 1 to genome_size", {
  df <- make_rolling_df()
  result <- methylRollingMedian(df, position_col = "Position", methyl_col = "beta",
                                w_size = 10L, genome_size = 100L, method = "exact")
  expect_equal(result$position[1],   1L)
  expect_equal(result$position[100], 100L)
})

test_that("method = 'exact' med_methyl values are in [0, 1]", {
  df <- make_rolling_df()
  result <- methylRollingMedian(df, position_col = "Position", methyl_col = "beta",
                                w_size = 10L, genome_size = 100L, method = "exact")
  med <- result$med_methyl
  expect_true(all(is.na(med) | (med >= 0 & med <= 1)))
})

test_that("method = 'exact' median at dense position equals known value", {
  # Single site at position 50 with beta 0.6; with a small window it should dominate.
  df <- data.frame(Position = 50L, beta = 0.6, stringsAsFactors = FALSE)
  result <- methylRollingMedian(df, position_col = "Position", methyl_col = "beta",
                                w_size = 1L, genome_size = 100L, method = "exact")
  # With w_size=1 the window spans 2 positions; only position 50 has data so
  # the rolling median (na.rm=TRUE) at that row should equal 0.6.
  expect_equal(result[result$position == 50, "med_methyl"], 0.6)
})

test_that("non-default genome_size is respected and returns that many rows", {
  df <- data.frame(Position = c(100L, 300L, 500L), beta = c(0.2, 0.5, 0.8),
                   stringsAsFactors = FALSE)
  result <- methylRollingMedian(df, position_col = "Position", methyl_col = "beta",
                                w_size = 50L, genome_size = 600L, method = "exact")
  expect_equal(nrow(result), 600L)
})

test_that("circular wraparound copies beginning sites to the end of the array", {
  # Use a tiny genome so we can inspect boundary effects directly.
  # Positions 1 and 2 have methylation; with wraparound they should influence
  # positions near genome_size even if no direct data exists there.
  df <- data.frame(Position = c(1L, 2L), beta = c(0.9, 0.9), stringsAsFactors = FALSE)
  result <- methylRollingMedian(df, position_col = "Position", methyl_col = "beta",
                                w_size = 5L, genome_size = 20L, method = "exact")
  # With wraparound the last few positions should have non-NA medians influenced
  # by the beginning-of-genome values. At minimum, the final row should not be NA.
  expect_false(is.na(result$med_methyl[20]))
})

test_that("KNOWN ISSUE: method = 'fast' uses wrong output column name 'mean_methyl'", {
  # The 'fast' method returns 'mean_methyl' instead of 'med_methyl'.
  # This test documents the naming inconsistency so it is caught if ever fixed.
  df <- make_rolling_df()
  result <- methylRollingMedian(df, position_col = "Position", methyl_col = "beta",
                                w_size = 10L, genome_size = 100L, method = "fast")
  expect_true("mean_methyl" %in% names(result))
  expect_false("med_methyl" %in% names(result))
})

test_that("KNOWN ISSUE: hardcoded default genome_size of 4641652 is E. coli-specific", {
  # Calling the function without genome_size on non-E.-coli data silently uses
  # the wrong genome size. This test documents the default so any change is visible.
  expect_equal(
    formals(methylRollingMedian)$genome_size,
    4641652
  )
})
