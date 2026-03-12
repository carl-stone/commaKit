# varByCoverage is an internal (non-exported) function.
# Access it via CoMMA:::varByCoverage.

test_that("returns a dataframe with 'coverage' and 'variance' columns", {
  result <- CoMMA:::varByCoverage(fixture_var_df)
  expect_true(is.data.frame(result))
  expect_true("coverage" %in% names(result))
  expect_true("variance" %in% names(result))
})

test_that("returns one row per unique coverage level", {
  # fixture_var_df has coverage levels 5 and 10
  result <- CoMMA:::varByCoverage(fixture_var_df)
  expect_equal(nrow(result), 2L)
})

test_that("coverage column contains all unique coverage levels from input", {
  result <- CoMMA:::varByCoverage(fixture_var_df)
  expect_setequal(as.integer(result$coverage), c(5L, 10L))
})

test_that("single coverage level produces one-row output", {
  single_cov <- data.frame(
    Coverage_Sample       = c(3L, 3L, 3L),
    Percent_Methyl_Sample = c(0.4, 0.5, 0.6),
    Ancestor_Mean         = c(0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
  result <- CoMMA:::varByCoverage(single_cov)
  expect_equal(nrow(result), 1L)
})

# --- Known issue test ---

test_that("KNOWN ISSUE: variance calculation subtracts dataframes producing list, not numeric", {
  # The implementation does: var(temp["Percent_Methyl_Sample"] - temp["Ancestor_Mean"])
  # Subtracting two single-column data.frames produces a data.frame, not a
  # numeric vector. var() on a data.frame returns a 1x1 matrix, not a scalar.
  # This documents the bug so it is caught if fixed.
  result <- CoMMA:::varByCoverage(fixture_var_df)
  # The variance column will contain matrices or lists rather than plain numerics.
  # We test that the output is NOT a plain numeric vector to document the bug.
  # (If this test starts failing, the bug has been fixed — delete this test and
  # replace with a correctness check.)
  variance_val <- result[result$coverage == 5, "variance"]
  expect_false(
    is.numeric(variance_val) && length(variance_val) == 1,
    label = "variance column should NOT be a plain numeric scalar (known bug)"
  )
})
