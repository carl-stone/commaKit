test_that("position inside a single feature gets that feature's column added", {
  result <- annotateMethylSites(fixture_methyl_df, fixture_meta_df, location = "Position")
  # Position 500 falls inside Gene [400, 700]
  expect_equal(result[result$Position == 500, "Gene"], "geneA")
})

test_that("position inside no feature gets No_Feature = '1'", {
  result <- annotateMethylSites(fixture_methyl_df, fixture_meta_df, location = "Position")
  # Positions 1500 and 3000 are outside all features
  expect_equal(result[result$Position == 1500, "No_Feature"], "1")
  expect_equal(result[result$Position == 3000, "No_Feature"], "1")
})

test_that("position inside multiple overlapping features gets all feature columns", {
  # Position 425 falls inside both Gene [400, 700] and Promoter [200, 450]
  methyl_overlap <- data.frame(Position = 425L, beta = 0.5, stringsAsFactors = FALSE)
  result <- annotateMethylSites(methyl_overlap, fixture_meta_df, location = "Position")
  expect_equal(result[result$Position == 425, "Gene"],     "geneA")
  expect_equal(result[result$Position == 425, "Promoter"], "promA")
})

test_that("position exactly on the left boundary of a feature is annotated", {
  methyl_at_left <- data.frame(Position = 400L, beta = 0.5, stringsAsFactors = FALSE)
  result <- annotateMethylSites(methyl_at_left, fixture_meta_df, location = "Position")
  expect_equal(result[result$Position == 400, "Gene"], "geneA")
  expect_true(is.na(result[result$Position == 400, "No_Feature"]))
})

test_that("position exactly on the right boundary of a feature is annotated", {
  methyl_at_right <- data.frame(Position = 700L, beta = 0.5, stringsAsFactors = FALSE)
  result <- annotateMethylSites(methyl_at_right, fixture_meta_df, location = "Position")
  expect_equal(result[result$Position == 700, "Gene"], "geneA")
})

test_that("annotated and unannotated positions coexist correctly in one call", {
  result <- annotateMethylSites(fixture_methyl_df, fixture_meta_df, location = "Position")
  # Position 500: annotated — Gene column set, No_Feature absent
  expect_false(is.na(result[result$Position == 500, "Gene"]))
  expect_true(is.na(result[result$Position == 500, "No_Feature"]))
  # Position 1500: unannotated — No_Feature set, Gene absent
  expect_equal(result[result$Position == 1500, "No_Feature"], "1")
})

test_that("input dataframe row count is preserved", {
  result <- annotateMethylSites(fixture_methyl_df, fixture_meta_df, location = "Position")
  expect_equal(nrow(result), nrow(fixture_methyl_df))
})

test_that("non-Position column name is accepted via 'location' argument", {
  methyl_renamed <- fixture_methyl_df
  names(methyl_renamed)[names(methyl_renamed) == "Position"] <- "Pos"
  result <- annotateMethylSites(methyl_renamed, fixture_meta_df, location = "Pos")
  expect_equal(nrow(result), nrow(methyl_renamed))
  expect_equal(result[result$Pos == 500, "Gene"], "geneA")
})

test_that("single-row methyl_df with overlap is annotated correctly", {
  single <- data.frame(Position = 600L, beta = 0.9, stringsAsFactors = FALSE)
  result <- annotateMethylSites(single, fixture_meta_df, location = "Position")
  expect_equal(nrow(result), 1L)
  expect_equal(result$Gene, "geneA")
})

test_that("single-row methyl_df with no overlap gets No_Feature", {
  single <- data.frame(Position = 9999L, beta = 0.1, stringsAsFactors = FALSE)
  result <- annotateMethylSites(single, fixture_meta_df, location = "Position")
  expect_equal(nrow(result), 1L)
  expect_equal(result$No_Feature, "1")
})

# --- Known issue tests: document current behavior so regressions are caught ---

test_that("KNOWN ISSUE: empty methyl_df does not error", {
  # The for-loop silently iterates zero times on an empty vector.
  # Document that the function survives and returns a zero-row frame.
  empty <- data.frame(Position = integer(0), beta = double(0), stringsAsFactors = FALSE)
  expect_no_error(
    result <- annotateMethylSites(empty, fixture_meta_df, location = "Position")
  )
  expect_equal(nrow(result), 0L)
})
