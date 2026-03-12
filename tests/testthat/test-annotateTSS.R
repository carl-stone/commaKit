test_that("position near sense strand TSS gets correct relative position", {
  # tuA: Left=1000, Strand="+". Window size 200.
  # Position 1100 is within [1000-200, 1000+200] = [800, 1200]. ✓
  # RelPos = 1100 - 1000 = 100
  methyl <- data.frame(Position = 1100L, beta = 0.5, stringsAsFactors = FALSE)
  result <- annotateTSS(methyl, fixture_tu_df, location = "Position", size = 200, long = FALSE)
  expect_true("RelPos_+1" %in% names(result))
  expect_equal(result$`RelPos_+1`, 100L)
})

test_that("position near antisense strand TSS gets correct relative position", {
  # tuB: Right=3500, Strand="-". Window size 200.
  # Position 3400 is within [3500-200, 3500+200] = [3300, 3700]. ✓
  # RelPos = 3500 - 3400 = 100
  methyl <- data.frame(Position = 3400L, beta = 0.5, stringsAsFactors = FALSE)
  result <- annotateTSS(methyl, fixture_tu_df, location = "Position", size = 200, long = FALSE)
  expect_true("RelPos_-1" %in% names(result))
  expect_equal(result$`RelPos_-1`, 100L)
})

test_that("position outside all TSS windows gets NoTSS = 'X'", {
  # Position 10 is far from both tuA (Left=1000) and tuB (Right=3500)
  methyl <- data.frame(Position = 10L, beta = 0.5, stringsAsFactors = FALSE)
  result <- annotateTSS(methyl, fixture_tu_df, location = "Position", size = 200, long = FALSE)
  expect_equal(result$NoTSS, "X")
})

test_that("position exactly at TSS (Left boundary) is annotated as RelPos = 0", {
  # Position 1000 == Left of tuA. RelPos = 1000 - 1000 = 0.
  methyl <- data.frame(Position = 1000L, beta = 0.5, stringsAsFactors = FALSE)
  result <- annotateTSS(methyl, fixture_tu_df, location = "Position", size = 200, long = FALSE)
  expect_true("RelPos_+1" %in% names(result))
  expect_equal(result$`RelPos_+1`, 0L)
})

test_that("long = FALSE returns wide format with one row per input site", {
  result <- annotateTSS(fixture_methyl_df, fixture_tu_df, location = "Position",
                        size = 200, long = FALSE)
  expect_equal(nrow(result), nrow(fixture_methyl_df))
})

test_that("long = TRUE pivots RelPos columns and removes NA rows", {
  # Position 1100 is near tuA only — one non-NA RelPos row should come back.
  methyl <- data.frame(Position = 1100L, beta = 0.5, stringsAsFactors = FALSE)
  result <- annotateTSS(methyl, fixture_tu_df, location = "Position", size = 200, long = TRUE)
  expect_true("RelPos" %in% names(result))
  expect_true("TSS_strand" %in% names(result))
  expect_false(any(is.na(result$RelPos)))
})

test_that("non-Transcription-Units features in meta_df are ignored", {
  # Add a Gene row to fixture_tu_df — it should not affect TSS annotation
  meta_with_extra <- rbind(
    fixture_tu_df,
    data.frame(Type = "Gene", Site = "geneX", Left = 1050L, Right = 1200L,
               Strand = "+", stringsAsFactors = FALSE)
  )
  methyl <- data.frame(Position = 1100L, beta = 0.5, stringsAsFactors = FALSE)
  result_extra <- annotateTSS(methyl, meta_with_extra, location = "Position",
                              size = 200, long = FALSE)
  result_base  <- annotateTSS(methyl, fixture_tu_df,    location = "Position",
                              size = 200, long = FALSE)
  expect_equal(result_extra$`RelPos_+1`, result_base$`RelPos_+1`)
})

test_that("window size = 0 only annotates position exactly at TSS", {
  # Position 1000 == Left of tuA; should be annotated.
  at_tss <- data.frame(Position = 1000L, beta = 0.5, stringsAsFactors = FALSE)
  result_at <- annotateTSS(at_tss, fixture_tu_df, location = "Position", size = 0L, long = FALSE)
  expect_false(is.na(result_at$`RelPos_+1`))

  # Position 1001 is one base away; should not be annotated.
  near_tss <- data.frame(Position = 1001L, beta = 0.5, stringsAsFactors = FALSE)
  result_near <- annotateTSS(near_tss, fixture_tu_df, location = "Position", size = 0L, long = FALSE)
  expect_equal(result_near$NoTSS, "X")
})

test_that("input row count is preserved when long = FALSE", {
  result <- annotateTSS(fixture_methyl_df, fixture_tu_df, location = "Position",
                        size = 500, long = FALSE)
  expect_equal(nrow(result), nrow(fixture_methyl_df))
})
